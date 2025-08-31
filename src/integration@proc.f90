! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Quadrature schemes and integration-related tools.
submodule (integration) proc
  use param, only: mmlen
  implicit none

  !xx! private procedures
  ! subroutine int_output_header(bas,res)
  ! subroutine int_output_fields(bas,res)
  ! subroutine int_reorder_gridout(ff,bas)
  ! subroutine intgrid_fields(bas,res)
  ! subroutine intgrid_hirshfeld_fields(bas,res)
  ! subroutine intgrid_hirshfeld_overlap(bas,res)
  ! subroutine intgrid_deloc(bas,res)
  ! subroutine write_sijchk(sijfname,nbnd,nbndw,nlat,nmo,nlattot,nspin,nattr,sijtype,kpt,occ,sij)
  ! subroutine write_sijrestart(sijfname,nbnd,nbndw,nwan,nmo,nspin,nattr,sijtype,is,ibnd1,ibnd2,sij)
  ! function read_sijrestart(sijfname,nbnd,nbndw,nwan,nmo,nspin,nattr,sijtype,is,ibnd1,ibnd2,sij) result(valid)
  ! subroutine write_fachk(fafname,nbnd,nbndw,nlat,nmo,nlattot,nspin,nattr,sijtype,fa)
  ! function read_chk_header(fname,nbnd,nbndw,nlat,nmo,nlattot,nspin,nattr,sijtype) result(haschk)
  ! subroutine read_sijchk_body(sijfname,kpt,occ,sij)
  ! subroutine read_fachk_body(fafname,fa)
  ! subroutine calc_sij_wannier(fid,wancut,useu,restart,imtype,natt1,iatt,ilvec,idg1,xattr,dat,luevc,luevc_ibnd,sij)
  ! subroutine calc_sij_psink(fid,imtype,natt1,iatt,ilvec,idg1,xattr,dat,sij)
  ! subroutine calc_fa_wannier(res,nmo,nbnd,nlat,nattr,nspin,di3)
  ! subroutine calc_di3_wannier(res,nmo,nbnd,nlat,nattr,nspin,atom1,atom2)
  ! subroutine calc_fa_psink(res,nmo,nbnd,nlat,nattr,nspin,kpt,occ)
  ! subroutine calc_di3_psink(res,nmo,nbnd,nlat,nattr,nspin,kpt,occ,atom1,atom2)
  ! subroutine find_sij_translations(res,nmo,nbnd,nlat,nlattot)
  ! subroutine check_sij_sanity(res,nspin,nmo,nbnd,nlat,nlattot)
  ! function quadpack_f(x,unit,xnuc) result(res)
  ! subroutine int_output_multipoles(bas,res)
  ! subroutine int_output_deloc(bas,res)
  ! subroutine int_output_hirshfeld_overlap(bas,res)
  ! subroutine int_output_json(jsonfile,bas,res)
  ! subroutine assign_strings(i,icp,usesym,scp,sncp,sname,smult,sz)
  ! subroutine int_gridbasins(bas)
  ! subroutine int_cubew(bas)
  ! subroutine unpackidx(idx,io,jo,ko,bo,nmo,nbnd,nlat)
  ! subroutine packidx(io,jo,ko,bo,idx,nmo,nbnd,nlat)

contains

  !> Driver for the integration in grids
  module subroutine intgrid_driver(line)
    use hirshfeld, only: hirsh_grid, voronoi_grid
    use bader, only: bader_integrate
    use yt, only: yt_integrate, yt_isosurface, yt_weights, ytdata, ytdata_clean
    use systemmod, only: sy
    use fieldmod, only: type_grid
    use grid3mod, only: grid3
    use global, only: eval_next, dunit0, iunit, iunitname0
    use tools_io, only: ferror, faterr, lgetword, equal, isexpression_or_word, uout,&
       string, fclose, isinteger, getword
    use types, only: basindat, int_result

    character*(*), intent(in) :: line

    real*8, parameter :: ratom_def0 = 1d0

    character(len=:), allocatable :: word, jsonfile
    integer :: n(3), ntot, idum
    integer :: lp, lp2
    logical :: ok, nonnm, noatoms
    real*8 :: ratom_def
    real*8, allocatable :: faux(:,:,:)
    type(int_result), allocatable :: res(:)
    type(basindat) :: bas

    ! only grids
    if (.not.sy%isinit) then
       call ferror("intgrid_integrate","system not initialized",faterr)
       return
    end if
    if (.not.allocated(sy%c)) then
       call ferror("intgrid_integrate","system does not have crystal",faterr)
       return
    end if
    if (.not.sy%goodfield(sy%iref)) then
       call ferror("intgrid_integrate","reference field not initalized",faterr)
       return
    end if
    if (sy%f(sy%iref)%type /= type_grid) then
       call ferror("intgrid_driver","BADER/YT/ISOSURFACE/HIRSHFELD/VORONOI can only be used with grids",faterr,line,syntax=.true.)
       return
    end if
    if (sy%npropi <= 0) then
       call ferror("intgrid_driver","no integrable properties",faterr,line,syntax=.true.)
       return
    end if

    ! method and header
    lp = 1
    word = lgetword(line,lp)
    bas%higher = .true.
    allocate(bas%docelatom(0:sy%f(sy%iref)%ncpcel))
    bas%docelatom = .true.
    if (equal(word,"yt")) then
       bas%imtype = imtype_yt
    elseif (equal(word,"bader")) then
       bas%imtype = imtype_bader
    elseif (equal(word,"isosurface")) then
       bas%imtype = imtype_isosurface
       word = lgetword(line,lp)
       if (equal(word,'higher')) then
          bas%higher = .true.
       elseif (equal(word,'lower')) then
          bas%higher = .false.
       else
          call ferror("intgrid_driver","ISOSURFACE must be followed by HIGHER or LOWER",faterr,line,syntax=.true.)
          return
       end if
       ok = eval_next(bas%isov,line,lp)
       if (.not.ok) then
          call ferror("intgrid_driver","ISOSURFACE HIGHER/LOWER must be followed by the contour value",faterr,line,syntax=.true.)
          return
       end if
    elseif (equal(word,"hirshfeld")) then
       bas%imtype = imtype_hirshfeld
    elseif (equal(word,"voronoi")) then
       bas%imtype = imtype_voronoi
    else
       call ferror("intgrid_driver","wrong method",faterr,line,syntax=.true.)
       return
    endif

    ! parse the input
    jsonfile = ""
    ratom_def = ratom_def0
    nonnm = .true.
    noatoms = .false.
    bas%wcube = .false.
    bas%ndrawbasin = -1
    bas%basinfmt = "obj"
    bas%expr = ""
    do while(.true.)
       word = lgetword(line,lp)
       if (equal(word,"nnm") .and. (bas%imtype == imtype_bader .or. bas%imtype == imtype_yt)) then
          nonnm = .false.
       elseif (equal(word,"noatoms") .and. (bas%imtype == imtype_bader .or. bas%imtype == imtype_yt)) then
          noatoms = .true.
       elseif (equal(word,"ratom") .and. (bas%imtype == imtype_bader .or. bas%imtype == imtype_yt)) then
          nonnm = .false.
          ok = eval_next(ratom_def,line,lp)
          if (.not.ok) then
             call ferror("intgrid_driver","wrong RATOM keyword",faterr,line,syntax=.true.)
             return
          end if
          ratom_def = ratom_def / dunit0(iunit)
       elseif (equal(word,"wcube") .and. bas%imtype /= imtype_voronoi) then
          bas%wcube = .true.
       elseif (equal(word,"basins") .and. bas%imtype /= imtype_hirshfeld) then
          lp2 = lp
          word = lgetword(line,lp)
          if (equal(word,"obj")) then
             bas%basinfmt = "obj"
          elseif (equal(word,"ply")) then
             bas%basinfmt = "ply"
          elseif (equal(word,"off")) then
             bas%basinfmt = "off"
          else
             lp = lp2
             bas%basinfmt = "obj"
          end if
          ok = isinteger(bas%ndrawbasin,line,lp)
          if (.not.ok) &
             bas%ndrawbasin = 0
       elseif (equal(word,"discard") .and.&
          (bas%imtype == imtype_bader .or. bas%imtype == imtype_yt .or. bas%imtype == imtype_isosurface)) then
          ok = isexpression_or_word(bas%expr,line,lp)
          if (.not. ok) then
             call ferror("intgrid_driver","wrong DISCARD keyword",faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"only") .and. (bas%imtype /= imtype_isosurface)) then
          bas%docelatom = .false.
          do while (.true.)
             lp2 = lp
             ok = isinteger(idum,line,lp)
             if (.not.ok) then
                lp = lp2
                exit
             end if
             if (idum < 0 .or. idum > sy%f(sy%iref)%ncpcel) then
                call ferror("intgrid_driver","Invalid complete CP list nuclear ID",faterr,line,syntax=.true.)
                return
             end if
             bas%docelatom(idum) = .true.
          end do
       elseif (equal(word,"json") .and. (bas%imtype == imtype_bader .or. bas%imtype == imtype_yt)) then
          jsonfile = getword(line,lp)
          if (len_trim(jsonfile) == 0) then
             call ferror("intgrid_driver","Invalid JSON file",faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror("intgrid_driver","Unknown extra keyword: " // word,faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! field number of grid points
    bas%n = sy%f(sy%iref)%grid%n
    ntot = n(1)*n(2)*n(3)

    ! distance for atom assignment
    if (noatoms) then
       bas%atexist = .false.
       bas%ratom = ratom_def
    elseif (nonnm) then
       bas%atexist = .true.
       bas%ratom = 1d40
    else
       bas%atexist = .true.
       bas%ratom = ratom_def
    end if

    ! prepare the array for the basin field
    allocate(bas%f(bas%n(1),bas%n(2),bas%n(3)))
    if (sy%f(sy%iref)%usecore .and. (bas%imtype == imtype_yt .or. bas%imtype == imtype_bader)) then
       ! yt,bader with usecore -> augment with core
       call sy%c%promolecular_array3(faux,sy%f(sy%iref)%grid%n,sy%f(sy%iref)%zpsp)
       bas%f = sy%f(sy%iref)%grid%f + faux
    elseif (bas%imtype == imtype_isosurface .and..not.bas%higher) then
       ! isosurface lower -> minus the reference field
       bas%f = -sy%f(sy%iref)%grid%f
       bas%isov = -bas%isov
    elseif (bas%imtype == imtype_hirshfeld) then
       ! hirshfeld -> the promolecular density
       call sy%c%promolecular_array3(faux,sy%f(sy%iref)%grid%n)
       bas%f = faux
    else
       ! rest -> reference field
       bas%f = sy%f(sy%iref)%grid%f
    end if

    ! call the integration method
    bas%luw = 0
    if (bas%imtype == imtype_bader) then
       write (uout,'("* Henkelman et al. integration ")')
       write (uout,'("  Please cite: ")')
       write (uout,'("    G. Henkelman, A. Arnaldsson, and H. Jonsson,")')
       write (uout,'("       Comput. Mater. Sci. 36, 254-360 (2006). (10.1016/j.commatsci.2005.04.010)")')
       write (uout,'("    E. Sanville, S. Kenny, R. Smith, and G. Henkelman,")')
       write (uout,'("       J. Comput. Chem. 28, 899-908 (2007). (10.1002/jcc.20575)")')
       write (uout,'("    W. Tang, E. Sanville, and G. Henkelman,")')
       write (uout,'("       J. Phys.: Condens. Matter 21, 084204 (2009). (10.1088/0953-8984/21/8/084204)")')
       write (uout,'("+ Distance for atomic assignment (",A,"): ",A)') iunitname0(iunit),&
          string(max(bas%ratom,0d0),'e',decimal=4)
       if (len_trim(bas%expr) > 0) &
          write (uout,'("+ Discard attractor expression: ",A)') trim(bas%expr)
       call bader_integrate(sy,bas,sy%iref)
       write (uout,'("+ Attractors in BADER: ",A)') string(bas%nattr)
    elseif (bas%imtype == imtype_yt) then
       write (uout,'("* Yu-Trinkle integration ")')
       write (uout,'("  Please cite: ")')
       write (uout,'("  Min Yu, Dallas Trinkle, J. Chem. Phys. 134 (2011) 064111. (10.1063/1.3553716)")')
       write (uout,'("+ Distance for atomic assignment (",A,"): ",A)') iunitname0(iunit),&
          string(max(bas%ratom,0d0),'e',decimal=4)
       if (len_trim(bas%expr) > 0) &
          write (uout,'("+ Discard attractor expression: ",A)') trim(bas%expr)
       call yt_integrate(sy,bas)
       write (uout,'("+ Attractors in YT: ",A)') string(bas%nattr)
    elseif (bas%imtype == imtype_isosurface) then
       write (uout,'("* Isosurface integration ")')
       if (bas%higher) then
          write (uout,'("+ Integrated regions with reference field HIGHER than ",A)') string(bas%isov,'e',8,4)
       else
          write (uout,'("+ Integrated regions with reference field LOWER than ",A)') string(-bas%isov,'e',8,4)
       end if
       if (len_trim(bas%expr) > 0) &
          write (uout,'("+ Discard expression: ",A)') trim(bas%expr)
       call yt_isosurface(sy,bas)
    elseif (bas%imtype == imtype_hirshfeld) then
       write (uout,'("* Hirshfeld integration ")')
       call hirsh_grid(sy,bas)
    elseif (bas%imtype == imtype_voronoi) then
       write (uout,'("* Voronoi integration ")')
       call voronoi_grid(sy,bas)
    endif

    ! Reorder the attractors
    call int_reorder_gridout(sy%f(sy%iref),bas)
    if (bas%imtype == imtype_yt .or. bas%imtype == imtype_bader) then
       write (uout,'("+ Attractors after reordering: ",A)') string(bas%nattr)
       write (uout,*)
    end if

    ! Write weight cubes
    call int_cubew(bas)

    ! Bains plotting
    call int_gridbasins(bas)

    ! Prepare for the calculation of properties
    allocate(res(sy%npropi))
    write (uout,'("+ Integrating atomic properties"/)')

    ! Integrate scalar fields and multipoles
    if (bas%imtype == imtype_hirshfeld) then
       call intgrid_hirshfeld_fields(bas,res)
    else
       call intgrid_fields(bas,res)
    end if

    ! localization and delocalization indices
    if (bas%imtype == imtype_bader .or. bas%imtype == imtype_yt) then
       call intgrid_deloc(bas,res)
    elseif (bas%imtype == imtype_hirshfeld) then
       call intgrid_hirshfeld_overlap(bas,res)
    end if

    ! deallocate the basin field
    deallocate(bas%f)

    ! header for integration output
    write (uout,*)
    call int_output_header(bas,res)

    ! output integrated scalar field properties
    call int_output_fields(bas,res)

    ! output multipoles
    call int_output_multipoles(bas,res)

    ! output localization and delocalization indices
    call int_output_deloc(bas,res)

    ! output hirshfeld overlap populations
    call int_output_hirshfeld_overlap(bas,res)

    ! write the json file
    if (len(jsonfile) > 0) then
       call int_output_json(jsonfile,bas,res)
    end if

    ! clean up
    if (bas%imtype == imtype_yt .and. bas%luw /= 0) then
       call fclose(bas%luw)
    endif
    deallocate(res)

  end subroutine intgrid_driver

  !> Do a radial numerical quadrature on the given ray, between the
  !> selected radii and return the properties. The r^2 factor is
  !> included.
  module subroutine int_radialquad(x,theta,phi,r0,rend,lprop,abserr,neval,iaserr,ierr)
    use systemmod, only: sy
    use quadpack, only: dqags, dqng, dqag
    use global, only: int_radquad_type, int_gauleg, int_radquad_nr, int_qags, &
       int_radquad_relerr, int_qng, int_radquad_abserr, int_qag, int_iasprec
    use tools_math, only: gauleg
    real*8, intent(in) :: x(3)  !< The center of the basin (cartesian coords)
    real*8, intent(in) :: theta !< Polar angle of the ray
    real*8, intent(in) :: phi   !< Azimuthal angle or the ray
    real*8, intent(in) :: r0    !< Left limit of the radial interval
    real*8, intent(in) :: rend  !< Right limit of the radial interval
    real*8, intent(out) :: lprop(sy%npropi) !< The integrated properties
    real*8, intent(out) :: abserr !< Estimated absolute error
    integer, intent(out) :: neval !< Number of evaluations of grdall
    real*8, intent(out) :: iaserr(sy%npropi) !< Estimated IAS precision error
    integer, intent(out) :: ierr

    integer :: k, isign
    real*8 :: unit(3), xaux(3), r1, r2
    real*8, dimension(sy%npropi) :: atprop
    integer :: ier
    ! quadpack internals
    integer, parameter :: limit = 10000
    integer, parameter :: lenw = 4*limit
    integer :: iwork(limit)
    real*8 :: work(lenw)
    integer :: last, key
    ! radial integration
    real*8, allocatable :: rpoints(:), rweights(:)

    unit = (/ sin(theta) * cos(phi),&
              sin(theta) * sin(phi),&
              cos(theta) /)
    ierr = 0

    if (rend > r0) then
       r1 = r0
       r2 = rend
       isign = 1
    else
       r1 = rend
       r2 = r0
       isign = -1
    end if

    select case(INT_radquad_type)
    case (INT_gauleg)
       ! allocate and set radial nodes
       allocate(rpoints(INT_radquad_nr),rweights(INT_radquad_nr))
       call gauleg(r1,r2,rpoints,rweights,INT_radquad_nr)

       ! do the quadrature
       lprop = 0d0
       do k = 1, INT_radquad_nr
          xaux = x + rpoints(k) * unit
          call sy%grdall(xaux,atprop)
          lprop = lprop + rpoints(k)**2 * rweights(k) * atprop
       end do
       lprop = lprop * isign
       abserr = 0d0
       neval = INT_radquad_nr

       ! deallocate
       deallocate(rpoints,rweights)

    case (INT_qags)
       call dqags(quadpack_f,unit,x,r0,rend,INT_radquad_abserr,&
          INT_radquad_relerr,lprop,abserr,neval,&
          ier,limit,lenw,last,iwork,work)
       ierr = max(ierr,ier)

    case (INT_qng)
       call dqng(quadpack_f,unit,x,r0,rend,INT_radquad_abserr,&
          INT_radquad_relerr,lprop,abserr,neval,ier)
       ierr = max(ierr,ier)

    case (INT_qag)
       if (INT_radquad_nr >= 30) then
          key = 6
       else if (INT_radquad_nr >= 25) then
          key = 5
       else if (INT_radquad_nr >= 20) then
          key = 4
       else if (INT_radquad_nr >= 15) then
          key = 3
       else if (INT_radquad_nr >= 10) then
          key = 2
       else
          key = 1
       end if
       call dqag(quadpack_f,unit,x,r0,rend,INT_radquad_abserr,&
          INT_radquad_relerr,key,lprop,abserr,neval,&
          ier,limit,lenw,last,iwork,work)
       ierr = max(ierr,ier)
    end select

    ! IAS precision error
    xaux = x + rend * unit
    call sy%grdall(xaux,atprop)
    iaserr = abs(atprop * rend**2 * INT_iasprec)
    neval = neval + 1

  end subroutine int_radialquad

  !> Sum the gauss-legendre 2d (ntheta * nphi) quadrature. The srf
  !> minisurface contains the rays generated by gauleg_nodes, with
  !> radii corresponding to the IAS. ntheta and nphi. The resulting
  !> integrated properties are returned in the lprop vector. abserr is
  !> the integrated error of the radial quadrature. iaserr is the error
  !> caused by IAS inaccuracies.  neval is the number of evaluations of
  !> grdall used. rbeta is a reference radius (beta-sphere).
  module subroutine gauleg_mquad(srf,ntheta,nphi,rbeta,lprop,abserr,neval,iaserr)
    use systemmod, only: sy
    use surface, only: minisurf
    use tools_math, only: gauleg
    use tools_io, only: ferror, warning, uout
    use param, only: pi
    type(minisurf), intent(inout) :: srf !< Surface representing the basin
    integer, intent(in) :: ntheta !< Number of polar points
    integer, intent(in) :: nphi   !< Number of azimuthal points
    real*8, intent(in) :: rbeta   !< Beta-spehre radius
    real*8, intent(out) :: lprop(sy%npropi) !< Properties vector
    real*8, intent(out) :: abserr !< Integrated radial quad. error
    integer, intent(out) :: neval !< Number of evaluations
    real*8, intent(out) :: iaserr(sy%npropi) !< Integrated IAS precision error

    integer :: i, j, v, leval, nneval
    real*8 :: rprop(sy%npropi), pprop(sy%npropi)
    real*8 :: rerr, perr, riaserr(sy%npropi), piaserr(sy%npropi)
    integer :: realnphi, ierr, err, vidx(ntheta,nphi+1), lerr
    ! Nodes and weights
    real*8, allocatable, dimension(:) :: tpoints
    real*8, allocatable, dimension(:) :: tweights
    real*8, allocatable, dimension(:) :: ppoints
    real*8, allocatable, dimension(:) :: pweights

    allocate(tpoints(ntheta),tweights(ntheta))
    allocate(ppoints(nphi+1),pweights(nphi+1))
    ppoints = 0d0
    pweights = 0d0

    ! initialize
    neval = 0
    lprop = 0d0
    abserr = 0d0
    iaserr = 0d0
    err = 0

    ! set theta-nodes
    call gauleg(0d0,pi,tpoints,tweights,ntheta)

    ! fill idx array
    v = 0
    do i = 1, ntheta
       do j = 1, int(nphi*abs(sin(tpoints(i))))+1
          v = v + 1
          vidx(i,j) = v
       end do
    end do

    !$omp parallel do &
    !$omp private(realnphi,nneval,pprop,perr,piaserr,lerr,v,rprop,rerr,leval,riaserr,ierr) &
    !$omp firstprivate(ppoints,pweights)
    do i = 1, ntheta
       realnphi = int(nphi*abs(sin(tpoints(i))))+1
       call gauleg(0d0,2d0*pi,ppoints,pweights,realnphi)
       nneval = 0
       pprop = 0d0
       perr = 0d0
       piaserr = 0d0
       lerr = 0
       do j = 1, realnphi
          v = vidx(i,j)
          call int_radialquad(srf%n,srf%th(v),srf%ph(v),rbeta,srf%r(v),&
             rprop,rerr,leval,riaserr,ierr)
          pprop = pprop + rprop * pweights(j)
          perr = perr + rerr * pweights(j)
          piaserr = piaserr + riaserr * pweights(j)
          nneval = nneval + leval
          lerr = max(lerr,ierr)
       end do
       !$omp critical (sum_gauleg_mquad)
       lprop = lprop + sin(srf%th(v)) * tweights(i) * pprop
       iaserr = iaserr + sin(srf%th(v)) * tweights(i) * piaserr
       abserr = abserr + sin(srf%th(v)) * tweights(i) * perr
       neval = neval + nneval
       err = max(err,lerr)
       !$omp end critical (sum_gauleg_mquad)
    end do
    !$omp end parallel do
    deallocate(tpoints,tweights,ppoints,pweights)

    if (err /= 0) then
       call ferror('gauleg_mquad','Radial integration had non-zero error code',warning)
       write (uout,'(a,I2)') " ier = ", err
       write (uout,'(a)') " Check the routine documentation for more info."
    end if

  end subroutine gauleg_mquad

  !> Sum the Lebedev 2d quadrature with npts points. The srf
  !> minisurface contains the rays generated by lebedev_nodes, with
  !> radii corresponding to the IAS. ntheta and nphi. The resulting
  !> integrated properties are returned in the lprop vector. abserr is
  !> the integrated error of the radial quadrature. iaserr is the error
  !> caused by IAS inaccuracies.  neval is the number of evaluations of
  !> grdall used. rbeta is a reference radius (beta-sphere).
  module subroutine lebedev_mquad(srf,npts,rbeta,lprop,abserr,neval,iaserr)
    use systemmod, only: sy
    use surface, only: minisurf
    use tools_math, only: select_lebedev
    use tools_io, only: uout, ferror, warning
    type(minisurf), intent(inout) :: srf !< Surface representing the basin
    integer, intent(in) :: npts   !< Number of points
    real*8, intent(in) :: rbeta   !< Beta-spehre radius
    real*8, intent(out) :: lprop(sy%npropi) !< Properties vector
    real*8, intent(out) :: abserr !< Integrated radial quad. error
    integer, intent(out) :: neval !< Number of evaluations
    real*8, intent(out) :: iaserr(sy%npropi) !< Integrated IAS precision error

    integer :: i, leval
    real*8 :: rprop(sy%npropi)
    real*8 :: rerr, riaserr(sy%npropi)
    integer :: ierr, err
    ! Nodes and weights
    real*8 :: xleb(5810)
    real*8 :: yleb(5810)
    real*8 :: zleb(5810)
    real*8 :: wleb(5810)

    call select_lebedev(npts,xleb,yleb,zleb,wleb)

    neval = 0
    lprop = 0d0
    abserr = 0d0
    iaserr = 0d0
    err = 0

    !$omp parallel do private(rprop,rerr,leval,riaserr,ierr)
    do i = 1, npts
       call int_radialquad(srf%n,srf%th(i),srf%ph(i),rbeta,srf%r(i),&
          rprop,rerr,leval,riaserr,ierr)
       !$omp critical (sum_lebedev_mquad)
       lprop = lprop + rprop * wleb(i)
       iaserr = iaserr + riaserr * wleb(i)
       abserr = abserr + rerr  * wleb(i)
       neval = neval + leval
       err = max(err,ierr)
       !$omp end critical (sum_lebedev_mquad)
    end do
    !$omp end parallel do

    if (err /= 0) then
       call ferror('lebedev_mquad','Radial integration had non-zero error code',warning)
       write (uout,'(a,I2)') " ier = ", err
       write (uout,'(a)') " Check the routine documentation for more info."
    end if

  end subroutine lebedev_mquad

  !xx! private procedures

  !> Write the header for the integration results section to the output. Includes
  !> a list of properties integrated in the basins. bas = integration driver
  !> data, res(1:npropi) = results. If nomol0, prevent the output of molecular
  !> properties. If usesym0, use multiplicities.
  module subroutine int_output_header(bas,res,nomol0,usesym0)
    use systemmod, only: sy, itype_v, itype_expr, itype_mpoles, itype_names,&
       itype_deloc_psink, itype_deloc_wnr
    use global, only: iunitname0, iunit, dunit0
    use tools_io, only: uout, string, ioj_left, ioj_center, ioj_right
    use types, only: basindat, int_result
    type(basindat), intent(in) :: bas
    type(int_result), intent(in) :: res(:)
    logical, intent(in), optional :: nomol0
    logical, intent(in), optional :: usesym0

    integer :: i, j, k, fid, id(3)
    character(len=:), allocatable :: saux, label, cini, itaux
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    real*8 :: x(3), xcm(3)
    integer, allocatable :: idxmol(:,:)
    logical :: nomol, usesym

    nomol = .false.
    if (present(nomol0)) nomol = nomol0
    usesym = .false.
    if (present(usesym0)) usesym = usesym0

    ! List of integrable properties accepted/rejected and why some
    ! were rejected
    if (bas%imtype == imtype_isosurface) then
       write (uout,'("* List of properties integrated in the isosurface regions")')
    elseif (bas%imtype == imtype_hirshfeld .or. bas%imtype == imtype_voronoi) then
       write (uout,'("* List of integrated atomic properties")')
    else
       write (uout,'("* List of properties integrated in the attractor basins")')
    end if
    write (uout,'("+ The ""Label"" entries will be used to identify the integrable property")')
    write (uout,'("  in the tables of integrated atomic properties. The entries with")')
    write (uout,'("  an ""x"" will not be integrated.")')
    write (uout,'("# id    Label      fid  Field       Additional")')
    do i = 1, sy%npropi
       fid = sy%propi(i)%fid
       if (.not.res(i)%done) then
          label = "--inactive--"
          cini = "x "
          if (len_trim(res(i)%reason) > 0) then
             saux = "Reason: " // string(res(i)%reason)
          else
             saux = ""
          end if
       else
          label = sy%propi(i)%prop_name
          cini = "  "
       end if
       if (sy%propi(i)%itype == itype_v) then
          if (res(i)%done) saux = ""
          itaux = "--"
       elseif (sy%propi(i)%itype == itype_expr) then
          if (res(i)%done) saux = "expr = " // string(sy%propi(i)%expr)
          itaux = "--"
       elseif (sy%propi(i)%itype == itype_mpoles) then
          if (res(i)%done) saux = "Lmax = " // string(sy%propi(i)%lmax)
          itaux = string(fid)
       elseif (sy%propi(i)%itype == itype_deloc_wnr .or. sy%propi(i)%itype == itype_deloc_psink) then
          if (res(i)%done .and. sy%propi(i)%di3) then
             if (sy%propi(i)%di3_atom1 < 0 .or. sy%propi(i)%di3_atom1 > bas%nattr) then
                saux = " 3-center indices calculated "
             elseif (sy%propi(i)%di3_atom2(1) < 0 .or. sy%propi(i)%di3_atom2(1) > bas%nattr) then
                saux = " 3-center indices calculated for attractor " // string(sy%propi(i)%di3_atom1) // " "
             else
                saux = " 3-center indices calculated for attractor pair " // string(sy%propi(i)%di3_atom1) //&
                   "," // string(sy%propi(i)%di3_atom2(1)) // "+(" //&
                   string(sy%propi(i)%di3_atom2(2)) // "," //&
                   string(sy%propi(i)%di3_atom2(3)) // "," //&
                   string(sy%propi(i)%di3_atom2(4)) // ") "
             end if
          end if
          itaux = string(fid)
       else
          if (res(i)%done) saux = ""
          itaux = string(fid)
       end if

       write (uout,'(A2,99(A," "))') cini, string(i,4,ioj_left), string(label,12,ioj_center), &
          string(itaux,2,ioj_right), string(itype_names(sy%propi(i)%itype),12,ioj_left),&
          string(saux)
    end do
    write (uout,*)

    ! key
    write (uout,'("* Key for the interpretation of table headings")')
    if (bas%imtype == imtype_isosurface) then
       write (uout,'("# Id = isosurface identifier")')
       if (bas%higher) then
          write (uout,'("# Position = region maximum position")')
          write (uout,'("# Field = value of the reference field at the maximum")')
       else
          write (uout,'("# Position = region minimum position")')
          write (uout,'("# Field = value of the reference field at the minimum")')
       end if
       write (uout,'("# Volume = volume of the isosurface region.")')
       write (uout,'("# Pop (population) = the reference field integrated in the isosurface region.")')
       write (uout,'("# Lap = Laplacian of the reference field integrated in the isosurface region.")')
       write (uout,'("# $xxx = field xxx integrated in the isosurface region.")')
    else
       write (uout,'("# Id = attractor identifier")')
       write (uout,'("# cp/at = critical point/atom from the complete CP/atom list")')
       write (uout,'("# (lvec) = lattice translation from CP/atom to the main cell''s representative")')
       write (uout,'("# ncp/nat = critical point/atom from the non-equivalent CP/atom list")')
       write (uout,'("# Name = atomic name")')
       write (uout,'("# Z = atomic number")')
       write (uout,'("# Position = attractor position")')
       write (uout,'("# Mol = molecule identifier (see corresponding table above)")')
       write (uout,'("# Volume = volume of the reference field basins.")')
       write (uout,'("# Pop (population) = the reference field integrated in its own basins.")')
       write (uout,'("# Lap = Laplacian of the reference field integrated in the reference field basins.")')
       write (uout,'("# $xxx = field xxx integrated in the basins of the reference field.")')
    end if
    write (uout,*)

    if (bas%imtype == imtype_isosurface) then
       ! isosurface output
       if (bas%higher) then
          write (uout,'("* List of isosurface regions and associated maxima")')
       else
          write (uout,'("* List of isosurface regions and associated minima")')
       end if
       if (.not.sy%c%ismolecule) then
          write (uout,'("# Id           Position (cryst.)                 Field")')
       else
          write (uout,'("# Id           Position (",A,")                  Field")') iunitname0(iunit)
       endif
       do i = 1, bas%nattr
          if (.not.sy%c%ismolecule) then
             x = bas%xattr(:,i)
          else
             x = (sy%c%x2c(bas%xattr(:,i)) + sy%c%molx0) * dunit0(iunit)
          endif
          id = nint(bas%xattr(:,i) * sy%f(sy%iref)%grid%n)
          id = modulo(id,sy%f(sy%iref)%grid%n) + 1
          write (uout,'("  ",99(A," "))') string(i,4,ioj_left), (string(x(j),'f',12,7,4),j=1,3), &
             string(sy%f(sy%iref)%grid%f(id(1),id(2),id(3)),'e',14,7)
       end do
       write (uout,*)
    else
       ! YT/BADER/HIRSHFELD/VORONOI output

       ! List of attractors and positions
       write (uout,'("* List of attractors integrated")')
       if (.not.sy%c%ismolecule) then
          write (uout,'("# Id   cp   ncp   Name  Z   mult           Position (cryst.) ")')
       else
          write (uout,'("# Id   cp   ncp   Name  Z   mult           Position (",A,") ")') iunitname0(iunit)
       endif
       do i = 1, bas%nattr
          call assign_strings(i,bas%icp(i),usesym,scp,sncp,sname,smult,sz)
          if (.not.sy%c%ismolecule) then
             x = bas%xattr(:,i)
          else
             x = (sy%c%x2c(bas%xattr(:,i)) + sy%c%molx0) * dunit0(iunit)
          endif
          write (uout,'("  ",99(A," "))') string(i,4,ioj_left), scp, sncp, sname, sz, &
             smult, (string(x(j),'f',12,7,4),j=1,3)
       end do
       write (uout,*)

       ! List of molecules and positions
       if (.not.nomol .and. .not.sy%c%ismolecule .and. all(sy%c%mol(1:sy%c%nmol)%discrete)) then
          allocate(idxmol(2,bas%nattr))
          ! Assign attractors to molecules
          idxmol = 0
          do i = 1, bas%nattr
             jlo: do j = 1, sy%c%nmol
                do k = 1, sy%c%mol(j)%nat
                   if (bas%icp(i) == sy%c%mol(j)%at(k)%cidx) then
                      idxmol(1,i) = j
                      idxmol(2,i) = k
                      exit jlo
                   end if
                end do
             end do jlo
          end do

          ! List of molecules and associated attractors
          write (uout,'("* List of molecules integrated")')
          write (uout,'("+ Id   at (lvec)   nat    Name  Z              Position (cryst.) ")')
          do i = 1, sy%c%nmol
             xcm = sy%c%mol(i)%cmass()
             xcm = sy%c%c2x(xcm)
             ! name the molecule
             write (uout,'("# Molecule ",A," with ",A," atoms at ",3(A," "))') string(i), &
                string(sy%c%mol(i)%nat), (string(xcm(j),'f',10,6,3),j=1,3)

             ! Atomic composition
             do j = 1, bas%nattr
                if (idxmol(1,j) /= i) cycle
                call assign_strings(j,bas%icp(j),usesym,scp,sncp,sname,smult,sz)
                x = sy%c%mol(idxmol(1,j))%at(idxmol(2,j))%x
                write (uout,'(A," ",A,"(",2(A," "),A,")"," ",99(A," "))') &
                   string(j,4,ioj_left), scp, (string(sy%c%mol(idxmol(1,j))%at(idxmol(2,j))%lvec(k),2,ioj_right),k=1,3),&
                   sncp, sname, sz, (string(x(k),'f',12,7,4),k=1,3)
             end do
          end do
          write (uout,*)
          deallocate(idxmol)
       end if
    end if

  end subroutine int_output_header

  !> Write to output the result of integrating scalar fields in the
  !> atomic basins. bas = integration driver data, res(1:npropi) =
  !> results. If nomol0, prevent the output of molecular
  !> properties. If usesym0, use multiplicities.
  module subroutine int_output_fields(bas,res,nomol0,usesym0)
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_left, ioj_center
    use types, only: basindat, int_result, out_field
    type(basindat), intent(in) :: bas
    type(int_result), intent(in) :: res(:)
    logical, intent(in), optional :: nomol0
    logical, intent(in), optional :: usesym0

    integer, parameter :: ncols = 5

    integer :: i, j, k
    integer :: ip, iplast, ipmax
    integer :: nacprop(ncols)
    real*8 :: sump(sy%npropi), xmult
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    integer, allocatable :: idxmol(:,:)
    logical :: nomol, usesym

    nomol = .false.
    if (present(nomol0)) nomol = nomol0
    usesym = .false.
    if (present(usesym0)) usesym = usesym0

    ! Integrated scalar fields, atomic properties
    if (bas%imtype == imtype_isosurface) then
       write (uout,'("* Integrated isosurface regions")')
    else
       write (uout,'("* Integrated atomic properties")')
    end if
    write (uout,'("# (See key above for interpretation of column headings.)")')
    iplast = 0
    do ip = 0, (count(res(1:sy%npropi)%outmode == out_field)-1)/ncols
       ! show only the properties that have been done
       nacprop = 0
       ipmax = 0
       do i = iplast+1, sy%npropi
          if (res(i)%done .and. res(i)%outmode == out_field) then
             ipmax = ipmax + 1
             nacprop(ipmax) = i
          end if
          if (ipmax == ncols) exit
       end do
       if (ipmax == 0) exit
       iplast = nacprop(ipmax)

       ! Table header for this set of properties
       write (uout,'("# Integrable properties ",A," to ",A)') string(nacprop(1)), string(nacprop(ipmax))
       if (bas%imtype == imtype_isosurface) then
          write (uout,'("# Id   ",5(A," "))') &
             (string(sy%propi(nacprop(j))%prop_name,15,ioj_center),j=1,ipmax)
       else
          write (uout,'("# Id   cp   ncp   Name  Z   mult ",5(A," "))') &
             (string(sy%propi(nacprop(j))%prop_name,15,ioj_center),j=1,ipmax)
       end if

       ! Table rows
       sump = 0d0
       do i = 1, bas%nattr
          if (.not.bas%docelatom(bas%icp(i))) cycle
          call assign_strings(i,bas%icp(i),usesym,scp,sncp,sname,smult,sz)
          if (bas%icp(i) > 0 .and. usesym) then
             xmult = sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(bas%icp(i))%idx)%mult
          else
             xmult = 1
          endif
          ! add to the sum
          do j = 1, ipmax
             sump(nacprop(j)) = sump(nacprop(j)) + res(nacprop(j))%psum(i) * xmult
          end do
          ! table entry
          if (bas%imtype == imtype_isosurface) then
             write (uout,'("  ",99(A," "))') string(i,4,ioj_left),&
                (string(res(nacprop(j))%psum(i),'e',15,8,4),j=1,ipmax)
          else
             write (uout,'("  ",99(A," "))') string(i,4,ioj_left), scp, sncp, sname, sz, smult, &
                (string(res(nacprop(j))%psum(i),'e',15,8,4),j=1,ipmax)
          end if
       end do
       if (bas%imtype == imtype_isosurface) then
          write (uout,'(7("-"),99(A))') ("----------------",j=1,ipmax)
          write (uout,'("  Sum  ",99(A," "))') &
             (string(sump(nacprop(j)),'e',15,8,4),j=1,ipmax)
       else
          write (uout,'(32("-"),99(A))') ("----------------",j=1,ipmax)
          write (uout,'("  Sum                           ",99(A," "))') &
             (string(sump(nacprop(j)),'e',15,8,4),j=1,ipmax)
       end if
       write (uout,*)
    end do

    if (bas%imtype /= imtype_isosurface .and. .not.nomol .and. .not.sy%c%ismolecule .and.&
       all(sy%c%mol(1:sy%c%nmol)%discrete)) then
       allocate(idxmol(2,bas%nattr))
       ! Assign attractors to molecules
       idxmol = 0
       do i = 1, bas%nattr
          jlo: do j = 1, sy%c%nmol
             do k = 1, sy%c%mol(j)%nat
                if (bas%icp(i) == sy%c%mol(j)%at(k)%cidx) then
                   idxmol(1,i) = j
                   idxmol(2,i) = k
                   exit jlo
                end if
             end do
          end do jlo
       end do

       ! List of integrated properties in the molecules
       write (uout,'("* Integrated molecular properties")')
       write (uout,'("# (See key above for interpretation of column headings.)")')
       iplast = 0
       do ip = 0, (count(res(1:sy%npropi)%outmode == out_field)-1)/ncols
          ! show only the properties that have been done
          nacprop = 0
          ipmax = 0
          do i = iplast+1, sy%npropi
             if (res(i)%done .and. res(i)%outmode == out_field) then
                ipmax = ipmax + 1
                nacprop(ipmax) = i
             end if
             if (ipmax == ncols) exit
          end do
          if (ipmax == 0) exit
          iplast = nacprop(ipmax)

          ! Table header for this set of properties
          write (uout,'("# Integrable properties ",A," to ",A)') string(nacprop(1)), string(nacprop(ipmax))
          write (uout,'("# Mol ",5(A," "))') (string(sy%propi(nacprop(j))%prop_name,15,ioj_center),j=1,ipmax)

          ! Table rows
          do k = 1, sy%c%nmol+1
             if (k == sy%c%nmol+1 .and. all(idxmol > 0)) cycle
             sump = 0d0
             do i = 1, bas%nattr
                if (.not.bas%docelatom(bas%icp(i))) cycle
                if (idxmol(1,i) /= mod(k,sy%c%nmol+1)) cycle
                ! add to the sum
                if (bas%icp(i) > 0 .and. usesym) then
                   xmult = sy%f(sy%iref)%cp(sy%f(sy%iref)%cpcel(bas%icp(i))%idx)%mult
                else
                   xmult = 1
                endif
                do j = 1, ipmax
                   sump(nacprop(j)) = sump(nacprop(j)) + res(nacprop(j))%psum(i) * xmult
                end do
             end do
             ! table entry
             if (k < sy%c%nmol+1) then
                write (uout,'("  ",99(A," "))') string(k,4,ioj_left), (string(sump(nacprop(j)),'e',15,8,4),j=1,ipmax)
             else
                write (uout,'("  ",99(A," "))') "????", (string(sump(nacprop(j)),'e',15,8,4),j=1,ipmax)
             end if
          end do
          write (uout,*)
       end do
       deallocate(idxmol)
    end if

  end subroutine int_output_fields

  !> The attractors coming out of YT and BADER are not in order
  !> compatible with the crystal structure. Reorder them, including
  !> the weights of the YT. On output, gives the identity of the
  !> attractors (icp) in the complete CP list. bas contains the
  !> integration information. Others: just set up the icp array
  !> and exit.
  subroutine int_reorder_gridout(ff,bas)
    use fieldmod, only: field, type_grid
    use tools_io, only: ferror, faterr, fopen_scratch, fclose
    use types, only: realloc, basindat
    use param, only: icrd_crys
    type(field), intent(inout) :: ff
    type(basindat), intent(inout) :: bas

    integer :: i, j
    integer, allocatable :: idgaux(:,:,:), assigned(:)
    integer :: nn, nid, nattr0, luw2, n(3), nbasin, nvec
    real*8 :: dist
    integer, allocatable :: nlo(:), ibasin(:), ibasin2(:), iio(:), inear(:,:)
    real*8, allocatable :: fnear(:,:), xattr(:,:)

    if (ff%type /= type_grid) &
       call ferror("int_reorder_gridout","BADER/YT/ISOSURFACE/... can only be used with grids",faterr)

    if (bas%imtype == imtype_isosurface) then
       ! isosurfaces not associated to atoms
       if (allocated(bas%icp)) deallocate(bas%icp)
       allocate(bas%icp(bas%nattr))
       bas%icp = 0
       return
    elseif (bas%imtype == imtype_hirshfeld .or. bas%imtype == imtype_voronoi) then
       ! hirshfeld and voronoi always associated to atoms
       if (allocated(bas%icp)) deallocate(bas%icp)
       allocate(bas%icp(bas%nattr))
       do i = 1, bas%nattr
          bas%icp(i) = i
       end do
       return
    end if
    n = ff%grid%n

    ! reorder the maxima and assign maxima to atoms according to ratom
    if (allocated(bas%icp)) deallocate(bas%icp)
    allocate(bas%icp(bas%nattr),xattr(3,bas%nattr),assigned(bas%nattr))
    bas%icp = 0
    assigned = 0
    nattr0 = bas%nattr
    ! assign attractors to atoms
    if (bas%atexist) then
       do i = 1, nattr0
          nid = ff%c%identify_atom(bas%xattr(:,i),icrd_crys,distmax=bas%ratom)
          if (nid > 0) then
             assigned(i) = nid
          else
             ! maybe the closest point is a known nnm
             call ff%nearest_cp(bas%xattr(:,i),nid,dist,type=ff%typnuc)
             if (dist < bas%ratom) then
                assigned(i) = nid
             end if
          end if
       end do
    end if
    ! create the new known attractors in the correct order
    bas%nattr = 0
    do i = 1, ff%ncpcel
       if (any(assigned == i)) then
          bas%nattr = bas%nattr + 1
          bas%icp(bas%nattr) = i
          xattr(:,bas%nattr) = ff%cpcel(i)%x
       endif
    end do
    ! the rest are their own nnm, add them to the CP list, accumulate
    ! using a radius equal to ratom.
    do i = 1, nattr0
       if (assigned(i) > 0) cycle
       bas%nattr = bas%nattr + 1
       assigned(i) = bas%nattr
       bas%icp(bas%nattr) = 0

       nn = 1
       xattr(:,bas%nattr) = bas%xattr(:,i)
       do j = i+1, nattr0
          if (assigned(j) > 0) cycle
          if (ff%c%are_lclose(bas%xattr(:,i),bas%xattr(:,j),bas%ratom)) then
             nn = nn + 1
             assigned(j) = bas%nattr
          end if
       end do
       call ff%addcp(ff%c%x2c(xattr(:,bas%nattr)),1d-2,1d-1,2d-1,huge(1d0),ff%typnuc)
    end do

    ! update the idg
    allocate(idgaux(size(bas%idg,1),size(bas%idg,2),size(bas%idg,3)))
    idgaux = bas%idg
    do i = 1, nattr0
       if (assigned(i) /= i) then
          where (bas%idg == i)
             idgaux = assigned(i)
          end where
       end if
    end do
    call move_alloc(idgaux,bas%idg)

    ! update the weights the YT file
    if (bas%luw /= 0) then
       ! read all the info from the scratch file
       rewind(bas%luw)
       read (bas%luw) nbasin, nn, nvec
       if (nattr0 /= nbasin) &
          call ferror('int_reorder_gridout','inconsistent number of attractors in yt checkpoint',faterr)
       allocate(ibasin(nn),ibasin2(nn),nlo(nn),inear(nvec,nn),fnear(nvec,nn),iio(nn))
       read (bas%luw) nlo
       read (bas%luw) ibasin
       read (bas%luw) iio
       read (bas%luw) inear
       read (bas%luw) fnear

       ! build a new ibasin array
       ibasin2 = ibasin
       do i = 1, nattr0
          where (ibasin == i)
             ibasin2 = assigned(i)
          end where
       end do

       ! write the new data to a new scratch file
       luw2 = fopen_scratch()
       write (luw2) nattr0, nn, nvec
       write (luw2) nlo
       write (luw2) ibasin2
       write (luw2) iio
       write (luw2) inear
       write (luw2) fnear
       call flush(luw2)
       rewind(luw2)
       call fclose(bas%luw)
       bas%luw = luw2
       deallocate(nlo,ibasin,ibasin2,iio,inear,fnear)
    end if
    deallocate(assigned)

    if (allocated(bas%xattr)) deallocate(bas%xattr)
    call realloc(xattr,3,bas%nattr)
    call move_alloc(xattr,bas%xattr)

  end subroutine int_reorder_gridout

  !> Integrate scalar fields in atomic basins. bas = integration driver
  !> data, res(1:npropi) = results.
  subroutine intgrid_fields(bas,res)
    use yt, only: ytdata, ytdata_clean, yt_weights
    use systemmod, only: sy, itype_v, itype_f, itype_fval, itype_gmod, &
       itype_lap, itype_lapval, itype_mpoles, itype_expr
    use grid3mod, only: grid3
    use fieldmod, only: type_grid
    use tools_math, only: tosphere, genrlm_real
    use tools_io, only: uout, string
    use types, only: basindat, int_result, out_mpoles, out_field
    use param, only: ifformat_as_ft_lap, ifformat_as_ft_grad
    type(basindat), intent(in) :: bas
    type(int_result), intent(inout) :: res(:)

    integer :: i1, i2, i3
    integer :: i, k, m, fid, ntot, lmax, ix
    real*8, allocatable :: fint(:,:,:), w(:,:,:), rrlm(:)
    type(grid3) :: faux
    logical :: ok, fillgrd
    logical :: plmask(sy%npropi), first
    real*8 :: lprop(sy%npropi), x(3), x2(3), padd, p(3), dv(3), r
    real*8 :: tp(2)
    type(ytdata) :: dat

    first = .true.
    ntot = bas%n(1)*bas%n(2)*bas%n(3)
    do k = 1, sy%npropi
       if (res(k)%done) cycle
       if (.not.sy%propi(k)%used) cycle
       if (sy%propi(k)%itype == itype_v) then
          ! integrate the basin volume
          if (allocated(res(k)%psum)) deallocate(res(k)%psum)
          allocate(res(k)%psum(bas%nattr))
          res(k)%psum = 0d0
          write (uout,'("+ Integrated property (number ",A,"): ",A)') string(k), string(sy%propi(k)%prop_name)

          if (first) call run_first_pass()

          ! compute weights and integrate the scalar field properties
          !$omp parallel do private(padd) firstprivate(w)
          do i = 1, bas%nattr
             if (.not.bas%docelatom(bas%icp(i))) cycle
             if (bas%imtype == imtype_yt) then
                call yt_weights(din=dat,idb=i,w=w)
                padd = sum(w) * sy%c%omega / real(ntot,8)
             else
                padd = count(bas%idg == i) * sy%c%omega / real(ntot,8)
             endif
             !$omp critical (accum)
             res(k)%psum(i) = res(k)%psum(i) + padd
             !$omp end critical (accum)
          end do
          !$omp end parallel do
          res(k)%outmode = out_field
       elseif (sy%propi(k)%itype == itype_f.or.sy%propi(k)%itype == itype_fval.or.&
          sy%propi(k)%itype == itype_gmod.or.sy%propi(k)%itype == itype_lap .or.&
          sy%propi(k)%itype == itype_lapval.or.sy%propi(k)%itype == itype_expr.or.&
          sy%propi(k)%itype == itype_mpoles) then
          ! integrate scalar fields other than volume
          fid = sy%propi(k)%fid
          if (.not.sy%goodfield(fid)) then
             res(k)%reason = "unknown or invalid field"
             cycle
          end if
          write (uout,'("+ Integrated property (number ",A,"): ",A)') string(k), string(sy%propi(k)%prop_name)

          if (first) call run_first_pass()

          if (.not.allocated(fint)) &
             allocate(fint(bas%n(1),bas%n(2),bas%n(3)))

          ! copy the grid if it is available
          fillgrd = .false.
          ok = (sy%f(fid)%type == type_grid)
          if (ok) ok = all(sy%f(fid)%grid%n == bas%n)
          if (ok) then
             if (sy%propi(k)%itype == itype_fval .or.&
                sy%propi(k)%itype == itype_f.and..not.sy%f(fid)%usecore.or.&
                sy%propi(k)%itype == itype_mpoles.and..not.sy%f(fid)%usecore) then
                fint = sy%f(fid)%grid%f
             elseif (sy%propi(k)%itype == itype_lapval .or.&
                sy%propi(k)%itype == itype_lap.and..not.sy%f(fid)%usecore) then
                call faux%fft(sy%f(fid)%grid,ifformat_as_ft_lap)
                fint = faux%f
             elseif (sy%propi(k)%itype == itype_gmod.and..not.sy%f(fid)%usecore) then
                call faux%fft(sy%f(fid)%grid,ifformat_as_ft_grad)
                fint = faux%f
             else
                fillgrd = .true.
             end if
             call faux%end()
          else
             fillgrd = .true.
          endif

          ! otherwise generate the field on the grid
          if (fillgrd) then
             plmask = .false.
             plmask(k) = .true.
             !$omp parallel do private(x,x2,lprop)
             do i1 = 1, bas%n(1)
                x(1) = real(i1-1,8) / bas%n(1)
                do i2 = 1, bas%n(2)
                   x(2) = real(i2-1,8) / bas%n(2)
                   do i3 = 1, bas%n(3)
                      x(3) = real(i3-1,8) / bas%n(3)
                      x2 = sy%c%x2c(x)
                      call sy%grdall(x2,lprop,plmask)
                      !$omp critical (write)
                      fint(i1,i2,i3) = lprop(k)
                      !$omp end critical (write)
                   end do
                end do
             end do
             !$omp end parallel do
          end if

          if (.not.sy%propi(k)%itype == itype_mpoles) then
             ! all types of integration except for multipoles
             if (allocated(res(k)%psum)) deallocate(res(k)%psum)
             allocate(res(k)%psum(bas%nattr))
             res(k)%psum = 0d0

             ! compute weights and integrate the scalar field properties
             !$omp parallel do private(padd) firstprivate(w)
             do i = 1, bas%nattr
                if (.not.bas%docelatom(bas%icp(i))) cycle
                if (bas%imtype == imtype_yt) then
                   call yt_weights(din=dat,idb=i,w=w)
                   padd = sum(w * fint) * sy%c%omega / real(ntot,8)
                else
                   padd = sum(fint,bas%idg==i) * sy%c%omega / real(ntot,8)
                endif
                !$omp critical (accum)
                res(k)%psum(i) = res(k)%psum(i) + padd
                !$omp end critical (accum)
             end do
             !$omp end parallel do
             res(k)%outmode = out_field
          else
             ! multipoles
             ! allocate temporary space
             lmax = sy%propi(k)%lmax
             allocate(rrlm((lmax+1)*(lmax+1)))
             rrlm = 0d0

             ! allocate result space
             if (allocated(res(k)%mpole)) deallocate(res(k)%mpole)
             allocate(res(k)%mpole((lmax+1)*(lmax+1),bas%nattr))
             res(k)%mpole = 0d0

             if (bas%imtype == imtype_yt) then
                w = 0d0
                !$omp parallel do private(p,dv,r,tp) firstprivate(w,rrlm)
                do m = 1, bas%nattr
                   if (.not.bas%docelatom(bas%icp(m))) cycle
                   call yt_weights(din=dat,idb=m,w=w)
                   do i1 = 1, bas%n(1)
                      do i2 = 1, bas%n(2)
                         do i3 = 1, bas%n(3)
                            if (abs(w(i1,i2,i3)) < 1d-15) cycle
                            p = real((/i1-1,i2-1,i3-1/),8)
                            dv = p/real(bas%n,8) - bas%xattr(:,m)
                            call sy%c%shortest(dv,r)
                            call tosphere(dv,r,tp)
                            call genrlm_real(lmax,r,tp,rrlm)

                            !$omp critical (accum)
                            res(k)%mpole(:,m) = res(k)%mpole(:,m) + rrlm * fint(i1,i2,i3) * w(i1,i2,i3)
                            !$omp end critical (accum)
                         end do
                      end do
                   end do
                end do
                !$omp end parallel do
                res(k)%outmode = out_mpoles
             else
                ! calcualate the multipoles, with bader and isosurface
                !$omp parallel do private(p,ix,dv,r,tp) firstprivate(rrlm)
                do i1 = 1, bas%n(1)
                   p(1) = real(i1-1,8)
                   do i2 = 1, bas%n(2)
                      p(2) = real(i2-1,8)
                      do i3 = 1, bas%n(3)
                         p(3) = real(i3-1,8)
                         ix = bas%idg(i1,i2,i3)
                         dv = p/real(bas%n,8) - bas%xattr(:,ix)
                         call sy%c%shortest(dv,r)
                         call tosphere(dv,r,tp)
                         call genrlm_real(lmax,r,tp,rrlm)

                         !$omp critical (accum)
                         res(k)%mpole(:,ix) = res(k)%mpole(:,ix) + rrlm * fint(i1,i2,i3)
                         !$omp end critical (accum)
                      end do
                   end do
                end do
                !$omp end parallel do
             endif
             res(k)%mpole = res(k)%mpole * sy%c%omega / real(ntot,8)
             deallocate(rrlm)
          end if
       else
          ! none of the above
          cycle
       end if
       res(k)%done = .true.
       res(k)%reason = ""
    end do

    ! clean up
    if (allocated(fint)) deallocate(fint)
    if (allocated(w)) deallocate(w)
    if (bas%imtype == imtype_yt) then
       call ytdata_clean(dat)
    end if

  contains
    subroutine run_first_pass()
      ! YT, get the data
      if (bas%imtype == imtype_yt) then
         allocate(w(bas%n(1),bas%n(2),bas%n(3)))
         w = 0d0
         call yt_weights(luw=bas%luw,dout=dat)
      end if

      ! do not run this again
      first = .false.
    end subroutine run_first_pass

  end subroutine intgrid_fields

  !> Integrate scalar fields in atomic basins, the outer loop goes
  !> over atoms instead of integrable properties, for hirshfeld
  !> integration only. bas = integration driver data, res(1:npropi) =
  !> results.
  subroutine intgrid_hirshfeld_fields(bas,res)
    use grid1mod, only: grid1, agrid
    use hirshfeld, only: hirsh_weights
    use systemmod, only: sy, itype_v, itype_f, itype_fval, itype_gmod, &
       itype_lap, itype_lapval, itype_mpoles, itype_expr
    use grid3mod, only: grid3
    use fieldmod, only: type_grid
    use global, only: cutrad
    use tools_io, only: uout, string, ferror, faterr
    use types, only: basindat, int_result, out_field, realloc
    use param, only: ifformat_as_ft_lap, ifformat_as_ft_grad, icrd_crys, maxzat, VSMALL
    type(basindat), intent(in) :: bas
    type(int_result), intent(inout) :: res(:)

    integer :: i, k, l, ntot, fid
    real*8, allocatable :: w(:,:,:)
    integer :: i1, i2, i3, iz
    type(grid3) :: faux
    logical :: ok, fillgrd
    logical :: plmask(sy%npropi)
    real*8 :: lprop(sy%npropi), x(3), x2(3), x0(3), xdelta(3,3)
    real*8 :: rhoa, raux1, raux2, fac, tosum
    integer :: nmap
    integer, allocatable :: idmap(:)
    real*8, allocatable :: fmap(:,:,:,:)
    integer :: nat
    integer, allocatable :: nid(:), lvec(:,:)
    real*8, allocatable :: dist(:), rcutmax(:,:)

    if (bas%imtype /= imtype_hirshfeld) return

    ! initialize
    allocate(w(bas%n(1),bas%n(2),bas%n(3)))
    ntot = bas%n(1)*bas%n(2)*bas%n(3)

    ! prepare results of integrable properties
    nmap = 0
    allocate(idmap(sy%npropi))
    idmap = 0
    do k = 1, sy%npropi
       if (res(k)%done.or..not.sy%propi(k)%used) cycle
       if (sy%propi(k)%itype == itype_v) then
          ! volume
          if (allocated(res(k)%psum)) deallocate(res(k)%psum)
          allocate(res(k)%psum(bas%nattr))
          res(k)%psum = 0d0
          res(k)%outmode = out_field
          write (uout,'("+ Integrated property (number ",A,"): ",A)') string(k), string(sy%propi(k)%prop_name)

       elseif (sy%propi(k)%itype == itype_mpoles) then
          ! multipoles (not implemented) yet
          call ferror('intgrid_hirshfeld_fields','multipoles not implemented with hirshfeld',faterr)

       elseif (sy%propi(k)%itype == itype_f.or.sy%propi(k)%itype == itype_fval.or.&
          sy%propi(k)%itype == itype_gmod.or.sy%propi(k)%itype == itype_lap .or.&
          sy%propi(k)%itype == itype_lapval.or.sy%propi(k)%itype == itype_expr) then
          ! other scalar fields

          ! check it is a good field
          fid = sy%propi(k)%fid
          if (.not.sy%goodfield(fid)) then
             res(k)%reason = "unknown or invalid field"
             cycle
          end if

          ! create space for the grid
          nmap = nmap + 1
          idmap(k) = nmap
          if (.not.allocated(fmap)) then
             allocate(fmap(bas%n(1),bas%n(2),bas%n(3),1))
          else
             call realloc(fmap,bas%n(1),bas%n(2),bas%n(3),nmap)
          end if

          ! copy the grid if it is available
          fillgrd = .false.
          ok = (sy%f(fid)%type == type_grid)
          if (ok) ok = all(sy%f(fid)%grid%n == bas%n)
          if (ok) then
             if (sy%propi(k)%itype == itype_fval .or.&
                sy%propi(k)%itype == itype_f.and..not.sy%f(fid)%usecore.or.&
                sy%propi(k)%itype == itype_mpoles.and..not.sy%f(fid)%usecore) then
                fmap(:,:,:,nmap) = sy%f(fid)%grid%f
             elseif (sy%propi(k)%itype == itype_lapval .or.&
                sy%propi(k)%itype == itype_lap.and..not.sy%f(fid)%usecore) then
                call faux%fft(sy%f(fid)%grid,ifformat_as_ft_lap)
                fmap(:,:,:,nmap) = faux%f
             elseif (sy%propi(k)%itype == itype_gmod.and..not.sy%f(fid)%usecore) then
                call faux%fft(sy%f(fid)%grid,ifformat_as_ft_grad)
                fmap(:,:,:,nmap) = faux%f
             else
                fillgrd = .true.
             end if
             call faux%end()
          else
             fillgrd = .true.
          endif

          ! otherwise generate the field on the grid
          if (fillgrd) then
             plmask = .false.
             plmask(k) = .true.
             !$omp parallel do private(x,x2,lprop)
             do i1 = 1, bas%n(1)
                x(1) = real(i1-1,8) / bas%n(1)
                do i2 = 1, bas%n(2)
                   x(2) = real(i2-1,8) / bas%n(2)
                   do i3 = 1, bas%n(3)
                      x(3) = real(i3-1,8) / bas%n(3)
                      x2 = sy%c%x2c(x)
                      call sy%grdall(x2,lprop,plmask)
                      !$omp critical (write)
                      fmap(i1,i2,i3,nmap) = lprop(k)
                      !$omp end critical (write)
                   end do
                end do
             end do
             !$omp end parallel do
          end if

          ! set up the result space
          if (allocated(res(k)%psum)) deallocate(res(k)%psum)
          allocate(res(k)%psum(bas%nattr))
          res(k)%psum = 0d0
          res(k)%outmode = out_field
          write (uout,'("+ Integrated property (number ",A,"): ",A)') string(k), string(sy%propi(k)%prop_name)
       else
          cycle
       end if
    end do

    ! grid dimensions
    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(bas%n(i),8)
    end do

    ! calculate the cutoffs
    if (allocated(rcutmax)) deallocate(rcutmax)
    allocate(rcutmax(sy%c%nspc,2))
    rcutmax = 0d0
    do i = 1, sy%c%nspc
       iz = sy%c%spc(i)%z
       if (iz == 0 .or. iz > maxzat) cycle
       if (agrid(iz)%isinit) then
          rcutmax(i,2) = min(cutrad(iz),agrid(iz)%rmax)
       else
          call ferror('intgrid_hirshfeld_overlap','hirshfeld requires atomic grids',faterr)
       end if
    end do

    ! run over grid points
    !$omp parallel do private(x0,nat,fac,rhoa,raux1,raux2,tosum) firstprivate(nid,dist,lvec)
    do i3 = 1, bas%n(3)
       do i2 = 1, bas%n(2)
          do i1 = 1, bas%n(1)
             x0 = (i1-1) * xdelta(:,1) + (i2-1) * xdelta(:,2) + (i3-1) * xdelta(:,3)
             call sy%c%list_near_atoms(x0,icrd_crys,.false.,nat,nid,dist,lvec,up2dsp=rcutmax)

             fac = 1d0 / max(bas%f(i1,i2,i3),VSMALL)

             ! run over atoms in the environment
             do i = 1, nat
                ! calculate density and accumulate
                call agrid(sy%c%spc(sy%c%atcel(nid(i))%is)%z)%interp(dist(i),rhoa,raux1,raux2)
                tosum = fac * rhoa

                !$omp critical (results)
                ! add result
                do l = 1, sy%npropi
                   if (res(l)%done.or..not.sy%propi(l)%used) cycle
                   if (sy%propi(l)%itype == itype_v) then
                      ! volume
                      res(l)%psum(nid(i)) = res(l)%psum(nid(i)) + tosum
                   elseif (idmap(l) > 0) then
                      res(l)%psum(nid(i)) = res(l)%psum(nid(i)) + tosum * fmap(i1,i2,i3,idmap(l))
                   end if
                end do
                !$omp end critical (results)
             end do ! i
          end do ! i1
       end do ! i2
    end do ! i3
    !$omp end parallel do

    ! wrap up
    do l = 1, sy%npropi
       if (res(l)%done.or..not.sy%propi(l)%used) cycle
       res(l)%psum = res(l)%psum * sy%c%omega / product(bas%n)

       res(l)%done = .true.
       res(l)%reason = ""
    end do

  end subroutine intgrid_hirshfeld_fields

  !> Calculate Hirshfeld overlap populations (bond orders):
  !>   B_AB = int wA * wB * rho(r) dr
  !> Only for grids. bas = integration driver data, res(1:npropi) = results.
  subroutine intgrid_hirshfeld_overlap(bas,res)
    use grid1mod, only: grid1, agrid
    use systemmod, only: sy, itype_hirshfeld_ovpop
    use fieldmod, only: type_grid
    use global, only: cutrad
    use tools_io, only: uout, string, ferror, faterr
    use types, only: basindat, int_result, realloc, out_hirsh_ovpop
    use param, only: maxzat, icrd_crys, VSMALL
    type(basindat), intent(in) :: bas
    type(int_result), intent(inout) :: res(:)

    integer :: i, j, l, i1, i2, i3, iz, fid
    logical :: ok, fillgrd
    integer, allocatable :: idmap(:)
    real*8, allocatable :: fmap(:,:,:,:), aux(:,:,:,:,:)
    real*8 :: x(3), x2(3), x0(3), xdelta(3,3), rhoa, rhob, raux1, raux2, tosum
    real*8 :: lprop(sy%npropi), fac
    integer :: lb(3), ub(3), lbn(3), ubn(3), nat, lt(3), nmap
    logical :: plmask(sy%npropi)
    integer, allocatable :: nid(:), lvec(:,:)
    real*8, allocatable :: dist(:), rcutmax(:,:)

    ! only for hirshfeld integration
    if (bas%imtype /= imtype_hirshfeld) return

    ! prepare results of integrable properties
    nmap = 0
    allocate(idmap(sy%npropi))
    idmap = 0
    do l = 1, sy%npropi
       if (res(l)%done) cycle
       if (.not.sy%propi(l)%used) cycle
       if (sy%propi(l)%itype /= itype_hirshfeld_ovpop) cycle

       ! check it is a good field
       fid = sy%propi(l)%fid
       if (.not.sy%goodfield(fid)) then
          res(l)%reason = "unknown or invalid field"
          cycle
       end if

       ! create space for the grid
       nmap = nmap + 1
       idmap(l) = nmap
       if (.not.allocated(fmap)) then
          allocate(fmap(bas%n(1),bas%n(2),bas%n(3),1))
       else
          call realloc(fmap,bas%n(1),bas%n(2),bas%n(3),nmap)
       end if

       ! copy the grid if it is available
       fillgrd = .false.
       ok = (sy%f(fid)%type == type_grid)
       if (ok) ok = all(sy%f(fid)%grid%n == bas%n)
       if (ok) then
          if (.not.sy%f(fid)%usecore) then
             fmap(:,:,:,nmap) = sy%f(fid)%grid%f
          else
             fillgrd = .true.
          end if
       else
          fillgrd = .true.
       endif

       ! otherwise generate the field on the grid
       if (fillgrd) then
          plmask = .false.
          plmask(l) = .true.
          !$omp parallel do private(x,x2,lprop)
          do i1 = 1, bas%n(1)
             x(1) = real(i1-1,8) / bas%n(1)
             do i2 = 1, bas%n(2)
                x(2) = real(i2-1,8) / bas%n(2)
                do i3 = 1, bas%n(3)
                   x(3) = real(i3-1,8) / bas%n(3)
                   x2 = sy%c%x2c(x)
                   call sy%grdall(x2,lprop,plmask)
                   !$omp critical (write)
                   fmap(i1,i2,i3,nmap) = lprop(l)
                   !$omp end critical (write)
                end do
             end do
          end do
          !$omp end parallel do
       end if

       ! prepare the output array
       if (allocated(res(l)%hirsh_op)) deallocate(res(l)%hirsh_op)
       allocate(res(l)%hirsh_op(bas%nattr,bas%nattr,0:0,0:0,0:0))
       res(l)%hirsh_op = 0d0
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(sy%propi(l)%prop_name)
    end do

    ! ! run the calculation
    ! do l = 1, sy%npropi
    !    if (res(l)%done) cycle
    !    if (.not.sy%propi(l)%used) cycle
    !    if (sy%propi(l)%itype /= itype_hirshfeld_ovpop) cycle

    ! prepare for output
    lb = 0
    ub = 0

    ! grid dimensions
    do i = 1, 3
       xdelta(:,i) = 0d0
       xdelta(i,i) = 1d0 / real(bas%n(i),8)
    end do

    ! calculate the cutoffs
    if (allocated(rcutmax)) deallocate(rcutmax)
    allocate(rcutmax(sy%c%nspc,2))
    rcutmax = 0d0
    do i = 1, sy%c%nspc
       iz = sy%c%spc(i)%z
       if (iz == 0 .or. iz > maxzat) cycle
       if (agrid(iz)%isinit) then
          rcutmax(i,2) = min(cutrad(iz),agrid(iz)%rmax)
       else
          call ferror('intgrid_hirshfeld_overlap','hirshfeld requires atomic grids',faterr)
       end if
    end do

    ! run over grid points
    !$omp parallel do private(x0,nat,fac,lt,rhoa,rhob,raux1,raux2,tosum) firstprivate(nid,dist,lvec)
    do i3 = 1, bas%n(3)
       do i2 = 1, bas%n(2)
          do i1 = 1, bas%n(1)
             x0 = (i1-1) * xdelta(:,1) + (i2-1) * xdelta(:,2) + (i3-1) * xdelta(:,3)
             call sy%c%list_near_atoms(x0,icrd_crys,.false.,nat,nid,dist,lvec,up2dsp=rcutmax)

             fac = 1d0 / max(bas%f(i1,i2,i3) * bas%f(i1,i2,i3),VSMALL)

             ! run over pairs of atoms in the environment
             do i = 1, nat
                do j = i, nat
                   lt = lvec(:,j) - lvec(:,i)

                   ! calculate densities and accumulate
                   call agrid(sy%c%spc(sy%c%atcel(nid(i))%is)%z)%interp(dist(i),rhoa,raux1,raux2)
                   call agrid(sy%c%spc(sy%c%atcel(nid(j))%is)%z)%interp(dist(j),rhob,raux1,raux2)
                   tosum = fac * rhoa * rhob

                   !$omp critical (results)
                   ! reallocate if necessary
                   if (-abs(lt(1))<lb(1).or.-abs(lt(2))<lb(2).or.-abs(lt(3))<lb(3).or.&
                      abs(lt(1))>ub(1).or.abs(lt(2))>ub(2).or.abs(lt(3))>ub(3)) then
                      lbn = min(lb,-abs(lt))
                      ubn = max(ub,abs(lt))
                      do l = 1, sy%npropi
                         if (idmap(l) == 0) cycle
                         allocate(aux(bas%nattr,bas%nattr,lbn(1):ubn(1),lbn(2):ubn(2),lbn(3):ubn(3)))
                         aux = 0d0
                         aux(:,:,lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)) = res(l)%hirsh_op
                         call move_alloc(aux,res(l)%hirsh_op)
                      end do
                      lb = lbn
                      ub = ubn
                   end if
                   ! add result
                   do l = 1, sy%npropi
                      if (idmap(l) == 0) cycle
                      res(l)%hirsh_op(nid(i),nid(j),lt(1),lt(2),lt(3)) = &
                         res(l)%hirsh_op(nid(i),nid(j),lt(1),lt(2),lt(3)) + tosum * fmap(i1,i2,i3,idmap(l))
                      if (i /= j) then
                         res(l)%hirsh_op(nid(j),nid(i),-lt(1),-lt(2),-lt(3)) = &
                            res(l)%hirsh_op(nid(j),nid(i),-lt(1),-lt(2),-lt(3)) + tosum * fmap(i1,i2,i3,idmap(l))
                      end if
                   end do
                   !$omp end critical (results)
                end do ! j
             end do ! i
          end do ! i1
       end do ! i2
    end do ! i3
    !$omp end parallel do

    ! wrap up
    do l = 1, sy%npropi
       if (idmap(l) == 0) cycle
       res(l)%hirsh_op(:,:,:,:,:) = res(l)%hirsh_op(:,:,:,:,:) * sy%c%omega / product(bas%n)

       res(l)%done = .true.
       res(l)%reason = ""
       res(l)%outmode = out_hirsh_ovpop
    end do

  end subroutine intgrid_hirshfeld_overlap

  !> Calculate localization and delocalization indices using
  !> grids. bas = integration driver data, res(1:npropi) = results.
  subroutine intgrid_deloc(bas,res)
    use bader, only: bader_remap
    use yt, only: yt_weights, ytdata, ytdata_clean, yt_remap
    use systemmod, only: sy, itype_deloc_wnr, itype_deloc_psink, itype_deloc_sijchk, itype_deloc_fachk
    use fieldmod, only: type_grid, type_wfn
    use crystalmod, only: crystal
    use global, only: fileroot
    use tools_io, only: uout, string, fclose
    use tools_math, only: matinv
    use types, only: basindat, realloc, int_result, out_deloc, sijtype_wnr, sijtype_psink,&
       sijtype_unknown
    type(basindat), intent(in) :: bas
    type(int_result), intent(inout) :: res(:)

    integer :: j, l, natt1
    logical :: calcsij, first, ok
    integer :: luevc(2), luevc_ibnd(2)
    integer :: fid
    integer :: nlat(3), nbnd, nbndw(2), nlattot, nmo, nspin, nattn
    integer, allocatable :: iatt(:), ilvec(:,:), idg1(:,:,:)
    type(ytdata) :: dat
    character(len=:), allocatable :: sijfname, fafname
    real*8, allocatable :: w(:,:,:), kpt(:,:), occ(:,:,:)
    integer :: isijtype

    ! only for bader integration
    if (bas%imtype /= imtype_bader .and. bas%imtype /= imtype_yt) return

    ! run over all integrable properties
    first = .true.
    do l = 1, sy%npropi
       if (res(l)%done) cycle
       if (.not.sy%propi(l)%used) cycle
       if (sy%propi(l)%itype /= itype_deloc_wnr.and.sy%propi(l)%itype /= itype_deloc_psink.and.&
          sy%propi(l)%itype /= itype_deloc_sijchk.and.sy%propi(l)%itype /= itype_deloc_fachk) cycle
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(sy%propi(l)%prop_name)

       ! check consistency of the field, if applicable
       ! assign checkpoints
       if (sy%propi(l)%itype == itype_deloc_wnr .or. sy%propi(l)%itype == itype_deloc_psink) then
          fid = sy%propi(l)%fid
          if (.not.sy%goodfield(fid)) then
             res(l)%reason = "unknown or invalid field"
             cycle
          end if
          sijfname = trim(sy%f(fid)%file) // "-sij"
          fafname = trim(sy%f(fid)%file) // "-fa"
       else
          sijfname = trim(fileroot) // ".chk-sij"
          fafname = trim(fileroot) // ".chk-fa"
       end if

       ! set the type of the sij matrix
       if (sy%propi(l)%itype == itype_deloc_wnr) then
          res(l)%sijtype = sijtype_wnr
      elseif (sy%propi(l)%itype == itype_deloc_psink) then
          res(l)%sijtype = sijtype_psink
       else
          res(l)%sijtype = sijtype_unknown
       end if

       ! maybe we can read the Fa information and jump to the end
       if (sy%propi(l)%itype == itype_deloc_fachk) then
          if (read_chk_header(sy%propi(l)%fachkfile,nbnd,nbndw,nlat,nmo,nlattot,nspin,natt1,res(l)%sijtype)) then
             if (natt1 == bas%nattr) then
                write (uout,'("# Reading Fa checkpoint file: ",A)') string(sy%propi(l)%fachkfile)
                if (allocated(res(l)%fa)) deallocate(res(l)%fa)
                allocate(res(l)%fa(bas%nattr,bas%nattr,nlattot,nspin))
                call read_fachk_body(sy%propi(l)%fachkfile,res(l)%fa)
                goto 999
             else
                res(l)%reason = "inconsistent number of attractors in Fa checkpoint"
                cycle
             end if
          else
             res(l)%reason = "Fa checkpoint file not found"
             cycle
          end if
       else if ((sy%propi(l)%itype == itype_deloc_wnr .or. sy%propi(l)%itype == itype_deloc_psink).and.&
          sy%propi(l)%fachk .and..not.sy%propi(l)%di3) then
          ! If DI3 calculation is requested, we need to read/calculate the Sij anyway
          if (read_chk_header(fafname,nbnd,nbndw,nlat,nmo,nlattot,nspin,natt1,isijtype)) then
             fid = sy%propi(l)%fid

             ok = (sy%propi(l)%itype == itype_deloc_psink) .or. all(nbndw == sy%f(fid)%grid%qe%nbndw)
             ok = ok .and. all(nlat == sy%f(fid)%grid%qe%nk)
             ok = ok .and. (nlattot == sy%f(fid)%grid%qe%nks)
             ok = ok .and. (nspin == sy%f(fid)%grid%qe%nspin)
             ok = ok .and. (nbnd == sy%f(fid)%grid%qe%nbnd)
             ok = ok .and. (natt1 == bas%nattr)
             ok = ok .and. (isijtype == res(l)%sijtype)
             if (ok) then
                write (uout,'("# Reading Fa checkpoint file: ",A)') string(fafname)
                if (allocated(res(l)%fa)) deallocate(res(l)%fa)
                allocate(res(l)%fa(bas%nattr,bas%nattr,nlattot,nspin))
                call read_fachk_body(fafname,res(l)%fa)
                goto 999
             end if
          end if
       end if

       ! maybe we can read the Sij information and bypass the Sij calculation
       calcsij = .true.
       if (sy%propi(l)%itype == itype_deloc_sijchk) then
          ! read the sij from a checkpoint file (without the corresponding field)
          if (read_chk_header(sy%propi(l)%sijchkfile,nbnd,nbndw,nlat,nmo,nlattot,nspin,natt1,res(l)%sijtype)) then
             if (natt1 == bas%nattr) then
                write (uout,'("# Reading Sij checkpoint file: ",A)') string(sy%propi(l)%sijchkfile)
                if (allocated(res(l)%sijc)) deallocate(res(l)%sijc)
                allocate(res(l)%sijc(nmo,nmo,bas%nattr,nspin),kpt(3,nlattot),occ(nbnd,nlattot,nspin))
                call read_sijchk_body(sy%propi(l)%sijchkfile,kpt,occ,res(l)%sijc)
                res(l)%nlat = nlat
                res(l)%nspin = nspin
                calcsij = .false.
             else
                res(l)%reason = "inconsistent number of attractors in Sij checkpoint"
                cycle
             end if
          else
             res(l)%reason = "Sij checkpoint file not found"
             cycle
          end if
       elseif ((sy%propi(l)%itype == itype_deloc_wnr .or. sy%propi(l)%itype == itype_deloc_psink).and.&
          sy%propi(l)%sijchk) then
          if (read_chk_header(sijfname,nbnd,nbndw,nlat,nmo,nlattot,nspin,natt1,isijtype)) then
             fid = sy%propi(l)%fid

             ok = (sy%propi(l)%itype == itype_deloc_psink) .or. all(nbndw == sy%f(fid)%grid%qe%nbndw)
             ok = ok .and. all(nlat == sy%f(fid)%grid%qe%nk)
             ok = ok .and. (nlattot == sy%f(fid)%grid%qe%nks)
             ok = ok .and. (nspin == sy%f(fid)%grid%qe%nspin)
             ok = ok .and. (nbnd == sy%f(fid)%grid%qe%nbnd)
             ok = ok .and. (natt1 == bas%nattr)
             ok = ok .and. (isijtype == res(l)%sijtype)
             if (ok) then
                write (uout,'("# Reading Sij checkpoint file: ",A)') string(sijfname)
                if (allocated(res(l)%sijc)) deallocate(res(l)%sijc)

                ! read the sij checkpoint body and check the kpts and occs
                allocate(res(l)%sijc(nmo,nmo,bas%nattr,nspin))
                allocate(kpt(3,nlattot),occ(nbnd,nlattot,nspin))
                call read_sijchk_body(sijfname,kpt,occ,res(l)%sijc)
                ok = ok .and. all(abs(kpt - sy%f(fid)%grid%qe%kpt) < 1d-4)
                ok = ok .and. all(abs(occ - sy%f(fid)%grid%qe%occ) < 1d-4)
                if (ok) calcsij = .false.
             end if
          end if
       end if ! sy%propi(l)%itype == itype_deloc, etc.

       ! assign values to some integers and check consistency of the input field
       if (sy%propi(l)%itype == itype_deloc_wnr .or.sy%propi(l)%itype == itype_deloc_psink) then
          ! check consistency of the input field
          if (sy%f(fid)%type /= type_grid) then
             if (sy%f(fid)%type /= type_wfn) &
                res(l)%reason = "cannot calculate delocalization indices with non-grid fields"
             cycle
          end if
          if (.not.sy%f(fid)%grid%isqe) then
             res(l)%reason = "QE data not available for this field"
             cycle
          end if
          if (sy%propi(l)%itype == itype_deloc_wnr.and..not.sy%f(fid)%grid%iswan.and.sy%propi(l)%useu) then
             res(l)%reason = "Wannier rotation data not available for this field"
             cycle
          end if
          if (.not.all(sy%f(fid)%grid%n == bas%n)) then
             res(l)%reason = "Wannier and reference grids have different number of points"
             cycle
          end if
          if (sy%propi(l)%itype == itype_deloc_psink.and.any(sy%f(fid)%grid%qe%nk == 0)) then
             res(l)%reason = "Cannot use a pwc from open_grid.x or gamma tricks with psink"
             cycle
          end if
          if (product(sy%f(fid)%grid%qe%nk) /= sy%f(fid)%grid%qe%nks) then
             res(l)%reason = "Cannot use a pwc with k-point symmetry"
             cycle
          end if
          if (allocated(sy%f(fid)%grid%qe%spread)) then
             if (sy%propi(l)%useu .and. any(sy%f(fid)%grid%qe%spread < 0d0)) then
                res(l)%reason = "Cannot use the U when Wannier functions have negative spread"
                cycle
             end if
          end if

          ! assign integers
          nbnd = sy%f(fid)%grid%qe%nbnd
          nlat = sy%f(fid)%grid%qe%nk
          nlattot = sy%f(fid)%grid%qe%nks
          nspin = sy%f(fid)%grid%qe%nspin
          nmo = nlattot * nbnd
          nbndw = sy%f(fid)%grid%qe%nbndw
       end if

       !!! calculate Sij !!!
       if (calcsij) then
          ! write header and allocate YT weights
          write (uout,'("# Calculating atomic overlap matrices")')

          ! get the data for YT
          if (first .and. bas%imtype == imtype_yt) then
             allocate(w(bas%n(1),bas%n(2),bas%n(3)))
             w = 0d0
             call yt_weights(luw=bas%luw,dout=dat)
             first = .false.
          end if

          ! Recalculate the number of attractors without cell translation symmetry.
          ! Calculate basin spreads.
          write (uout,'(99(A," "))') "  Attractors before remapping =", string(bas%nattr)
          if (bas%imtype == imtype_bader) then
             call bader_remap(sy,bas,nattn,idg1,ilvec,iatt)
          else
             call yt_remap(sy,bas,dat,nattn,ilvec,iatt)
          end if
          write (uout,'(99(A," "))') "  Attractors after remapping =", string(nattn)

          ! allocate the sij
          if (allocated(res(l)%sijc)) deallocate(res(l)%sijc)
          allocate(res(l)%sijc(nmo,nmo,bas%nattr,nspin))
          res(l)%sijc = 0d0

          ! write out some info
          write (uout,'(99(A," "))') "  Number of bands (nbnd) =", string(nbnd)
          write (uout,'(99(A," "))') "  ... lattice translations (nlat) =", (string(nlat(j)),j=1,3)
          write (uout,'(99(A," "))') "  ... unique grid functions (nbnd x nlat) =", string(nlattot)
          write (uout,'(99(A," "))') "  ... spin channels =", string(nspin)
          if (sy%propi(l)%itype == itype_deloc_wnr) then
             write (uout,'(99(A," "))') "  Calculating overlaps using Wannier functions (wan)"
             if (sy%propi(l)%wancut > 0d0 .and. sy%propi(l)%useu) then
                write (uout,'(99(A," "))') "  Discarding overlaps if (spr(w1)+spr(w2)) * cutoff > d(cen(w1),cen(w2)), cutoff = ",&
                   string(sy%propi(l)%wancut,'f',5,2)
             else
                write (uout,'(99(A," "))') "  Discarding no overlaps."
             end if
          elseif (sy%propi(l)%itype == itype_deloc_psink) then
             write (uout,'(99(A," "))') "  Calculating overlaps using Bloch states (psink)"
          end if
          if (sy%propi(l)%di3) then
             write (uout,'(99(A," "))') "  Calculating three-center delocalization indices"
          end if
          if (sy%propi(l)%sijrestart) &
             write (uout,'(99(A," "))') "  Restart file read/written: ", trim(sy%f(fid)%file) // "-sijrestart"

          if (sy%propi(l)%itype == itype_deloc_wnr) then
             !!! using wannier functions !!!
             ! prepare the transformed files
             luevc = -1
             luevc_ibnd = 0
             write (uout,'(99(A," "))') "  Writing temporary evc files..."
             call sy%f(fid)%grid%rotate_qe_evc(luevc,luevc_ibnd,sy%propi(l)%useu)

             ! calculate overlaps
             write (uout,'(99(A," "))') "# Calculating overlaps..."
             call calc_sij_wannier(fid,sy%propi(l)%wancut,sy%propi(l)%useu,sy%propi(l)%sijrestart,&
                bas%imtype,nattn,iatt,bas%docelatom,ilvec,idg1,bas%xattr,dat,luevc,luevc_ibnd,res(l)%sijc)

             ! close the rotated evc scratch files
             if (luevc(1) >= 0) call fclose(luevc(1))
             if (luevc(2) >= 0) call fclose(luevc(2))

          elseif (sy%propi(l)%itype == itype_deloc_psink) then
             !!! using bloch functions !!!
             write (uout,'(99(A," "))') "# Calculating overlaps..."
             call calc_sij_psink(fid,sy%propi(l)%sijrestart,bas%imtype,nattn,iatt,bas%docelatom,ilvec,&
                idg1,bas%xattr,dat,res(l)%sijc)

          end if
          deallocate(iatt,ilvec)
          if (allocated(idg1)) deallocate(idg1)

          ! write the checkpoint
          if (sy%propi(l)%sijchk) then
             write (uout,'("# Writing Sij checkpoint file: ",A)') trim(sijfname)
             call write_sijchk(sijfname,nbnd,nbndw,nlat,nmo,nlattot,nspin,bas%nattr,res(l)%sijtype,&
                sy%f(fid)%grid%qe%kpt,sy%f(fid)%grid%qe%occ,res(l)%sijc)
          end if
       end if ! calcsij

       ! calculate the Sij translation mappings
       if (sy%propi(l)%itype == itype_deloc_wnr .or. &
          sy%propi(l)%itype == itype_deloc_sijchk .and. res(l)%sijtype == sijtype_wnr) then
          write (uout,'(99(A," "))') "# Calculating translation mappings..."
          call find_sij_translations(res(l),nmo,nbnd,nlat,nlattot)
       end if

       ! ! check the sanity of the Sij matrix
       ! write (uout,'(99(A," "))') "# Checking the sanity of the Sij matrix..."
       ! call check_sij_sanity(res(l),sy%f(fid)%grid%qe,nspin,nmo,nbnd,nlat,nlattot)

       !!! calculate Fa !!!
       write (uout,'("# Calculating Fa")')

       ! calculate the fa
       if (allocated(res(l)%fa)) deallocate(res(l)%fa)
       if (res(l)%sijtype == sijtype_wnr) then
          call calc_fa_wannier(res(l),nmo,nbnd,nlat,bas%nattr,nspin)
       elseif (res(l)%sijtype == sijtype_psink) then
          call calc_fa_psink(res(l),nmo,nbnd,nlat,bas%nattr,nspin,sy%f(fid)%grid%qe%kpt,sy%f(fid)%grid%qe%occ)
       end if

       ! write (*,*) "Check Fa...", sum(res(l)%fa(:,:,:,:))
       ! write (*,*) "Atomic charges..."
       ! do i = 1, bas%nattr
       !    write (*,*) i, sum(res(l)%fa(i,:,:,:)), sum(res(l)%fa(:,i,:,:))
       ! end do
       ! write (*,*) "Spin charges..."
       ! do is = 1, nspin
       !    do i = 1, bas%nattr
       !       write (*,*) i, sum(res(l)%fa(i,:,:,is)), sum(res(l)%fa(:,i,:,is))
       !    end do
       ! end do

       ! write the checkpoint file
       if (sy%propi(l)%fachk) then
          write (uout,'("# Writing Fa checkpoint file: ",A)') trim(fafname)
          call write_fachk(fafname,nbnd,nbndw,nlat,nmo,nlattot,nspin,bas%nattr,res(l)%sijtype,res(l)%fa)
       end if

       ! calculate three-body DIs, if requested
       if (sy%propi(l)%di3) then
          write (uout,'("# Calculating three-atom delocalization indices")')
          if (res(l)%sijtype == sijtype_wnr) then
             call calc_di3_wannier(res(l),nmo,nbnd,nlat,bas%nattr,nspin,&
                sy%propi(l)%di3_atom1,sy%propi(l)%di3_atom2)
          else
             call calc_di3_psink(res(l),nmo,nbnd,nlat,bas%nattr,nspin,&
                sy%f(fid)%grid%qe%kpt,sy%f(fid)%grid%qe%occ,&
                sy%propi(l)%di3_atom1,sy%propi(l)%di3_atom2)
          end if
       end if

       ! finished the Sij and Fa successfully
999    continue

       ! deallocate
       if (allocated(res(l)%sijc)) deallocate(res(l)%sijc)
       if (allocated(res(l)%sij_wnr_imap)) deallocate(res(l)%sij_wnr_imap)

       ! wrap up
       res(l)%done = .true.
       res(l)%reason = ""
       res(l)%outmode = out_deloc
       res(l)%nlat = nlat
       res(l)%nspin = nspin
   end do ! l = 1, sy%npropi

    ! clean up
    if (bas%imtype == imtype_yt) then
       if (allocated(w)) deallocate(w)
       call ytdata_clean(dat)
    end if

  end subroutine intgrid_deloc

  !> Write the Sij checkpoint file (DI integration).
  subroutine write_sijchk(sijfname,nbnd,nbndw,nlat,nmo,nlattot,nspin,nattr,sijtype,&
     kpt,occ,sij)
    use tools_io, only: fopen_write, fclose, ferror, faterr
    character(len=*), intent(in) :: sijfname
    integer, intent(in) :: nbnd, nbndw(2), nlat(3), nmo, nlattot, nspin, nattr, sijtype
    real*8, intent(in) :: kpt(3,nlattot)
    real*8, intent(in) :: occ(nbnd,nlattot,nspin)
    complex*16, intent(in) :: sij(:,:,:,:)

    integer :: lu

    lu = fopen_write(sijfname,"unformatted")
    if (lu < 0) goto 999
    write (lu,err=999) nbnd, nbndw, nlat, nmo, nlattot, nspin, nattr, sijtype
    write (lu,err=999) kpt
    write (lu,err=999) occ
    write (lu,err=999) sij
    call fclose(lu)

    return
999 continue
    call ferror("write_sijchk","error writing sijchk file",faterr)

  end subroutine write_sijchk

  !> Write the Sij restart file (DI integration).
  subroutine write_sijrestart(sijfname,nbnd,nbndw,nwan,nmo,nspin,nattr,sijtype,is,ibnd1,ibnd2,sij)
    use tools_io, only: fopen_write, fclose, ferror, faterr
    character(len=*), intent(in) :: sijfname
    integer, intent(in) :: nbnd, nbndw(2), nwan(3), nmo, nspin, nattr, sijtype
    integer, intent(in) :: is, ibnd1, ibnd2
    complex*16, intent(in) :: sij(:,:,:,:)

    integer :: lu

    lu = fopen_write(sijfname,"unformatted")
    if (lu < 0) goto 999
    write (lu,err=999) nbnd, nbndw, nwan, nmo, nspin, nattr, sijtype
    write (lu,err=999) is, ibnd1, ibnd2
    write (lu,err=999) sij
    call fclose(lu)

    return
999 continue
    call ferror("write_sijrestart","error writing sij restart file",faterr)

  end subroutine write_sijrestart

  !> Read the Sij restart file (DI integration). Returns true if
  !> usable. If the file is valid, also return the spin and band
  !> indices for the last calculated band pair as well as the current
  !> sij.
  function read_sijrestart(sijfname,nbnd,nbndw,nwan,nmo,nspin,nattr,sijtype,is,ibnd1,ibnd2,sij) result(valid)
    use tools_io, only: fopen_read, fclose
    character(len=*), intent(in) :: sijfname
    integer, intent(in) :: nbnd, nbndw(2), nwan(3), nmo, nspin, nattr, sijtype
    integer, intent(out) :: is, ibnd1, ibnd2
    complex*16, intent(inout) :: sij(:,:,:,:)
    logical :: valid

    integer :: lu
    integer :: nbnd_, nbndw_(2), nwan_(3), nmo_, nspin_, nattr_, sijtype_

    is = -1
    ibnd1 = -1
    ibnd2 = -1
    valid = .true.
    lu = fopen_read(sijfname,"unformatted",errstop=.false.)
    if (lu < 0) goto 999
    read (lu,err=999,end=999) nbnd_, nbndw_, nwan_, nmo_, nspin_, nattr_, sijtype_
    valid = (nbnd_ == nbnd) .and. all(nbndw_ == nbndw) .and. all(nwan_ == nwan) .and.&
       (nmo_ == nmo) .and. (nspin_ == nspin) .and. (nattr_ == nattr) .and.&
       (sijtype_ == sijtype)
    if (valid) then
       read (lu,err=999,end=999) is, ibnd1, ibnd2
       read (lu,err=999) sij
    end if
    call fclose(lu)

    return
999 continue
    if (lu > 0) call fclose(lu)
    valid = .false.

  end function read_sijrestart

  !> Write the Fa checkpoint file (DI integration).
  subroutine write_fachk(fafname,nbnd,nbndw,nlat,nmo,nlattot,nspin,nattr,sijtype,fa)
    use tools_io, only: fopen_write, fclose, ferror, faterr
    character(len=*), intent(in) :: fafname
    integer, intent(in) :: nbnd, nbndw(2), nlat(3), nmo, nlattot, nspin, nattr, sijtype
    real*8, intent(in) :: fa(:,:,:,:)

    integer :: lu

    lu = fopen_write(fafname,"unformatted")
    if (lu < 0) goto 999
    write (lu,err=999) nbnd, nbndw, nlat, nmo, nlattot, nspin, nattr, sijtype
    write (lu,err=999) fa
    call fclose(lu)

    return
999 continue
    call ferror("write_fachk","error writing fachk file",faterr)

  end subroutine write_fachk

  !> Read the header for the Sij/Fa checkpoint file (DI
  !> integration). If found and read, return .true.
  function read_chk_header(fname,nbnd,nbndw,nlat,nmo,nlattot,nspin,nattr,sijtype) result(haschk)
    use tools_io, only: fopen_read, fclose
    character(len=*), intent(in) :: fname
    integer, intent(out) :: nbnd, nbndw(2), nlat(3), nmo, nlattot, nspin, nattr, sijtype
    logical :: haschk

    integer :: lu

    inquire(file=fname,exist=haschk)
    if (.not.haschk) return
    lu = fopen_read(fname,"unformatted",errstop=.false.)
    if (lu < 0) goto 999
    read (lu,err=999,end=999) nbnd, nbndw, nlat, nmo, nlattot, nspin, nattr, sijtype
    call fclose(lu)
    return

999 continue
    if (lu > 0) call fclose(lu)
    haschk = .false.

  end function read_chk_header

  !> Read the body of the Sij checkpoint file (DI integration).
  subroutine read_sijchk_body(sijfname,kpt,occ,sij)
    use tools_io, only: fopen_read, fclose, ferror, faterr
    character(len=*), intent(in) :: sijfname
    real*8, intent(inout) :: kpt(:,:)
    real*8, intent(inout) :: occ(:,:,:)
    complex*16, intent(inout) :: sij(:,:,:,:)

    integer :: lu

    lu = fopen_read(sijfname,"unformatted")
    if (lu < 0) goto 999
    read (lu,err=999,end=999)
    read (lu,err=999,end=999) kpt
    read (lu,err=999,end=999) occ
    read (lu,err=999) sij
    call fclose(lu)

    return
999 continue
    call ferror("read_sijchk_body","error reading sijchk file",faterr)

  end subroutine read_sijchk_body

  !> Read the body of the Fa checkpoint file (DI integration).
  subroutine read_fachk_body(fafname,fa)
    use tools_io, only: fopen_read, fclose, ferror, faterr
    character(len=*), intent(in) :: fafname
    real*8, intent(inout) :: fa(:,:,:,:)

    integer :: lu

    lu = fopen_read(fafname,"unformatted")
    if (lu < 0) goto 999
    read (lu,err=999,end=999)
    read (lu,err=999) fa
    call fclose(lu)

    return
999 continue
    call ferror("read_fachk_body","error reading fachk file",faterr)

  end subroutine read_fachk_body

  !> Calculate the atomic overlap matrices (sij) from the complex
  !> Wannier functions in field fid. Use integration method imtype
  !> (bader/yt), with natt1 remapped attractors, iatt = attractor
  !> mapping, ilvec = attractor lattice vector, idg1 = grid assignment
  !> to attractors in Bader, xattr = attractor position, dat = YT data
  !> type. luevc are the two scratch files for the rotated evc and
  !> luevc_ibnd are the band pointers in those files.
  subroutine calc_sij_wannier(fid,wancut,useu,restart,imtype,natt1,iatt,doatcel,&
     ilvec,idg1,xattr,dat,luevc,luevc_ibnd,sij)
    use systemmod, only: sy
    use yt, only: yt_weights, ytdata, ytdata_clean
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use tools_io, only: ferror, faterr, uout, string
    use types, only: sijtype_wnr
    integer, intent(in) :: fid
    real*8, intent(in) :: wancut
    logical, intent(in) :: useu
    logical, intent(in) :: restart
    integer, intent(in) :: imtype
    integer, intent(in) :: natt1
    integer, intent(in) :: iatt(natt1)
    logical, intent(in) :: doatcel(0:)
    integer, intent(in) :: ilvec(3,natt1)
    integer, intent(in), allocatable :: idg1(:,:,:)
    real*8, intent(in) :: xattr(:,:)
    type(ytdata), intent(in) :: dat
    integer, intent(in) :: luevc(2)
    integer, intent(inout) :: luevc_ibnd(2)
    complex*16, intent(out) :: sij(:,:,:,:)

    integer :: i, is, ibnd1, ibnd2, is_restart, ibnd1_restart, ibnd2_restart
    integer :: imo, imo1, ia, ja, ka, iba, ilata, jmo, jmo1, ib, jb, kb, ibb, ilatb
    integer :: n(3), nbnd, nbndw(2), nwan(3), nlat, nmo, nspin
    integer :: ncalc, m1, m2, m3, p(3)
    real*8 :: d0, d2, x(3), xs(3)
    logical :: okrestart
    type(crystalseed) :: ncseed
    type(crystal) :: nc
    logical, allocatable :: lovrlp(:,:,:,:,:,:)
    complex*16, allocatable :: psic(:,:,:), psic2(:,:,:)
    complex*16 :: padd
    complex*16, allocatable :: f1(:,:,:,:), f2(:,:,:,:)
    real*8, allocatable :: w(:,:,:)
    logical, allocatable :: wmask(:,:,:)
    character(len=:), allocatable :: restartname

    sij = 0d0
    n = sy%f(sy%iref)%grid%n
    nwan = sy%f(fid)%grid%qe%nk
    nlat = sy%f(fid)%grid%qe%nks
    nspin = sy%f(fid)%grid%qe%nspin
    nbnd = sy%f(fid)%grid%qe%nbnd
    nbndw = sy%f(fid)%grid%qe%nbndw
    nmo = nlat * nbnd
    restartname = trim(sy%f(fid)%file) // "-sijrestart"

    ! build the supercell
    ncseed%isused = .true.
    do i = 1, 3
       ncseed%m_x2c(:,i) = sy%c%m_x2c(:,i) * nwan(i)
    end do
    ncseed%useabr = 2
    ncseed%nat = 0
    ncseed%havesym = 0
    ncseed%findsym = 0
    ncseed%ismolecule = sy%c%ismolecule
    call nc%struct_new(ncseed,.true.)

    ! restart
    if (restart) then
       okrestart = read_sijrestart(restartname,nbnd,nbndw,nwan,nmo,nspin,natt1,sijtype_wnr,&
          is_restart,ibnd1_restart,ibnd2_restart,sij)
       if (okrestart) then
          write (uout,'("    Restarting at band (",A,",",A,") of total ",A,". Spin ",A,"/",A)') &
             string(ibnd1_restart), string(ibnd2_restart), string(nbndw(is_restart)), string(is_restart),&
             string(nspin)
       else
          sij = 0d0
       end if
    else
       okrestart = .false.
       is_restart = -1
       ibnd1_restart = -1
       ibnd2_restart = -1
    end if

    if (any(n /= sy%f(fid)%grid%n)) &
       call ferror("calc_sij_wannier","inconsistent grid sizes",faterr)
    allocate(psic(n(1),n(2),n(3)))
    allocate(f1(n(1),n(2),n(3),nlat))
    allocate(f2(n(1),n(2),n(3),nlat))
    allocate(lovrlp(0:nwan(1)-1,0:nwan(2)-1,0:nwan(3)-1,0:nwan(1)-1,0:nwan(2)-1,0:nwan(3)-1))

    ! the big loop
    if (imtype == imtype_yt) &
       allocate(w(n(1),n(2),n(3)),wmask(n(1),n(2),n(3)),psic2(n(1),n(2),n(3)))
    do is = 1, nspin
       if (okrestart .and. is < is_restart) cycle
       do ibnd1 = 1, nbndw(is)
          if (okrestart .and. is == is_restart .and. ibnd1 < ibnd1_restart) cycle
          ! first wannier function
          call sy%f(fid)%grid%get_qe_wnr(sy%c%omega,ibnd1,is,luevc,luevc_ibnd,f1)

          do ibnd2 = ibnd1, nbndw(is)
             if (okrestart .and. is == is_restart .and. ibnd1 == ibnd1_restart .and. ibnd2 <= ibnd2_restart) cycle
             ! second wannier function
             if (ibnd1 == ibnd2) then
                f2 = f1
             else
                call sy%f(fid)%grid%get_qe_wnr(sy%c%omega,ibnd2,is,luevc,luevc_ibnd,f2)
             endif

             ! lovrlp
             lovrlp = .true.
             if (wancut > 0d0 .and. useu) then
                d0 = (sy%f(fid)%grid%qe%spread(ibnd1,is)+sy%f(fid)%grid%qe%spread(ibnd2,is)) * wancut
                do imo = 1, nmo
                   call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nwan)
                   if (iba /= ibnd1) cycle
                   do jmo = 1, nmo
                      call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nwan)
                      if (ibb /= ibnd2) cycle
                      x = (sy%f(fid)%grid%qe%center(:,ibnd1,is) + (/ia,ja,ka/)&
                         - (sy%f(fid)%grid%qe%center(:,ibnd2,is) + (/ib,jb,kb/))) / real(nwan,8)
                      call nc%shortest(x,d2)
                      if (d2 > d0) &
                         lovrlp(ia,ja,ka,ib,jb,kb) = .false.
                   end do
                end do
             end if

             ncalc = 0
             if (imtype == imtype_bader) then
                ! bader integration
                psic = 0d0
                !$omp parallel do private(ia,ja,ka,iba,ib,jb,kb,ibb,padd,imo1,jmo1,ilata,ilatb) firstprivate(psic)
                do imo = 1, nmo
                   call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nwan)
                   if (iba /= ibnd1) cycle
                   ilata = 1 + ka + nwan(3) * (ja + nwan(2) * ia)
                   do jmo = 1, nmo
                      call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nwan)
                      if (ibb /= ibnd2) cycle
                      if (.not.lovrlp(ia,ja,ka,ib,jb,kb)) cycle
                      ilatb = 1 + kb + nwan(3) * (jb + nwan(2) * ib)

                      psic = conjg(f1(:,:,:,ilata)) * f2(:,:,:,ilatb)
                      do i = 1, natt1
                         if (.not.doatcel(iatt(i))) cycle
                         padd = sum(psic,idg1==i)
                         call packidx(ia-ilvec(1,i),ja-ilvec(2,i),ka-ilvec(3,i),iba,imo1,nmo,nbnd,nwan)
                         call packidx(ib-ilvec(1,i),jb-ilvec(2,i),kb-ilvec(3,i),ibb,jmo1,nmo,nbnd,nwan)
                         !$omp critical (add)
                         ncalc = ncalc + 1
                         sij(imo1,jmo1,iatt(i),is) = sij(imo1,jmo1,iatt(i),is) + padd
                         if (ibnd1 /= ibnd2) then
                            sij(jmo1,imo1,iatt(i),is) = sij(jmo1,imo1,iatt(i),is) + conjg(padd)
                         end if
                         !$omp end critical (add)
                      end do ! natt1
                   end do ! jmo
                end do ! imo
                !$omp end parallel do
             else
                ! yt integration
                psic = 0d0
                psic2 = 0d0
                w = 0d0
                wmask = .false.
                !$omp parallel do private(p,x,xs,d2,ia,ja,ka,iba,ib,jb,kb,ibb,padd,imo1,jmo1,ilata,ilatb)&
                !$omp firstprivate(psic,psic2,w,wmask)
                do i = 1, natt1
                   if (.not.doatcel(iatt(i))) cycle
                   call yt_weights(din=dat,idb=iatt(i),w=w)
                   wmask = .false.
                   do m3 = 1, n(3)
                      do m2 = 1, n(2)
                         do m1 = 1, n(1)
                            if (abs(w(m1,m2,m3)) < 1d-15) cycle
                            p = (/m1,m2,m3/)
                            x = real(p-1,8) / n - xattr(:,iatt(i))
                            xs = x
                            call sy%c%shortest(xs,d2)
                            p = nint(x - sy%c%c2x(xs))
                            wmask(m1,m2,m3) = all(p == ilvec(:,i))
                         end do
                      end do
                   end do

                   psic2 = 0d0
                   do imo = 1, nmo
                      call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nwan)
                      if (iba /= ibnd1) cycle
                      ilata = 1 + ka + nwan(3) * (ja + nwan(2) * ia)
                      where (wmask)
                         psic2 = conjg(f1(:,:,:,ilata)) * w
                      end where

                      do jmo = 1, nmo
                         call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nwan)
                         if (ibb /= ibnd2) cycle
                         if (.not.lovrlp(ia,ja,ka,ib,jb,kb)) cycle
                         ilatb = 1 + kb + nwan(3) * (jb + nwan(2) * ib)
                         where (wmask)
                            psic =  psic2 * f2(:,:,:,ilatb)
                         end where

                         padd = sum(psic,wmask)
                         call packidx(ia-ilvec(1,i),ja-ilvec(2,i),ka-ilvec(3,i),iba,imo1,nmo,nbnd,nwan)
                         call packidx(ib-ilvec(1,i),jb-ilvec(2,i),kb-ilvec(3,i),ibb,jmo1,nmo,nbnd,nwan)
                         !$omp critical (add)
                         ncalc = ncalc + 1
                         sij(imo1,jmo1,iatt(i),is) = sij(imo1,jmo1,iatt(i),is) + padd
                         if (ibnd1 /= ibnd2) then
                            sij(jmo1,imo1,iatt(i),is) = sij(jmo1,imo1,iatt(i),is) + conjg(padd)
                         end if
                         !$omp end critical (add)
                      end do ! jmo
                   end do ! imo
                end do ! i = 1, natt1
                !$omp end parallel do
             end if ! imtype == bader/yt


             if (restart) &
                call write_sijrestart(restartname,nbnd,nbndw,nwan,nmo,nspin,natt1,sijtype_wnr,is,ibnd1,ibnd2,sij)
             write (uout,'("    Bands (",A,",",A,") of total ",A,". Spin ",A,"/",A,". Overlaps: ",A,"/",A)') &
                string(ibnd1), string(ibnd2), string(nbndw(is)), string(is), string(nspin),&
                string(ncalc), string(natt1*nlat*nlat)
          end do ! ibnd2
       end do ! ibnd1
    end do ! is

    ! clean up
    deallocate(f1,f2,psic,lovrlp)
    if (imtype == imtype_yt) &
       deallocate(w,wmask,psic2)

    ! scale
    sij = sij * sy%c%omega / (n(1)*n(2)*n(3))

  end subroutine calc_sij_wannier

  !> Calculate the atomic overlap matrices (sij) from the Bloch
  !> functions in field fid. Use integration method imtype (bader/yt),
  !> with natt1 remapped attractors, iatt = attractor mapping, ilvec =
  !> attractor lattice vector, idg1 = grid assignment to attractors in
  !> Bader, xattr = attractor position, dat = YT data type.
  subroutine calc_sij_psink(fid,restart,imtype,natt1,iatt,doatcel,ilvec,idg1,xattr,dat,sij)
    use systemmod, only: sy
    use yt, only: yt_weights, ytdata, ytdata_clean
    use crystalmod, only: crystal
    use tools_io, only: ferror, faterr, uout, string, fopen_read, fclose
    use types, only: sijtype_psink
    use param, only: tpi, img
    integer, intent(in) :: fid
    logical, intent(in) :: restart
    integer, intent(in) :: imtype
    integer, intent(in) :: natt1
    integer, intent(in) :: iatt(natt1)
    logical, intent(in) :: doatcel(0:)
    integer, intent(in) :: ilvec(3,natt1)
    integer, intent(in) :: idg1(:,:,:)
    real*8, intent(in) :: xattr(:,:)
    type(ytdata), intent(in) :: dat
    complex*16, intent(out) :: sij(:,:,:,:)

    integer :: n(3), is, ik, ibnd, ik1, ibnd1, ik2, ibnd2, nks, nbnd, nspin
    integer :: imo1, imo2, nmo, nbndw(2)
    integer :: i, m1, m2, m3, p(3), nlat(3)
    real*8 :: x(3), xs(3), d2, kdif(3)
    complex*16, allocatable :: psi1(:,:,:), psi2(:,:,:), psic(:,:,:), evcall(:,:,:,:), rseq(:)
    real*8, allocatable :: w(:,:,:)
    logical, allocatable :: wmask(:,:,:)
    integer :: luc, ireg
    logical :: okrestart
    integer :: is_restart, imo1_restart, idum
    character(len=:), allocatable :: restartname

    ! initialize
    sij = 0d0
    n = sy%f(sy%iref)%grid%n
    nlat = sy%f(fid)%grid%qe%nk
    nks = sy%f(fid)%grid%qe%nks
    nspin = sy%f(fid)%grid%qe%nspin
    nbnd = sy%f(fid)%grid%qe%nbnd
    nbndw = sy%f(fid)%grid%qe%nbndw
    nmo = nks * nbnd
    restartname = trim(sy%f(fid)%file) // "-sijrestart"

    ! check consistency
    if (any(n /= sy%f(fid)%grid%n)) &
       call ferror("calc_sij_psink","inconsistent grid sizes",faterr)

    ! open the pwc file
    luc = fopen_read(sy%f(fid)%grid%qe%fpwc,form="unformatted")

    ! skip the header of the pwc file
    if (sy%f(fid)%grid%qe%gamma_only) then
       ireg = 18
    else
       ireg = 17
    end if

    ! read all the evc from the pwc file and close it
    allocate(evcall(maxval(sy%f(fid)%grid%qe%ngk(1:nks)),nbnd,nks,nspin))
    evcall = 0d0
    do i = 1, ireg
       read (luc)
    end do
    do is = 1, nspin
       do ik = 1, nks
          do ibnd = 1, nbnd
             read (luc) evcall(1:sy%f(fid)%grid%qe%ngk(ik),ibnd,ik,is)
          end do
       end do
    end do
    call fclose(luc)

    ! allocate work space
    allocate(psi1(n(1),n(2),n(3)),psi2(n(1),n(2),n(3)),rseq(n(1)*n(2)*n(3)))

    ! restart
    if (restart) then
       okrestart = read_sijrestart(restartname,nbnd,nbndw,nlat,nmo,nspin,natt1,sijtype_psink,&
          is_restart,imo1_restart,idum,sij)
       if (okrestart) then
          ibnd1 = modulo(imo1_restart-1,nbnd) + 1
          ik1 = (imo1_restart-1) / nbnd + 1
          write (uout,'("    Restarting at state k=",A," n=",A," of total (k=",A,",n=",A,"). Spin ",A,"/",A)') &
             string(ik1), string(ibnd1), string(nks), string(nbndw(is_restart)), string(is_restart), string(nspin)
       else
          sij = 0d0
       end if
    else
       okrestart = .false.
       is_restart = -1
       imo1_restart = -1
    end if

    ! the big loop
    if (imtype == imtype_yt) &
       allocate(w(n(1),n(2),n(3)),wmask(n(1),n(2),n(3)),psic(n(1),n(2),n(3)))
    do is = 1, nspin
       if (okrestart .and. is < is_restart) cycle
       do imo1 = 1, nmo
          if (okrestart .and. is == is_restart .and. imo1 <= imo1_restart) cycle
          ibnd1 = modulo(imo1-1,nbnd) + 1
          ik1 = (imo1-1) / nbnd + 1
          if (ibnd1 > nbndw(is)) cycle

          rseq = 0d0
          rseq(sy%f(fid)%grid%qe%nl(sy%f(fid)%grid%qe%igk_k(1:sy%f(fid)%grid%qe%ngk(ik1),ik1))) = &
             evcall(1:sy%f(fid)%grid%qe%ngk(ik1),ibnd1,ik1,is)
          psi1 = reshape(rseq,shape(psi1))
          call cfftnd(3,n,+1,psi1)

          !$omp parallel do private(ibnd2,ik2,kdif,p,x,xs,d2) &
          !$omp             firstprivate(rseq,psi2,psic,w,wmask)
          do imo2 = imo1, nmo
             ibnd2 = modulo(imo2-1,nbnd) + 1
             ik2 = (imo2-1) / nbnd + 1
             if (ibnd2 > nbndw(is)) cycle

             rseq = 0d0
             rseq(sy%f(fid)%grid%qe%nl(sy%f(fid)%grid%qe%igk_k(1:sy%f(fid)%grid%qe%ngk(ik2),ik2))) = &
                evcall(1:sy%f(fid)%grid%qe%ngk(ik2),ibnd2,ik2,is)
             psi2 = reshape(rseq,shape(psi2))
             call cfftnd(3,n,+1,psi2)

             psic = conjg(psi1) * psi2
             kdif = sy%f(fid)%grid%qe%kpt(:,ik1) - sy%f(fid)%grid%qe%kpt(:,ik2)
             do m3 = 1, n(3)
                do m2 = 1, n(2)
                   do m1 = 1, n(1)
                      psic(m1,m2,m3) = psic(m1,m2,m3) * exp(-tpi*img*(&
                         kdif(1)*(real(m1-1,8)/real(n(1),8))+&
                         kdif(2)*(real(m2-1,8)/real(n(2),8))+&
                         kdif(3)*(real(m3-1,8)/real(n(3),8))))
                   end do
                end do
             end do

             if (imtype == imtype_bader) then
                do i = 1, natt1
                   if (doatcel(iatt(i))) then
                      sij(imo1,imo2,iatt(i),is) = sij(imo1,imo2,iatt(i),is) + &
                         sum(psic,idg1==i) * exp(tpi*img*(kdif(1)*ilvec(1,i)+kdif(2)*ilvec(2,i)+kdif(3)*ilvec(3,i)))
                   end if
                end do ! natt1
             else
                ! yt integration
                w = 0d0
                wmask = .false.
                do i = 1, natt1
                   if (.not.doatcel(iatt(i))) cycle
                   call yt_weights(din=dat,idb=iatt(i),w=w)
                   wmask = .false.
                   do m3 = 1, n(3)
                      do m2 = 1, n(2)
                         do m1 = 1, n(1)
                            if (abs(w(m1,m2,m3)) < 1d-15) cycle
                            p = (/m1,m2,m3/)
                            x = real(p-1,8) / n - xattr(:,iatt(i))
                            xs = x
                            call sy%c%shortest(xs,d2)
                            p = nint(x - sy%c%c2x(xs))
                            wmask(m1,m2,m3) = all(p == ilvec(:,i))
                         end do
                      end do
                   end do

                   where (wmask)
                      psi2 = w * psic
                   end where

                   sij(imo1,imo2,iatt(i),is) = sij(imo1,imo2,iatt(i),is) + &
                      sum(psi2,wmask) * exp(tpi*img*(kdif(1)*ilvec(1,i)+kdif(2)*ilvec(2,i)+kdif(3)*ilvec(3,i)))
                end do
             end if ! imtype == bader/yt
          end do ! imo2
          !$omp end parallel do
          write (uout,'("    State k=",A," n=",A," of total (k=",A,",n=",A,"). Spin ",A,"/",A)') &
             string(ik1), string(ibnd1), string(nks), string(nbndw(is)), string(is), string(nspin)
          if (restart) &
             call write_sijrestart(restartname,nbnd,nbndw,nlat,nmo,nspin,natt1,sijtype_psink,is,imo1,-1,sij)
       end do ! imo1
    end do ! is

    ! clean up
    deallocate(evcall,psi1,psi2)
    if (allocated(w)) deallocate(w)
    if (allocated(wmask)) deallocate(wmask)
    if (allocated(psic)) deallocate(psic)

    ! other half of the sij matrix
    do is = 1, nspin
       do imo1 = 1, nmo
          do imo2 = imo1+1, nmo
             sij(imo2,imo1,:,is) = conjg(sij(imo1,imo2,:,is))
          end do
       end do
    end do

    ! scale
    sij = sij / (n(1)*n(2)*n(3))

  end subroutine calc_sij_psink

  !> Calculate the Fa matrix from the Sij matrix, wannier version
  subroutine calc_fa_wannier(res,nmo,nbnd,nlat,nattr,nspin)
    use types, only: int_result
    type(int_result), intent(inout) :: res
    integer, intent(in) :: nmo, nbnd, nlat(3), nattr, nspin

    integer :: iat_a, iat_b, rat_b, nlattot
    integer :: is, imo, jmo, ia, ja, ka, iba, ib, jb, kb, ibb
    real*8 :: fac

    real*8, allocatable :: ftemp(:,:,:)
    logical, allocatable :: useovrlp(:,:,:)
    real*8, parameter :: ovrlp_thr = 1d-50

    ! number of lattice vectors
    nlattot = nlat(1)*nlat(2)*nlat(3)

    ! flag which overlaps will be used
    allocate(useovrlp(nmo,nmo,nspin))
    useovrlp = .true.
    do is = 1, nspin
       do imo = 1, nmo
          do jmo = imo, nmo
             useovrlp(imo,jmo,is) = (maxval(abs(res%sijc(imo,jmo,:,is))) >= ovrlp_thr)
          end do
       end do
    end do

    ! calculate Fa matrix
    allocate(res%fa(nattr,nattr,nlat(1)*nlat(2)*nlat(3),nspin),ftemp(nattr,nattr,nlattot))
    ftemp = 0d0
    res%fa = 0d0
    do is = 1, nspin
       !$omp parallel do private(ia,ja,ka,iba,ib,jb,kb,ibb,fac) firstprivate(ftemp)
       do imo = 1, nmo
          call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nlat)
          do jmo = imo, nmo
             call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nlat)
             if (.not.useovrlp(imo,jmo,is)) cycle

             ftemp = 0d0
             do iat_a = 1, nattr
                do iat_b = 1, nattr
                   do rat_b = 1, nlattot
                      ftemp(iat_a,iat_b,rat_b) = &
                         real(res%sijc(jmo,imo,iat_a,is) * &
                         res%sijc(res%sij_wnr_imap(imo,rat_b),res%sij_wnr_imap(jmo,rat_b),iat_b,is),8)
                   end do
                end do
             end do

             if (imo == jmo) then
                fac = 1d0
             else
                fac = 2d0
             end if
             !$omp critical (addfa)
             res%fa(:,:,:,is) = res%fa(:,:,:,is) + fac * ftemp
             !$omp end critical (addfa)
          end do
       end do
       !$omp end parallel do
    end do

  end subroutine calc_fa_wannier

  !> Calculate the three-body DIs, wannier version
  subroutine calc_di3_wannier(res,nmo,nbnd,nlat,nattr,nspin,atom1,atom2)
    use types, only: int_result
    type(int_result), intent(inout) :: res
    integer, intent(in) :: nmo, nbnd, nlat(3), nattr, nspin
    integer, intent(in) :: atom1, atom2(4)

    integer :: i, j, k, n
    integer :: iat_a, iat_b, iat_c, rat_b, rat_c
    integer :: iat_a_ini, iat_a_end, iat_b_ini, iat_b_end, rat_b_ini, rat_b_end
    integer :: is, imo, jmo, kmo, ia, ja, ka, iba, ib, jb, kb, ibb
    integer :: ic, jc, kc, icc, nlattot
    real*8 :: fac

    real*8, allocatable :: f3temp(:,:,:,:,:)
    logical, allocatable :: useovrlp(:,:,:)
    real*8, parameter :: ovrlp_thr = 1d-50

    ! number of lattice vectors
    nlattot = nlat(1)*nlat(2)*nlat(3)

    ! set limits for the sums
    iat_a_ini = 1
    iat_a_end = nattr
    if (atom1 > 0 .and. atom1 <= nattr) then
       iat_a_ini = atom1
       iat_a_end = atom1
    end if
    iat_b_ini = 1
    iat_b_end = nattr
    rat_b_ini = 1
    rat_b_end = nlattot
    if (atom2(1) > 0 .and. atom2(1) <= nattr) then
       iat_b_ini = atom2(1)
       iat_b_end = atom2(1)
       n = 0
       loopi: do i = 0, nlat(1)-1
          do j = 0, nlat(2)-1
             do k = 0, nlat(3)-1
                n = n + 1
                if (i == modulo(atom2(2),nlat(1)) .and. j == modulo(atom2(3),nlat(2)) .and.&
                   k == modulo(atom2(4),nlat(3))) then
                   rat_b_ini = n
                   rat_b_end = n
                   exit loopi
                end if
             end do
          end do
       end do loopi
    end if

    ! flag which overlaps will be used
    allocate(useovrlp(nmo,nmo,nspin))
    useovrlp = .true.
    do is = 1, nspin
       do imo = 1, nmo
          do jmo = imo, nmo
             useovrlp(imo,jmo,is) = (maxval(abs(res%sijc(imo,jmo,:,is))) >= ovrlp_thr)
          end do
       end do
    end do

    ! three-atom indices
    allocate(res%fa3(nattr,nattr,nlattot,nattr,nlattot,nspin))
    allocate(f3temp(nattr,nattr,nlattot,nattr,nlattot))
    res%fa3 = 0d0
    f3temp = 0d0
    do is = 1, nspin
       !$omp parallel do private(ia,ja,ka,iba,ib,jb,kb,ibb,ic,jc,kc,icc,fac) firstprivate(f3temp)
       do imo = 1, nmo
          call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nlat)
          do jmo = 1, nmo
             call unpackidx(jmo,ib,jb,kb,ibb,nmo,nbnd,nlat)
             if (.not.useovrlp(imo,jmo,is)) cycle
             do kmo = 1, nmo
                call unpackidx(jmo,ic,jc,kc,icc,nmo,nbnd,nlat)
                if (.not.useovrlp(imo,kmo,is)) cycle
                if (.not.useovrlp(jmo,kmo,is)) cycle

                f3temp = 0d0
                do iat_a = iat_a_ini, iat_a_end
                   do iat_b = iat_b_ini, iat_b_end
                      do rat_b = rat_b_ini, rat_b_end
                         do iat_c = 1, nattr
                            do rat_c = 1, nlattot
                               f3temp(iat_a,iat_b,rat_b,iat_c,rat_c) = &
                                  real(res%sijc(imo,kmo,iat_a,is) * &
                                  res%sijc(res%sij_wnr_imap(jmo,rat_b),res%sij_wnr_imap(imo,rat_b),iat_b,is) * &
                                  res%sijc(res%sij_wnr_imap(kmo,rat_c),res%sij_wnr_imap(jmo,rat_c),iat_c,is),8) + &
                                  real(res%sijc(imo,kmo,iat_a,is) * &
                                  res%sijc(res%sij_wnr_imap(kmo,rat_b),res%sij_wnr_imap(jmo,rat_b),iat_b,is) * &
                                  res%sijc(res%sij_wnr_imap(jmo,rat_c),res%sij_wnr_imap(imo,rat_c),iat_c,is),8)
                            end do ! rat_c
                         end do ! iat_c
                      end do ! rat_b
                   end do ! iat_b
                end do ! iat_a

                !$omp critical (addfa)
                res%fa3(:,:,:,:,:,is) = res%fa3(:,:,:,:,:,is) + f3temp
                !$omp end critical (addfa)
             end do ! kmo
          end do ! jmo
       end do ! imo
       !$omp end parallel do
    end do ! is
    res%fa3 = 0.5d0 * res%fa3

    write (*,*) "xx checking individual values of DI3"
    do is = 1, nspin
       do iat_a = iat_a_ini, iat_a_end
          do iat_b = iat_b_ini, iat_b_end
             do rat_b = rat_b_ini, rat_b_end
                do iat_c = 1, nattr
                   do rat_c = 1, nlattot
                      write (*,'("is=",I1," A:",I2," B:",I2,"/",I2," C:",I2,"/",I2," = ",F15.8)') &
                         is, iat_a, iat_b, rat_b, iat_c, rat_c, res%fa3(iat_a,iat_b,rat_b,iat_c,rat_c,is)
                   end do
                end do
             end do
          end do
       end do
    end do
    write (*,*) "fin"

    ! write (*,*) "xx checking consistency between fa and fa3"
    ! do is = 1, nspin
    !    do iat_a = iat_a_ini, iat_a_end
    !       do iat_b = iat_b_ini, iat_b_end
    !          do rat_b = rat_b_ini, rat_b_end
    !             write (*,'(I3," ",I3," ",I3," ",4(F20.12," "))') iat_a, iat_b, rat_b,&
    !                abs(res%fa(iat_a,iat_b,rat_b,is)-sum(res%fa3(iat_a,iat_b,rat_b,:,:,is))),&
    !                abs(res%fa(iat_a,iat_b,rat_b,is)/sum(res%fa3(iat_a,iat_b,rat_b,:,:,is))),&
    !                abs(res%fa(iat_a,iat_b,rat_b,is)),sum(res%fa3(iat_a,iat_b,rat_b,:,:,is))
    !          end do
    !       end do
    !    end do
    ! end do
    ! stop 1

  end subroutine calc_di3_wannier

  !> Calculate the Fa matrix from the Sij matrix, psink version
  subroutine calc_fa_psink(res,nmo,nbnd,nlat,nattr,nspin,kpt,occ)
    use types, only: int_result
    use param, only: tpi, img
    type(int_result), intent(inout) :: res
    integer, intent(in) :: nmo, nbnd, nlat(3), nattr, nspin
    real*8, intent(in) :: kpt(:,:), occ(:,:,:)

    integer :: i, j, is, k, ia, ja, ka, imo1, ik1, ibnd1, imo2, ik2, ibnd2
    real*8 :: fatemp, kdif(3), fspin
    integer :: nlattot

    ! set the spin multiplier
    if (nspin == 1) then
       fspin = 2d0
    else
       fspin = 1d0
    end if

    nlattot = nlat(1)*nlat(2)*nlat(3)
    allocate(res%fa(nattr,nattr,nlattot,nspin))
    res%fa = 0d0
    do is = 1, nspin
       !$omp parallel do private(k,fatemp,imo1,ibnd1,ik1,imo2,ibnd2,ik2,kdif)
       do i = 1, nattr
          do j = 1, nattr

             k = 0
             do ia = 0, nlat(1)-1
                do ja = 0, nlat(2)-1
                   do ka = 0, nlat(3)-1
                      k = k + 1
                      fatemp = 0d0

                      do imo1 = 1, nmo
                         ibnd1 = modulo(imo1-1,nbnd) + 1
                         ik1 = (imo1-1) / nbnd + 1
                         do imo2 = 1, nmo
                            ibnd2 = modulo(imo2-1,nbnd) + 1
                            ik2 = (imo2-1) / nbnd + 1

                            kdif = kpt(:,ik2) - kpt(:,ik1)
                            res%fa(i,j,k,is) = res%fa(i,j,k,is) + occ(ibnd1,ik1,is) * occ(ibnd2,ik2,is) * &
                               real(res%sijc(imo2,imo1,i,is) * res%sijc(imo1,imo2,j,is) * &
                               exp(tpi*img*(kdif(1)*ia+kdif(2)*ja+kdif(3)*ka)),8)
                         end do ! imo2 = imo1, nmo
                      end do ! imo1 = 1, nmo
                   end do ! ka = 0, nlat(3)-1
                end do ! ja = 0, nlat(2)-1
             end do ! ia = 0, nlat(1)-1
          end do ! j = 1, nattr
       end do ! i = 1, nattr
       !$omp end parallel do
    end do ! is = 1, nspin

    res%fa = res%fa / (fspin*fspin)

  end subroutine calc_fa_psink

  !> Calculate the Fa matrix from the Sij matrix, psink version
  subroutine calc_di3_psink(res,nmo,nbnd,nlat,nattr,nspin,kpt,occ,atom1,atom2)
    use types, only: int_result
    use param, only: tpi, img
    type(int_result), intent(inout) :: res
    integer, intent(in) :: nmo, nbnd, nlat(3), nattr, nspin
    real*8, intent(in) :: kpt(:,:), occ(:,:,:)
    integer, intent(in) :: atom1, atom2(4)

    integer :: i, j, is, k, n
    real*8 :: kdifb(3), kdifc(3), fspin, pfac
    integer :: nlattot
    integer :: imo, jmo, kmo
    integer :: ibnd, jbnd, kbnd, ik, jk, kk
    integer :: iat_a, iat_b, iat_c, rat_b, rat_c
    integer :: iat_a_ini, iat_a_end, iat_b_ini, iat_b_end, rat_b_ini, rat_b_end

    real*8, allocatable :: f3temp(:,:,:,:,:)
    logical, allocatable :: useovrlp(:,:,:)
    real*8, parameter :: ovrlp_thr = 1d-50

    real*8, allocatable :: lvec(:,:)

    ! number of lattice vectors
    nlattot = nlat(1)*nlat(2)*nlat(3)

    ! set limits for the sums
    iat_a_ini = 1
    iat_a_end = nattr
    if (atom1 > 0 .and. atom1 <= nattr) then
       iat_a_ini = atom1
       iat_a_end = atom1
    end if
    iat_b_ini = 1
    iat_b_end = nattr
    rat_b_ini = 1
    rat_b_end = nlattot
    if (atom2(1) > 0 .and. atom2(1) <= nattr) then
       iat_b_ini = atom2(1)
       iat_b_end = atom2(1)
       n = 0
       loopi: do i = 0, nlat(1)-1
          do j = 0, nlat(2)-1
             do k = 0, nlat(3)-1
                n = n + 1
                if (i == modulo(atom2(2),nlat(1)) .and. j == modulo(atom2(3),nlat(2)) .and.&
                   k == modulo(atom2(4),nlat(3))) then
                   rat_b_ini = n
                   rat_b_end = n
                   exit loopi
                end if
             end do
          end do
       end do loopi
    end if

    ! list of lattice vectors
    allocate(lvec(3,nlattot))
    nlattot = 0
    do i = 0, nlat(1)-1
       do j = 0, nlat(2)-1
          do k = 0, nlat(3)-1
             nlattot = nlattot + 1
             lvec(1,nlattot) = real(i,8)
             lvec(2,nlattot) = real(j,8)
             lvec(3,nlattot) = real(k,8)
          end do
       end do
    end do

    ! flag which overlaps will be used
    allocate(useovrlp(nmo,nmo,nspin))
    useovrlp = .true.
    do is = 1, nspin
       do imo = 1, nmo
          do jmo = imo, nmo
             useovrlp(imo,jmo,is) = (maxval(abs(res%sijc(imo,jmo,:,is))) >= ovrlp_thr)
          end do
       end do
    end do

    ! three-atom index calculation
    allocate(res%fa3(nattr,nattr,nlattot,nattr,nlattot,nspin))
    allocate(f3temp(nattr,nattr,nlattot,nattr,nlattot))
    res%fa3 = 0d0
    f3temp = 0d0
    do is = 1, nspin
       do imo = 1, nmo
          ibnd = modulo(imo-1,nbnd) + 1
          ik = (imo-1) / nbnd + 1
          do jmo = 1, nmo
             if (.not.useovrlp(imo,jmo,is)) cycle
             jbnd = modulo(jmo-1,nbnd) + 1
             jk = (jmo-1) / nbnd + 1
             do kmo = 1, nmo
                if (.not.useovrlp(imo,kmo,is)) cycle
                if (.not.useovrlp(jmo,kmo,is)) cycle
                kbnd = modulo(kmo-1,nbnd) + 1
                kk = (kmo-1) / nbnd + 1


                f3temp = 0d0
                do iat_a = iat_a_ini, iat_a_end
                   do iat_b = iat_b_ini, iat_b_end
                      do rat_b = rat_b_ini, rat_b_end
                         do iat_c = 1, nattr
                            do rat_c = 1, nlattot
                               kdifb = kpt(:,ik) - kpt(:,jk)
                               kdifc = kpt(:,jk) - kpt(:,kk)
                               pfac = kdifb(1)*lvec(1,rat_b)+kdifb(2)*lvec(2,rat_b)+kdifb(3)*lvec(3,rat_b) +&
                                  kdifc(1)*lvec(1,rat_c)+kdifc(2)*lvec(2,rat_c)+kdifc(3)*lvec(3,rat_c)
                               f3temp(iat_a,iat_b,rat_b,iat_c,rat_c) = f3temp(iat_a,iat_b,rat_b,iat_c,rat_c) +&
                                  real(res%sijc(imo,kmo,iat_a,is)*res%sijc(jmo,imo,iat_b,is)*res%sijc(kmo,jmo,iat_c,is)*&
                                  exp(tpi*img*pfac),8)

                               kdifb = kpt(:,jk) - kpt(:,kk)
                               kdifc = kpt(:,jk) - kpt(:,ik)
                               pfac = kdifb(1)*lvec(1,rat_b)+kdifb(2)*lvec(2,rat_b)+kdifb(3)*lvec(3,rat_b) +&
                                  kdifc(1)*lvec(1,rat_c)+kdifc(2)*lvec(2,rat_c)+kdifc(3)*lvec(3,rat_c)
                               f3temp(iat_a,iat_b,rat_b,iat_c,rat_c) = f3temp(iat_a,iat_b,rat_b,iat_c,rat_c) +&
                                  real(res%sijc(imo,kmo,iat_a,is)*res%sijc(kmo,jmo,iat_b,is)*res%sijc(jmo,imo,iat_c,is)*&
                                  exp(tpi*img*pfac),8)
                            end do ! rat_c
                         end do ! iat_c
                      end do ! rat_b
                   end do ! iat_b
                end do ! iat_a
                f3temp = f3temp * occ(ibnd,ik,is) * occ(jbnd,jk,is) * occ(kbnd,kk,is)

                !$omp critical (addfa)
                res%fa3(:,:,:,:,:,is) = res%fa3(:,:,:,:,:,is) + f3temp
                !$omp end critical (addfa)
             end do ! kmo = 1, nmo
          end do ! jmo = 1, nmo
       end do ! imo = 1, nmo
    end do ! is = 1, nspin

    ! set the spin multiplier
    if (nspin == 1) then
       fspin = 2d0
    else
       fspin = 1d0
    end if
    res%fa3 = 0.5d0 * res%fa3 / (fspin*fspin*fspin)

    ! do is = 1, nspin
    !    do iat_a = iat_a_ini, iat_a_end
    !       do iat_b = iat_b_ini, iat_b_end
    !          do rat_b = rat_b_ini, rat_b_end
    !             do iat_c = 1, nattr
    !                do rat_c = 1, nlattot
    !                   write (*,'("is=",I1," A:",I2," B:",I2,"/",I2," C:",I2,"/",I2," = ",F15.8)') &
    !                      is, iat_a, iat_b, rat_b, iat_c, rat_c, res%fa3(iat_a,iat_b,rat_b,iat_c,rat_c,is)
    !                end do
    !             end do
    !          end do
    !       end do
    !    end do
    ! end do
    ! write (*,*) "fin"

    ! write (*,*) "xx checking consistency between fa and fa3"
    ! do is = 1, nspin
    !    do iat_a = iat_a_ini, iat_a_end
    !       do iat_b = iat_b_ini, iat_b_end
    !          do rat_b = rat_b_ini, rat_b_end
    !             write (*,'(I3," ",I3," ",I3," ",4(F20.12," "))') iat_a, iat_b, rat_b,&
    !                abs(res%fa(iat_a,iat_b,rat_b,is)-sum(res%fa3(iat_a,iat_b,rat_b,:,:,is))),&
    !                abs(res%fa(iat_a,iat_b,rat_b,is)/sum(res%fa3(iat_a,iat_b,rat_b,:,:,is))),&
    !                abs(res%fa(iat_a,iat_b,rat_b,is)),sum(res%fa3(iat_a,iat_b,rat_b,:,:,is))
    !          end do
    !       end do
    !    end do
    ! end do
    ! stop 1

  end subroutine calc_di3_psink

  ! Calculate the information necessary to relate Sij^{A+R} to Sij^A.
  ! This information is different if Sij is calculated in terms of
  ! wannier functions (i,j=nR) or bloch functions (i,j=nk).
  subroutine find_sij_translations(res,nmo,nbnd,nlat,nlattot)
    use types, only: int_result, sijtype_wnr
    type(int_result), intent(inout) :: res
    integer, intent(in) :: nmo, nbnd, nlat(3), nlattot

    integer :: imo, ia, ja, ka, iba, k, ic, jc, kc, idx(3)

    if (res%sijtype /= sijtype_wnr) return

    ! translations for Sij from wannier functions

    ! Wannier functions:
    !     S_{nR,n'R'}^{A+R''} = S_{nR-R'',n'R'-R''}^{A}
    ! !! S(imo,1) translated by lattice vector k is S(imap(imo,k),1) !!
    if (allocated(res%sij_wnr_imap)) deallocate(res%sij_wnr_imap)
    allocate(res%sij_wnr_imap(nmo,nlattot))
    do imo = 1, nmo
       call unpackidx(imo,ia,ja,ka,iba,nmo,nbnd,nlat)
       k = 0
       do ic = 0, nlat(1)-1
          do jc = 0, nlat(2)-1
             do kc = 0, nlat(3)-1
                k = k + 1
                idx = (/ia-ic, ja-jc, ka-kc/)
                idx = modulo(idx,nlat)
                call packidx(idx(1),idx(2),idx(3),iba,res%sij_wnr_imap(imo,k),nmo,nbnd,nlat)
             end do
          end do
       end do
    end do

  end subroutine find_sij_translations

  ! Check the sanity of the sij matrix using the orthogonality and
  ! number of electron relations.
  subroutine check_sij_sanity(res,qe,nspin,nmo,nbnd,nlat,nlattot)
    use types, only: int_result, sijtype_wnr
    use grid3mod, only: qedat
    use param, only: tpi, img
    type(int_result), intent(in) :: res
    type(qedat), intent(in) :: qe
    integer, intent(in) :: nspin, nmo, nbnd, nlat(3), nlattot

    integer :: is, imo, jmo, k, l1, l2, l3
    integer :: ik1, ibnd1, ik2, ibnd2
    real*8 :: fspin, kdif(3), delta(3)
    complex*16 :: asum, saux

    ! set the spin multiplier
    if (nspin == 1) then
       fspin = 2d0
    else
       fspin = 1d0
    end if

    if (res%sijtype == sijtype_wnr) then
       ! checking sij from wannier functions

       ! sum_AR Sij^(A+R) = delta_ij
       do is = 1, nspin
          do imo = 1, nmo
             do jmo = 1, nmo
                asum = 0d0
                do k = 1, nlattot
                   asum = asum + sum(res%sijc(res%sij_wnr_imap(imo,k),res%sij_wnr_imap(jmo,k),:,is))
                end do
                write (*,*) is, imo, jmo, asum
             end do
          end do
       end do

       ! sum_i sum_A Sii^A = N
       do is = 1, nspin
          asum = 0d0
          do imo = 1, nmo
             asum = asum + sum(res%sijc(imo,imo,:,is))
          end do
          write (*,*) is, asum * fspin
       end do
    else
       ! sum_AR Sij^(A+R) = delta_ij
       do is = 1, nspin
          imo = 0
          do ik1 = 1, nlattot
             do ibnd1 = 1, nbnd
                imo = imo + 1

                jmo = 0
                do ik2 = 1, nlattot
                   do ibnd2 = 1, nbnd
                      jmo = jmo + 1

                      kdif = qe%kpt(:,ik2) - qe%kpt(:,ik1)

                      saux = 0d0
                      do l1 = 1, nlat(1)
                         do l2 = 1, nlat(2)
                            do l3 = 1, nlat(3)
                               delta = real((/l1,l2,l3/)-1,8)
                               saux = saux + exp(tpi*img*(kdif(1)*l1+kdif(2)*l2+kdif(3)*l3))
                            end do
                         end do
                      end do
                      saux = saux

                      asum = sum(res%sijc(imo,jmo,:,is)) * saux / real(nlat(1)*nlat(2)*nlat(3),8)
                      write (*,*) is, imo, jmo, asum
                   end do
                end do
             end do
          end do
       end do

       ! sum_i sum_A Sii^A = N
       do is = 1, nspin
          asum = 0d0
          imo = 0
          do ik1 = 1, nlattot
             do ibnd1 = 1, nbnd
                imo = imo + 1
                asum = asum + qe%occ(ibnd1,ik1,is) * sum(res%sijc(imo,imo,:,is))
             end do
          end do
          write (*,*) is, asum
       end do
    end if

  end subroutine check_sij_sanity

  !> Dummy function for quadpack integration
  function quadpack_f(x,unit,xnuc) result(res)
    use systemmod, only: sy

    real*8, intent(in) :: x
    real*8, intent(in) :: unit(3)
    real*8, intent(in) :: xnuc(3)
    real*8 :: res(sy%npropi)

    real*8 :: xaux(3)

    xaux = xnuc + x * unit
    call sy%grdall(xaux,res)
    res = res * x * x

  end function quadpack_f

  !> Write calculated multipole moments to the output. bas =
  !> integration driver data, res(1:npropi) = results.
  subroutine int_output_multipoles(bas,res)
    use systemmod, only: sy
    use tools_io, only: uout, string, ioj_left, ioj_center
    use types, only: basindat, int_result, out_mpoles
    type(basindat), intent(in) :: bas
    type(int_result), intent(in) :: res(:)

    integer :: i, j, k, l, nn, fid, lmax
    ! integer :: i, j, k, nn, n, l, lmax, fid
    integer, allocatable :: l_(:), m_(:)
    character(len=:), allocatable :: lbl
    character*3 :: ls, ms
    character(len=:), allocatable :: sncp, scp, sname, sz, smult

    if (bas%imtype == imtype_isosurface) return

    do l = 1, sy%npropi
       if (.not.res(l)%done .or. res(l)%outmode /= out_mpoles) cycle
       fid = sy%propi(l)%fid
       write (uout,'("* Basin multipole moments (using real solid harmonics)")')
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(sy%propi(l)%prop_name)
       write (uout,'("+ The calculated multipoles are: ")')
       write (uout,'("    Q_lm^A = int_A rho(r) * Rlm(r) dr ")')
       write (uout,'("  where the integral is over the basin of A, and Rlm is a real solid harmonic.")')
       write (uout,'("  The coordinates are referred to the attractor of the A basin. All quantities")')
       write (uout,'("  in atomic units.")')
       write (uout,*)

       ! figure out the indices for the labels
       lmax = sy%propi(l)%lmax
       allocate(l_((lmax+1)*(lmax+1)),m_((lmax+1)*(lmax+1)))
       nn = 0
       do i = 0, lmax
          do j = -i, i
             nn = nn + 1
             l_(nn) = i
             m_(nn) = j
          end do
       end do

       ! pretty output
       nn = (lmax+1)**2/5
       if (mod((lmax+1)**2,5) /= 0) nn = nn + 1
       do i = 1, nn
          ! header
          lbl = "# Id   cp   ncp   Name  Z   "
          do j = (i-1)*5+1,min(5*i,(lmax+1)**2)
             if (j == 1) then
                lbl = lbl // " " // string("1",15,justify=ioj_center)
             elseif (j == 2) then
                lbl = lbl // " " // string("x",15,justify=ioj_center)
             elseif (j == 3) then
                lbl = lbl // " " // string("z",15,justify=ioj_center)
             elseif (j == 4) then
                lbl = lbl // " " // string("y",15,justify=ioj_center)
             elseif (j == 5) then
                lbl = lbl // " " // string("sq(3)/2(x^2-y^2)",15,justify=ioj_center)
             elseif (j == 6) then
                lbl = lbl // " " // string("sq(3)xz",15,justify=ioj_center)
             elseif (j == 7) then
                lbl = lbl // " " // string("(3z^2-r^2)/2",15,justify=ioj_center)
             elseif (j == 8) then
                lbl = lbl // " " // string("sq(3)yz",15,justify=ioj_center)
             elseif (j == 9) then
                lbl = lbl // " " // string("sq(3)xy",15,justify=ioj_center)
             else
                write (ls,'(I3)') l_(j)
                write (ms,'(I3)') abs(m_(j))
                if (m_(j) <= 0) then
                   lbl = lbl // " " // string("C("//string(ls)//","//string(ms)//")",15,justify=ioj_center)
                else
                   lbl = lbl // " " // string("S("//string(ls)//","//string(ms)//")",15,justify=ioj_center)
                end if
             endif
          end do
          write (uout,'(A)') lbl

          ! body
          do j = 1, bas%nattr
             if (.not.bas%docelatom(bas%icp(i))) cycle
             call assign_strings(j,bas%icp(j),.false.,scp,sncp,sname,smult,sz)
             write (uout,'("  ",99(A," "))') &
                string(j,4,ioj_left), scp, sncp, sname, sz, &
                (string(res(l)%mpole((i-1)*5+1+k,j),'e',15,8,4),k=0,min(4,size(res(l)%mpole,1)-(i-1)*5-1))
          enddo
          if (i < nn) then
             write (uout,*)
          else
             write (uout,'(32("-"),99(A))') ("----------------",j=1,5)
          end if
       end do
       deallocate(l_,m_)
       write (uout,*)
    end do

  end subroutine int_output_multipoles

  !> Output the delocalization indices.  bas = integration driver
  !> data, res(1:npropi) = results.
  subroutine int_output_deloc(bas,res)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use global, only: iunit, iunitname0, dunit0
    use tools, only: qcksort
    use tools_io, only: uout, string, ioj_left, ioj_right, fopen_read,&
       fopen_write, fclose
    use types, only: basindat, int_result, out_deloc
    type(basindat), intent(in) :: bas
    type(int_result), intent(in) :: res(:)

    integer :: i, j, k, l, m, n, jo, ko, kk
    integer :: fid, natt, nlat(3), nspin, nlattot
    real*8 :: fspin, xli, xnn, r1(3), r2(3), asum, d2, raux, sfac
    real*8, allocatable :: dimol(:,:,:,:,:), limol(:), namol(:)
    real*8, allocatable :: dist(:), diout(:), xcm(:,:), diout_i(:), diout_o(:)
    logical, allocatable :: diout_ioprint1(:), diout_ioprint2(:)
    integer, allocatable :: io(:), ilvec(:,:), ilat(:), idat(:), idxmol(:,:)
    integer :: ia, ja
    integer :: ic, jc, kc, lvec1(3), lvec2(3), lvec3(3)
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    type(crystal) :: cr1
    type(crystalseed) :: ncseed
    logical :: di3, ok
    real*8 :: di3aux, daux

    if (bas%imtype == imtype_isosurface) return

    do l = 1, sy%npropi
       if (.not.res(l)%done) cycle
       if (res(l)%outmode /= out_deloc) cycle

       write (uout,'("* Localization and delocalization indices")')

       ! header
       fid = sy%propi(l)%fid
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(sy%propi(l)%prop_name)

       ! some integers for the run
       nlat = res(l)%nlat
       nspin = res(l)%nspin
       nlattot = nlat(1)*nlat(2)*nlat(3)
       natt = bas%nattr

       ! spin factor
       if (nspin == 1) then
          fspin = 2d0
       else
          fspin = 1d0
       end if

       ! localization indices
       write (uout,'("+ Localization (LI) and total population (N)")')
       write (uout,'("# Id   cp   ncp   Name  Z       LI(A)           N(A)")')
       do i = 1, natt
          if (.not.bas%docelatom(bas%icp(i))) cycle
          call assign_strings(i,bas%icp(i),.false.,scp,sncp,sname,smult,sz)
          xli = sum(abs(res(l)%fa(i,i,1,:))) * fspin
          xnn = sum(abs(res(l)%fa(i,:,:,:))) * fspin
          write (uout,'("  ",99(A," "))') string(i,4,ioj_left), scp, sncp, sname, sz, &
             string(xli,'f',15,8,4), string(xnn,'f',12,8,4)
       end do
       write (uout,*)

       ! build the supercell
       ncseed%isused = .true.
       do i = 1, 3
          ncseed%m_x2c(:,i) = sy%c%m_x2c(:,i) * nlat(i)
       end do
       ncseed%useabr = 2
       ncseed%nat = 0
       ncseed%havesym = 0
       ncseed%findsym = 0
       ncseed%ismolecule = sy%c%ismolecule
       call cr1%struct_new(ncseed,.true.)

       ! header
       di3 = sy%propi(l)%di3
       write (uout,'("+ Delocalization indices")')
       write (uout,'("  Each block gives information about a single atom in the main cell.")')
       write (uout,'("  First line: localization index. Next lines: delocalization index")')
       write (uout,'("  with all atoms in the environment. Last line: sum of LI + 0.5 * DIs,")')
       write (uout,'("  equal to the atomic population. Distances are in ",A,".")') iunitname0(iunit)
       if (di3) then
          write (uout,'("  Three-center: 3c(in) contribution to DI(AB) = sum of N(ABC) with C=A,B.")')
          write (uout,'("  3c(in) contribution to LI(A) = N(AAA). 3c(out) contribution = the DI or LI")')
          write (uout,'("  minus the 3c(in).")')
       end if
       write (uout,*)
       allocate(dist(natt*nlattot),io(natt*nlattot),diout(natt*nlattot),ilvec(3,natt*nlattot),idat(natt*nlattot))

       ! prepare for the di3 output
       ! diout_ioprint1 = print inner and outer contribution to DI3, filtered with di3_atom1
       ! diout_ioprint2 = same, but filtered with di3_atom2
       if (di3) then
          allocate(diout_i(natt*nlattot),diout_o(natt*nlattot),diout_ioprint1(natt),diout_ioprint2(natt*nlattot))
          k = 0
          m = 0
          do ic = 0, nlat(1)-1
             do jc = 0, nlat(2)-1
                do kc = 0, nlat(3)-1
                   k = k + 1
                   do j = 1, natt
                      m = m + 1
                      ok = sy%propi(l)%di3_atom2(1) < 0 .or. sy%propi(l)%di3_atom2(1) > natt .or.&
                         (sy%propi(l)%di3_atom2(1) == j) .and.&
                         ic == modulo(sy%propi(l)%di3_atom2(2),nlat(1)) .and.&
                         jc == modulo(sy%propi(l)%di3_atom2(3),nlat(2)) .and.&
                         kc == modulo(sy%propi(l)%di3_atom2(4),nlat(3))
                      diout_ioprint2(m) = ok
                   end do
                end do
             end do
          end do
          if (sy%propi(l)%di3_atom1 < 0 .or. sy%propi(l)%di3_atom1 > natt) then
             diout_ioprint1 = .true.
          else
             diout_ioprint1 = .false.
             diout_ioprint1(sy%propi(l)%di3_atom1) = .true.
          end if
       end if

       ! two-body indices table
       do i = 1, natt
          if (.not.bas%docelatom(bas%icp(i))) cycle
          call assign_strings(i,bas%icp(i),.false.,scp,sncp,sname,smult,sz)
          write (uout,'("# Attractor ",A," (cp=",A,", ncp=",A,", name=",A,", Z=",A,") at: ",3(A,"  "))') &
             string(i), trim(scp), trim(sncp), trim(adjustl(sname)), trim(sz), (trim(string(bas%xattr(j,i),'f',12,7)),j=1,3)
          ok = di3
          if (ok) ok = diout_ioprint1(i)
          if (ok) then
             write (uout,'("# Id   cp   ncp   Name  Z    Latt. vec.     &
               &----  Cryst. coordinates ----       Distance        LI/DI        3c(in)       3c(out)")')
          else
             write (uout,'("# Id   cp   ncp   Name  Z    Latt. vec.     &
               &----  Cryst. coordinates ----       Distance        LI/DI")')
          end if

          ! Precompute the localization/delocalization indices for this atom and
          ! calculate location and distance information. Obtain inner and outer
          ! contributions to DI3, if requested.
          dist = 0d0
          k = 0
          m = 0
          do ic = 0, nlat(1)-1
             do jc = 0, nlat(2)-1
                do kc = 0, nlat(3)-1
                   k = k + 1
                   do j = 1, natt
                      m = m + 1
                      io(m) = m
                      r1 = (bas%xattr(:,j) + (/ic,jc,kc/) - bas%xattr(:,i)) / real(nlat,8)
                      call cr1%shortest(r1,d2)
                      dist(m) = d2 * dunit0(iunit)
                      sfac = 1d0
                      if (dist(m) < 1d-5) sfac = 0.5d0

                      diout(m) = 2d0 * sum(res(l)%fa(i,j,k,:)) * fspin * sfac
                      if (di3) then
                         if (diout_ioprint1(i).and.diout_ioprint2(m)) then
                            if (j == i .and. k == 1) then
                               diout_i(m) = sum(res(l)%fa3(i,j,k,i,1,:))
                            else
                               diout_i(m) = sum(res(l)%fa3(i,j,k,i,1,:)) + sum(res(l)%fa3(i,j,k,j,k,:))
                            end if
                            diout_o(m) = sum(res(l)%fa3(i,j,k,:,:,:)) - diout_i(m)
                            diout_i(m) = diout_i(m) * 2d0 * fspin * sfac
                            diout_o(m) = diout_o(m) * 2d0 * fspin * sfac
                         end if
                      end if
                      idat(m) = j
                      ilvec(:,m) = nint(bas%xattr(:,i) + cr1%c2x(r1) * nlat - bas%xattr(:,j))
                   end do
                end do
             end do
          end do

          ! sort by increasing distance and output for this atom
          call qcksort(dist,io,1,natt*nlattot)

          asum = 0d0
          do m = 1, natt*nlattot
             j = io(m)
             if (.not.bas%docelatom(bas%icp(idat(j)))) cycle
             if (dist(j) < 1d-5) then
                ok = di3
                if (ok) ok = diout_ioprint1(i).and.diout_ioprint2(j)
                if (ok) then
                   write (uout,'("  Localization index",71("."),3(A," "))') string(diout(j),'f',12,8,4),&
                      string(diout_i(j),'f',12,8,4), string(diout_o(j),'f',12,8,4)
                else
                   write (uout,'("  Localization index",71("."),A)') string(diout(j),'f',12,8,4)
                end if
                asum = asum + diout(j)
             else
                call assign_strings(j,bas%icp(idat(j)),.false.,scp,sncp,sname,smult,sz)
                r1 = bas%xattr(:,idat(j)) + ilvec(:,j)
                ok = di3
                if (ok) ok = diout_ioprint1(i).and.diout_ioprint2(j)
                if (ok) then
                   write (uout,'("  ",99(A," "))') string(j,4,ioj_left), scp, sncp, sname, sz,&
                      (string(ilvec(k,j),3,ioj_right),k=1,3), (string(r1(k),'f',12,7,4),k=1,3),&
                      string(dist(j),'f',12,7,4), string(diout(j),'f',12,8,4),&
                      string(diout_i(j),'f',12,8,4), string(diout_o(j),'f',12,8,4)
                else
                   write (uout,'("  ",99(A," "))') string(j,4,ioj_left), scp, sncp, sname, sz,&
                      (string(ilvec(k,j),3,ioj_right),k=1,3), (string(r1(k),'f',12,7,4),k=1,3),&
                      string(dist(j),'f',12,7,4), string(diout(j),'f',12,8,4)
                end if
                asum = asum + 0.5d0 * diout(j)
             end if
          end do
          write (uout,'("  Total (atomic population)",64("."),A)') string(asum,'f',12,8,4)
          write (uout,*)
       end do

       ! three-body indices
       if (di3) then
          allocate(ilat(natt*nlattot))
          write (uout,'("* Three-center delocalization indices")')
          write (uout,'("  Each block shows the three-center contribution N(ABC) to each calculated")')
          write (uout,'("  localization and delocalization indices, sorted by distance to atom A.")')
          write (uout,'("  Distances are in ",A,".")') iunitname0(iunit)
          write (uout,*)

          do i = 1, natt
             if (.not.bas%docelatom(bas%icp(i))) cycle
             if (.not.diout_ioprint1(i)) cycle
             call assign_strings(i,bas%icp(i),.false.,scp,sncp,sname,smult,sz)
             write (uout,'("+ Three-center decomposition of delocalization indices related to:")')
             write (uout,'("+ Attractor ",A," (cp=",A,", ncp=",A,", name=",A,", Z=",A,") at: (",2(A,","),A,")")') &
                string(i), trim(scp), trim(sncp), trim(adjustl(sname)), trim(sz),&
                (trim(string(bas%xattr(k,i),'f',decimal=7)),k=1,3)
             write (uout,*)

             ! Sort the second atom by distance to i
             ! calculate location and distance information. Obtain inner and outer
             ! contributions to DI3, if requested.
             dist = 0d0
             k = 0
             m = 0
             do ic = 0, nlat(1)-1
                do jc = 0, nlat(2)-1
                   do kc = 0, nlat(3)-1
                      k = k + 1
                      do j = 1, natt
                         m = m + 1
                         io(m) = m
                         r1 = (bas%xattr(:,j) + (/ic,jc,kc/) - bas%xattr(:,i)) / real(nlat,8)
                         call cr1%shortest(r1,d2)
                         dist(m) = d2 * dunit0(iunit)
                         idat(m) = j
                         ilvec(:,m) = nint(bas%xattr(:,i) + cr1%c2x(r1) * nlat - bas%xattr(:,j))
                         ilat(m) = k
                      end do
                   end do
                end do
             end do
             call qcksort(dist,io,1,natt*nlattot)

             ! run over the second atom
             do m = 1, natt*nlattot
                jo = io(m)
                if (.not.bas%docelatom(bas%icp(idat(jo)))) cycle
                if (.not.diout_ioprint2(m)) cycle

                ! write A again
                call assign_strings(i,bas%icp(i),.false.,scp,sncp,sname,smult,sz)
                write (uout,'("# A: Attractor ",A," (cp=",A,", ncp=",A,", name=",A,", Z=",A,") at: (",2(A,","),A,")")') &
                   string(i), trim(scp), trim(sncp), trim(adjustl(sname)), trim(sz),&
                   (trim(string(bas%xattr(k,i),'f',decimal=7)),k=1,3)

                ! write B
                call assign_strings(jo,bas%icp(idat(jo)),.false.,scp,sncp,sname,smult,sz)
                write (uout,'("# B: Attractor ",A," (cp=",A,", ncp=",A,", name=",A,", Z=",A,") &
                   &at: (",2(A,","),A,")")') &
                   string(idat(jo)), trim(scp), trim(sncp), trim(adjustl(sname)), trim(sz),&
                   (trim(string(bas%xattr(kk,idat(jo)) + ilvec(kk,jo),'f',decimal=7)),kk=1,3)
                write (uout,'("# C: each row in the following table.")')

                if (dist(jo) < 1d-5) then
                   sfac = 0.5d0
                   write (uout,'("# A same as B, LI =",A)') string(diout(jo),'f',decimal=8)
                else
                   sfac = 1d0
                   write (uout,'("# dist(AB) = ",A,", DI = ",A)') &
                      string(dist(jo),'f',decimal=7), string(diout(jo),'f',decimal=8)
                end if

                ! run over the third atom
                write (uout,'("# Id   cp   ncp   Name  Z    Latt. vec.     &
                   &----  Cryst. coordinates ----        dist(AB)     dist(AC)     N(ABC)")')
                do n = 1, natt*nlattot
                   ko = io(n)
                   if (.not.bas%docelatom(bas%icp(idat(ko)))) cycle
                   call assign_strings(ko,bas%icp(idat(ko)),.false.,scp,sncp,sname,smult,sz)
                   r1 = bas%xattr(:,idat(ko)) + ilvec(:,ko)
                   r2 = (bas%xattr(:,idat(jo)) + ilvec(:,jo) - r1) / real(nlat,8)
                   call cr1%shortest(r2,d2)
                   daux = d2 * dunit0(iunit)
                   di3aux = 2d0 * sum(res(l)%fa3(i,idat(jo),ilat(jo),idat(ko),ilat(ko),:)) * fspin * sfac
                   write (uout,'("  ",99(A," "))') string(idat(ko),4,ioj_left), scp, sncp, sname, sz,&
                      (string(ilvec(kk,ko),3,ioj_right),kk=1,3), (string(r1(kk),'f',12,7,4),kk=1,3),&
                      string(dist(ko),'f',12,7,4),&
                      string(daux,'f',12,7,4),&
                      string(di3aux,'f',12,8,4)
                end do
                asum = 2d0 * sum(res(l)%fa3(i,idat(jo),ilat(jo),:,:,:)) * fspin * sfac
                write (uout,'("  Total (2-center LI/DI)",80("."),A)') string(asum,'f',12,8,4)
                write (uout,*)
             end do ! m = 1, natt*nlattot
          end do ! i = 1, natt
       end if

       ! Integrated molecular LI/DI
       if (.not.sy%c%ismolecule .and. all(sy%c%mol(1:sy%c%nmol)%discrete)) then
          allocate(idxmol(2,natt))
          ! Assign attractors to molecules
          idxmol = 0
          do i = 1, natt
             jlo: do j = 1, sy%c%nmol
                do k = 1, sy%c%mol(j)%nat
                   if (bas%icp(i) == sy%c%mol(j)%at(k)%cidx) then
                      idxmol(1,i) = j
                      idxmol(2,i) = k
                      exit jlo
                   end if
                end do
             end do jlo
          end do

          ! assign DIs to molecules
          allocate(dimol(sy%c%nmol,sy%c%nmol,0:nlat(1)-1,0:nlat(2)-1,0:nlat(3)-1),limol(sy%c%nmol),namol(sy%c%nmol))
          dimol = 0d0
          limol = 0d0
          namol = 0d0
          do i = 1, natt
             if (.not.bas%docelatom(bas%icp(i))) cycle
             ia = idxmol(1,i)
             if (ia == 0) cycle
             limol(ia) = limol(ia) + sum(abs(res(l)%fa(i,i,1,:))) * fspin
             namol(ia) = namol(ia) + sum(abs(res(l)%fa(i,:,:,:))) * fspin
             lvec1 = sy%c%mol(ia)%at(idxmol(2,i))%lvec
             k = 0
             m = 0
             do ic = 0, nlat(1)-1
                do jc = 0, nlat(2)-1
                   do kc = 0, nlat(3)-1
                      k = k + 1
                      do j = 1, natt
                         if (.not.bas%docelatom(bas%icp(j))) cycle
                         m = m + 1
                         if (idxmol(1,j) == 0) cycle
                         if (i == j .and. k == 1) cycle
                         ja = idxmol(1,j)
                         if (ja == 0) cycle
                         lvec2 = sy%c%mol(ja)%at(idxmol(2,j))%lvec
                         lvec3 = lvec1 - lvec2 + (/ic,jc,kc/)
                         lvec3 = modulo(lvec3,nlat)
                         raux = 2d0 * sum(abs(res(l)%fa(i,j,k,:))) * fspin
                         if (ia == ja .and. all(lvec3 == 0)) then
                            limol(ia) = limol(ia) + 0.5d0 * raux
                         else
                            dimol(ia,ja,lvec3(1),lvec3(2),lvec3(3)) = &
                               dimol(ia,ja,lvec3(1),lvec3(2),lvec3(3)) + raux
                         endif
                      end do
                   end do
                end do
             end do
          end do

          ! localization indices
          write (uout,'("* Integrated molecular properties")')
          write (uout,'("+ Localization indices")')
          write (uout,'("# Mol       LI(A)           N(A)")')
          do i = 1, sy%c%nmol
             write (uout,'("  ",99(A," "))') &
                string(i,4,ioj_left), string(limol(i),'f',15,8,4), string(namol(i),'f',12,8,4)
          end do
          write (uout,*)

          ! centers of mass
          allocate(xcm(3,sy%c%nmol))
          do i = 1, sy%c%nmol
             xcm(:,i) = sy%c%mol(i)%cmass()
             xcm(:,i) = sy%c%c2x(xcm(:,i))
          end do

          write (uout,'("+ Delocalization indices")')
          do i = 1, sy%c%nmol
             ! name the molecule
             write (uout,'("# Molecule ",A," with ",A," atoms at ",3(A," "))') string(i), string(sy%c%mol(i)%nat),&
                (string(xcm(j,i),'f',10,6,3),j=1,3)
             write (uout,'("# Mol   Latt. vec.    ---- Center of mass (cryst) ----      Distance      LI/DI")')

             m = 0
             do ic = 0, nlat(1)-1
                do jc = 0, nlat(2)-1
                   do kc = 0, nlat(3)-1
                      do j = 1, sy%c%nmol
                         m = m + 1
                         io(m) = m
                         r1 = (xcm(:,j) + (/ic,jc,kc/) - xcm(:,i)) / real(nlat,8)
                         call cr1%shortest(r1,d2)
                         dist(m) = d2 * dunit0(iunit)
                         diout(m) = dimol(i,j,ic,jc,kc)
                         idat(m) = j
                         ilvec(:,m) = nint(xcm(:,i) + cr1%c2x(r1) * nlat - xcm(:,j))
                      end do
                   end do
                end do
             end do

             ! sort by increasing distance and output for this molecule
             call qcksort(dist,io,1,sy%c%nmol*nlattot)
             asum = 0d0
             do m = 1, sy%c%nmol*nlattot
                j = io(m)
                if (dist(j) < 1d-5) then
                   write (uout,'("  Localization index",51("."),A)') string(limol(i),'f',12,8,4)
                   asum = asum + limol(i)
                else
                   r1 = xcm(:,idat(j)) + ilvec(:,j)
                   write (uout,'("  ",99(A," "))') string(idat(j),4,ioj_left), &
                      (string(ilvec(k,j),3,ioj_right),k=1,3), (string(r1(k),'f',12,7,4),k=1,3),&
                      string(dist(j),'f',12,7,4), string(diout(j),'f',12,8,4)
                   asum = asum + 0.5d0 * diout(j)
                end if
             end do
             write (uout,'("  Total (atomic population)",44("."),A)') string(asum,'f',12,8,4)
             write (uout,*)
          end do
          deallocate(dimol,limol,namol,idxmol,xcm)
       end if

       ! clean up
       call cr1%end()
       deallocate(dist,io,diout,ilvec,idat)
       if (allocated(ilat)) deallocate(ilat)
       if (allocated(diout_i)) deallocate(diout_i)
       if (allocated(diout_o)) deallocate(diout_o)
       if (allocated(diout_ioprint1)) deallocate(diout_ioprint1)
       if (allocated(diout_ioprint2)) deallocate(diout_ioprint2)
    end do

  end subroutine int_output_deloc

  !> Output the Hirshfeld overlap populations. bas = integration
  !> driver data, res(1:npropi) = results.
  subroutine int_output_hirshfeld_overlap(bas,res)
    use global, only: iunit, iunitname0
    use tools, only: qcksort
    use tools_io, only: uout, string
    use tools_io, only: ioj_left, ioj_right
    use types, only: basindat, int_result, out_hirsh_ovpop
    type(basindat), intent(in) :: bas
    type(int_result), intent(in) :: res(:)

    integer :: i, j, k, l, m, n
    integer :: fid
    real*8, allocatable :: dist(:), ovout(:)
    integer, allocatable :: io(:), ilvec(:,:), idat(:)
    real*8 :: x1(3), x2(3), asum
    integer :: ic, jc, kc
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    logical :: hascycled

    real*8, parameter :: epsthr = 1d-8 ! threshold for showing population

    if (bas%imtype == imtype_isosurface) return

    do l = 1, sy%npropi
       if (.not.res(l)%done) cycle
       if (res(l)%outmode /= out_hirsh_ovpop) cycle

       write (uout,'("* Hirshfeld overlap populations (bond orders)")')

       ! header
       fid = sy%propi(l)%fid
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(sy%propi(l)%prop_name)

       ! header
       write (uout,'("  Each block gives information about a single atom in the main cell.")')
       write (uout,'("  Each line gives the overlap population (bond order) with all atoms")')
       write (uout,'("  in the environment that have a significant contribution (> ",A)') string(epsthr,'e',decimal=2)
       write (uout,'("  Last line: sum of all overlaps, must equal the Hirshfeld population.")')
       write (uout,'("  Distances are in ",A,".")') iunitname0(iunit)
       write (uout,*)

       n = bas%nattr * size(res(l)%hirsh_op,1) * size(res(l)%hirsh_op,2) * size(res(l)%hirsh_op,3)
       allocate(dist(n),io(n),ovout(n),ilvec(3,n),idat(n))

       ! two-body indices table
       do i = 1, bas%nattr
          if (.not.bas%docelatom(bas%icp(i))) cycle
          call assign_strings(i,bas%icp(i),.false.,scp,sncp,sname,smult,sz)
          write (uout,'("# Atom ",A," (cp=",A,", ncp=",A,", name=",A,", Z=",A,") at: ",3(A,"  "))') &
             string(i), trim(scp), trim(sncp), trim(adjustl(sname)), trim(sz), (trim(string(bas%xattr(j,i),'f',12,7)),j=1,3)
          write (uout,'("# Id   cp   ncp   Name  Z    Latt. vec.     &
             &----  Cryst. coordinates ----       Distance        Overlap")')

          ! Organize the data we will present in the table, for sorting
          n = 0
          do j = 1, bas%nattr
             do ic = lbound(res(l)%hirsh_op,3), ubound(res(l)%hirsh_op,3)
                do jc = lbound(res(l)%hirsh_op,4), ubound(res(l)%hirsh_op,4)
                   do kc = lbound(res(l)%hirsh_op,5), ubound(res(l)%hirsh_op,5)
                      if (res(l)%hirsh_op(i,j,ic,jc,kc) > epsthr) then
                         n = n + 1
                         x1 = bas%xattr(:,i)
                         x2 = bas%xattr(:,j) + (/ic,jc,kc/)
                         dist(n) = sy%c%distance(x1,x2)
                         io(n) = n
                         ovout(n) = res(l)%hirsh_op(i,j,ic,jc,kc)
                         idat(n) = j
                         ilvec(:,n) = (/ic,jc,kc/)
                      end if
                   end do
                end do
             end do
          end do

          ! sort by increasing distance and output for this atom
          call qcksort(dist,io,1,n)

          ! write the table
          hascycled = .false.
          asum = 0d0
          do m = 1, n
             j = io(m)
             if (i == bas%icp(idat(j)).and..not.sy%propi(l)%hirsh_show_self_overlaps) then
                hascycled = .true.
                cycle
             end if
             if (ovout(j) < sy%propi(l)%hirsh_cutoff) then
                hascycled = .true.
                cycle
             end if
             if (.not.bas%docelatom(bas%icp(idat(j)))) cycle
             call assign_strings(j,bas%icp(idat(j)),.false.,scp,sncp,sname,smult,sz)
             x1 = bas%xattr(:,idat(j)) + ilvec(:,j)
             write (uout,'("  ",99(A," "))') string(j,4,ioj_left), scp, sncp, sname, sz,&
                (string(ilvec(k,j),3,ioj_right),k=1,3), (string(x1(k),'f',12,7,4),k=1,3),&
                string(dist(j),'f',12,7,4), string(ovout(j),'f',12,8,4)
             asum = asum + ovout(j)
          end do
          if (.not.hascycled) &
             write (uout,'("  Total (Hirshfeld population)",61("."),A)') string(asum,'f',12,8,4)
          write (uout,*)
       end do ! i

       ! clean up
       deallocate(dist,io,ovout,ilvec,idat)
    end do ! l

  end subroutine int_output_hirshfeld_overlap

  !> Write a JSON file containing the structure, the reference
  !> field details, and the results of the YT/BADER integration.
  subroutine int_output_json(file,bas,res)
    use json_module, only: json_core, json_value
    use systemmod, only: sy, itype_v, itype_mpoles, itype_expr, itype_names
    use tools_io, only: uout, string, fopen_write, fclose, ferror, faterr
    use types, only: basindat, int_result, out_field
    character*(*), intent(in) :: file
    type(basindat), intent(in) :: bas
    type(int_result), intent(in) :: res(:)

    integer :: fid, nprop
    integer :: i, j, k, l
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    real*8 :: x(3), xcm(3)
    integer, allocatable :: idxmol(:,:), nacprop(:)
    real*8, allocatable :: sump(:), xres(:)
    type(json_core) :: json
    type(json_value), pointer :: p, o1, a1, o2, a2, o3

    if (bas%imtype == imtype_isosurface) return

    ! open and write some info about the structure and the reference field
    write (uout,'("* WRITE JSON file: ",A/)') string(file)

    call json%initialize()
    call json%create_object(p,'')
    call json%add(p,'units','bohr')
    call sy%c%struct_write_json(json,p)
    call sy%f(sy%iref)%write_json(json,p)

    ! the JSON report, beginning
    call json%create_object(o1,'integration')
    call json%add(p,o1)
    if (bas%imtype == imtype_yt) then
       call json%add(o1,'type','yt')
    elseif (bas%imtype == imtype_bader) then
       call json%add(o1,'type','bader')
    end if
    call json%add(o1,'grid_npoints',bas%n)

    ! list of properties
    nprop = count(res(1:sy%npropi)%done .and. res(1:sy%npropi)%outmode == out_field)
    allocate(nacprop(nprop))
    call json%add(o1,'number_of_properties',nprop)
    call json%create_array(a1,'properties')
    call json%add(o1,a1)
    l = 0
    do i = 1, sy%npropi
       fid = sy%propi(i)%fid
       if (.not.res(i)%done) cycle
       if (.not.res(i)%outmode == out_field) cycle
       l = l + 1
       nacprop(l) = i

       call json%create_object(o2,'')
       call json%add(a1,o2)
       call json%add(o2,'id',i)
       call json%add(o2,'label',trim(sy%propi(i)%prop_name))
       if (sy%propi(i)%itype == itype_expr) then
          call json%add(o2,'expression',trim(sy%propi(i)%expr))
       elseif (sy%propi(i)%itype == itype_mpoles) then
          call json%add(o2,'field_id',fid)
          call json%add(o2,'lmax',sy%propi(i)%lmax)
       elseif (sy%propi(i)%itype /= itype_v) then
          call json%add(o2,'field_id',fid)
       end if
       call json%add(o2,'field_name',trim(itype_names(sy%propi(i)%itype)))
       nullify(o2)
    end do
    nullify(a1)

    ! List of attractors
    call json%add(o1,'number_of_attractors',bas%nattr)
    call json%create_array(a1,'attractors')
    call json%add(o1,a1)
    allocate(xres(nprop))
    do i = 1, bas%nattr
       call assign_strings(i,bas%icp(i),.false.,scp,sncp,sname,smult,sz)
       x = bas%xattr(:,i)
       call json%create_object(o2,'')
       call json%add(a1,o2)
       call json%add(o2,'id',i)
       call json%add(o2,'cell_cp',trim(scp))
       call json%add(o2,'nonequivalent_cp',trim(sncp))
       call json%add(o2,'name',trim(adjustl(sname)))
       call json%add(o2,'atomic_number',trim(sz))
       call json%add(o2,'fractional_coordinates',x)

       xres = 0d0
       do j = 1, nprop
          xres(j) = res(nacprop(j))%psum(i)
       end do
       call json%add(o2,'integrals',xres)
       nullify(o2)
    end do
    nullify(a1)
    deallocate(xres)

    ! list of molecules and positions
    if (.not.sy%c%ismolecule .and. all(sy%c%mol(1:sy%c%nmol)%discrete)) then
       allocate(idxmol(2,bas%nattr))
       ! assign attractors to molecules
       idxmol = 0
       do i = 1, bas%nattr
          jlo: do j = 1, sy%c%nmol
             do k = 1, sy%c%mol(j)%nat
                if (bas%icp(i) == sy%c%mol(j)%at(k)%cidx) then
                   idxmol(1,i) = j
                   idxmol(2,i) = k
                   exit jlo
                end if
             end do
          end do jlo
       end do

       allocate(sump(nprop))
       do k = 1, sy%c%nmol
          do i = 1, bas%nattr
             if (idxmol(1,i) /= k) cycle
          end do
       end do

       ! list of molecules and associated attractors
       call json%add(o1,'number_of_molecules',sy%c%nmol)
       call json%create_array(a1,'molecules')
       call json%add(o1,a1)

       do i = 1, sy%c%nmol
          xcm = sy%c%mol(i)%cmass()
          xcm = sy%c%c2x(xcm)

          call json%create_object(o2,'')
          call json%add(a1,o2)

          call json%add(o2,'id',i)
          call json%add(o2,'center_of_mass',xcm)
          call json%add(o2,'number_of_atoms',sy%c%mol(i)%nat)

          call json%create_array(a2,'atoms')
          call json%add(o2,a2)
          l = 0
          sump = 0d0
          do j = 1, bas%nattr
             if (idxmol(1,j) /= i) cycle
             l = l + 1
             call assign_strings(j,bas%icp(j),.false.,scp,sncp,sname,smult,sz)
             x = sy%c%mol(idxmol(1,j))%at(idxmol(2,j))%x

             call json%create_object(o3,'')
             call json%add(a2,o3)
             call json%add(o3,'id',l)
             call json%add(o3,'attractor_id',j)
             call json%add(o3,'cell_cp',trim(scp))
             call json%add(o3,'lvec',sy%c%mol(idxmol(1,j))%at(idxmol(2,j))%lvec(k))
             call json%add(o3,'nonequivalent_cp',trim(sncp))
             call json%add(o3,'name',trim(sname))
             call json%add(o3,'atomic_number',trim(sz))
             call json%add(o3,'fractional_coordinates',x)
             do k = 1, nprop
                sump(k) = sump(k) + res(nacprop(k))%psum(j)
             end do
          end do
          call json%add(o2,'integrals',sump)
          nullify(a2)
          nullify(o2)
       end do
       nullify(a1)
       deallocate(nacprop,idxmol,sump)
    end if

    ! integration options
    call json%add(o1,'noatoms',.not.bas%atexist)
    if (len_trim(bas%expr) > 0) &
       call json%add(o1,'discard',trim(bas%expr))
    call json%add(o1,'ratom',bas%ratom)

    ! print to file and finish
    call json%print(p,file)
    if (json%failed()) &
       call ferror("int_output_json","error writing JSON file",faterr)
    call json%destroy(p)

  end subroutine int_output_json

  !> Assign strings for attractor i. icp is the CP identifier for the
  !> attractor. usesym controls whether the multiplicity will be
  !> written to the output. The output strings are scp (the cell CP
  !> identifier, or -- if not a CP attractor), sncp (the
  !> non-equivalent CP identifier, or -- if not a CP attractor), sname
  !> (the name of the CP), smult (CP multiplicity), and sz (atomic
  !> number, if the CP is an atom).
  subroutine assign_strings(i,icp,usesym,scp,sncp,sname,smult,sz)
    use systemmod, only: sy
    use tools_io, only: string, ioj_left, ioj_center
    integer, intent(in) :: i
    integer, intent(in) :: icp
    logical, intent(in) :: usesym
    character(len=:), allocatable, intent(out) :: sncp, scp, sname, sz, smult

    integer :: idx

    if (icp > 0) then
       ! this is a cp
       idx = sy%f(sy%iref)%cpcel(icp)%idx
       scp = string(icp,4,ioj_left)
       sncp = string(idx,4,ioj_left)
       sname = string(sy%f(sy%iref)%cp(idx)%name,6,ioj_center)
       smult = string(sy%f(sy%iref)%cp(idx)%mult,4,ioj_center)
       if (sy%f(sy%iref)%cp(idx)%isnuc) then
          sz = string(sy%c%spc(sy%c%at(idx)%is)%z,2,ioj_left)
       else
          sz = "--"
       endif
    else
       ! this is an unknown nnm
       scp = " -- "
       sncp = " -- "
       sname = "  ??  "
       smult = string(1,4,ioj_center)
       sz = "--"
    end if
    if (.not.usesym) smult = " -- "

  end subroutine assign_strings

  !> Plot the atomic basins found by YT or BADER to graphical files.
  !> bas = integration driver data.
  subroutine int_gridbasins(bas)
    use iso_c_binding, only: c_ptr
    use yt, only: yt_weights, ytdata_clean, ytdata
    use systemmod, only: sy
    use global, only: fileroot
    use graphics, only: grhandle
    use tools_math, only: m_x2c_from_cellpar, matinv, cross
    use tools_io, only: string, uout, ferror, noerr, string
    use types, only: realloc, basindat
    type(basindat), intent(in) :: bas

    character(len=:), allocatable :: str
    integer :: i, j, i1, i2, i3, iaux, n(3), q(3), p(3)
    real*8 :: d2, x(3), x1(3), x2(3), xd(3), xrhotmp
    real*8, allocatable :: w(:,:,:)
    integer, allocatable :: idg0(:,:,:)
    type(ytdata) :: dat
    type(grhandle) :: gr
    integer :: nvert, nf
    real*8, allocatable :: xvert(:,:), xrho(:)
    integer, allocatable :: iface(:,:)
    type(c_ptr), target :: fid

    interface
       ! The definitions and documentation for these functions are in doqhull.c
       subroutine runqhull_basintriangulate_step1(n,x0,xvert,nf,fid) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
         integer(c_int), value :: n
         real(c_double) :: x0(3)
         real(c_double) :: xvert(3,n)
         integer(c_int) :: nf
         type(c_ptr) :: fid
       end subroutine runqhull_basintriangulate_step1
       subroutine runqhull_basintriangulate_step2(nf,iface,fid) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
         integer(c_int), value :: nf
         integer(c_int) :: iface(3,nf)
         type(c_ptr), value :: fid
       end subroutine runqhull_basintriangulate_step2
    end interface

    if (bas%ndrawbasin < 0) return
    if (bas%imtype == imtype_hirshfeld) return

    ! prepare wigner-seitz tetrahedra
    do i = 1, 3
       n(i) = size(bas%idg,i)
    end do

    ! output
    write (uout,'("+ Basins written to ",A,"_basins-*.",A/)') trim(fileroot), bas%basinfmt

    ! write the unit cell
    str = trim(fileroot) // "_basins-cell." // bas%basinfmt
    call sy%c%write_3dmodel(str,bas%basinfmt,(/1,1,1/),.true.,.false.,sy%c%ismol3d,&
       .true.,.true.,-1d0,(/0d0,0d0,0d0/),-1d0,(/0d0,0d0,0d0/))

    ! prepare the idg array
    allocate(idg0(n(1),n(2),n(3)))
    idg0 = 0
    if (bas%imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
       call yt_weights(luw=bas%luw,dout=dat)
       do i = 1, bas%nattr
          call yt_weights(din=dat,idb=i,w=w)
          where(w >= 0.5d0)
             idg0 = i
          end where
       end do
       call ytdata_clean(dat)
       deallocate(w)
    else
       idg0 = bas%idg
    endif

    ! write the basins
    allocate(xvert(3,10),xrho(10))
    do i = 1, bas%nattr
       if (bas%ndrawbasin > 0 .and. bas%ndrawbasin /= i) cycle

       nvert = 0
       !$omp parallel do private(p,q,x,xd,d2,xrhotmp)
       do i1 = 1, n(1)
          do i2 = 1, n(2)
             do i3 = 1, n(3)
                if (idg0(i1,i2,i3) /= i) cycle

                ! is this point on the border of the basin?
                p = (/i1,i2,i3/)
                do j = 1, sy%f(sy%iref)%grid%nvec
                   q = modulo(p + sy%f(sy%iref)%grid%vec(:,j) - 1,n) + 1
                   if (idg0(q(1),q(2),q(3)) /= i) then
                      ! move the point to the WS of the attractor and convert to Cartesian
                      x = real(p-1,8) / n
                      xd = x - bas%xattr(:,i)
                      call sy%c%shortest(xd,d2)
                      x = sy%c%x2c(bas%xattr(:,i)) + xd
                      xrhotmp = bas%f(q(1),q(2),q(3))

                      ! add to the list of basin points
                      !$omp critical (addvertex)
                      nvert = nvert + 1
                      if (nvert > size(xvert,2)) then
                         call realloc(xvert,3,2*nvert)
                         call realloc(xrho,2*nvert)
                      end if
                      xvert(:,nvert) = x
                      xrho(nvert) = xrhotmp
                      !$omp end critical (addvertex)
                   end if
                end do
             end do
          end do
       end do
       !$omp end parallel do

       if (nvert > 0) then
          ! run qhull
          x = sy%c%x2c(bas%xattr(:,i))
          call runqhull_basintriangulate_step1(nvert,x,xvert,nf,fid)
          allocate(iface(3,nf))
          call runqhull_basintriangulate_step2(nf,iface,fid)

          ! orient the faces
          !$omp parallel do private(x1,x2,iaux)
          do j = 1, nf
             x1 = xvert(:,iface(2,j)) - xvert(:,iface(1,j))
             x2 = xvert(:,iface(3,j)) - xvert(:,iface(1,j))
             x1 = cross(x1,x2)

             x2 = (xvert(:,iface(1,j)) + xvert(:,iface(2,j)) + xvert(:,iface(3,j))) / 3d0 - x
             if (dot_product(x1,x2) < 0d0) then
                iaux = iface(1,j)
                iface(1,j) = iface(2,j)
                iface(2,j) = iaux
             end if
          end do
          !$omp end parallel do

          ! write the triangulation to a file
          str = trim(fileroot) // "_basins-" // string(i) // "." // bas%basinfmt
          call gr%open(bas%basinfmt,str)
          if (bas%imtype == imtype_isosurface) then
             call gr%triangulation(nvert,xvert,nf,iface)
          else
             call gr%triangulation(nvert,xvert,nf,iface,xrho)
          end if
          call gr%close()
          deallocate(iface)
       else
          call ferror("int_gridbasins","Basin " // string(i) // " has zero volume.",noerr)
       end if
    end do
    deallocate(idg0,xvert,xrho)

  end subroutine int_gridbasins

  !> Write cube files with the integration weights.  bas = integration
  !> driver data.
  subroutine int_cubew(bas)
    use hirshfeld, only: hirsh_weights
    use systemmod, only: sy
    use yt, only: ytdata, yt_weights, ytdata_clean
    use global, only: fileroot
    use tools_io, only: uout, string
    use types, only: basindat
    type(basindat), intent(in) :: bas

    type(ytdata) :: dat
    real*8, allocatable :: w(:,:,:)
    character(len=:), allocatable :: file
    integer :: i

    if (.not.bas%wcube) return
    allocate(w(bas%n(1),bas%n(2),bas%n(3)))
    if (bas%imtype == imtype_yt) then
       w = 0d0
       call yt_weights(luw=bas%luw,dout=dat)
    end if

    do i = 1, bas%nattr
       if (bas%imtype == imtype_yt) then
          call yt_weights(din=dat,idb=i,w=w)
       elseif (bas%imtype == imtype_hirshfeld) then
          call hirsh_weights(sy,bas,i,w)
       else
          w = 0d0
          where (bas%idg == i)
             w = 1d0
          end where
       endif
       file = trim(fileroot) // "_wcube_" // string(i,2,pad0=.true.) // ".cube"
       call sy%c%writegrid_cube(w,file,.false.,.false.)
    end do
    write (uout,'("+ Weights written to ",A,"_wcube_*.cube"/)') trim(fileroot)

    if (bas%imtype /= imtype_yt .and. bas%imtype /= imtype_hirshfeld) then
       w = bas%idg
       file = trim(fileroot) // "_wcube_all.cube"
       call sy%c%writegrid_cube(w,file,.false.,.false.)
       write (uout,'("+ All weights written to ",A,"_wcube_*.cube"/)') trim(fileroot)
    end if

    deallocate(w)

    if (bas%imtype == imtype_yt) then
       call ytdata_clean(dat)
    end if

  end subroutine int_cubew

  !> Unpacking routine for use in delocalization index calculations.
  subroutine unpackidx(idx,io,jo,ko,bo,nmo,nbnd,nlat)
    integer, intent(in) :: idx, nmo, nbnd, nlat(3)
    integer, intent(out) :: io, jo, ko, bo

    integer :: iaux

    ! unpack
    iaux = modulo(idx-1,nmo)
    bo = modulo(iaux,nbnd)
    iaux = (idx-1 - bo) / nbnd
    ko = modulo(iaux,nlat(3))
    iaux = (iaux - ko) / nlat(3)
    jo = modulo(iaux,nlat(2))
    iaux = (iaux - jo) / nlat(2)
    io = modulo(iaux,nlat(1))
    bo = bo + 1

  end subroutine unpackidx

  !> Packing routine for use in delocalization index calculations.
  subroutine packidx(io,jo,ko,bo,idx,nmo,nbnd,nlat)
    integer, intent(in) :: io, jo, ko, bo, nmo, nbnd, nlat(3)
    integer, intent(out) :: idx

    integer :: zio, zjo, zko

    ! transformed indices
    zio = modulo(io,nlat(1))
    zjo = modulo(jo,nlat(2))
    zko = modulo(ko,nlat(3))

    ! translate and pack
    idx = bo + nbnd * (zko + nlat(3) * (zjo + nlat(2) * zio))

  end subroutine packidx

end submodule proc
