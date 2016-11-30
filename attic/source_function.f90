! in: fields.f90

  integer, parameter, public :: itype_source = 10
  character*10, parameter, public :: itype_names(10) = (/&
     "Volume    ","Field     ","Field (v) ","Gradnt mod","Laplacian ",&
     "Laplcn (v)","Expression","Multipoles","Deloc indx","Source fun"/)

! [...]

    ! add property
    nprops = nprops + 1
    if (nprops > size(integ_prop)) &
       call realloc(integ_prop,2*nprops)
    if (useexpr) then
       integ_prop(nprops)%used = .true.
       integ_prop(nprops)%fid = 0
       integ_prop(nprops)%itype = itype_expr
       integ_prop(nprops)%prop_name = trim(adjustl(expr))
       integ_prop(nprops)%expr = expr
    else
       integ_prop(nprops)%used = .true.
       integ_prop(nprops)%fid = id
       integ_prop(nprops)%itype = itype_f
       integ_prop(nprops)%prop_name = ""
       integ_prop(nprops)%lmax = 5
       do while (.true.)
          word = lgetword(line,lp)
          if (equal(word,"f")) then
             integ_prop(nprops)%itype = itype_f
             str = trim(str) // "#f"
          elseif (equal(word,"fval")) then
             integ_prop(nprops)%itype = itype_fval
             str = trim(str) // "#fval"
          elseif (equal(word,"gmod")) then
             integ_prop(nprops)%itype = itype_gmod
             str = trim(str) // "#gmod"
          elseif (equal(word,"lap")) then
             integ_prop(nprops)%itype = itype_lap
             str = trim(str) // "#lap"
          elseif (equal(word,"lapval")) then
             integ_prop(nprops)%itype = itype_lapval
             str = trim(str) // "#lapval"
          elseif (equal(word,"multipoles") .or. equal(word,"multipole")) then
             integ_prop(nprops)%itype = itype_mpoles
             str = trim(str) // "#mpol"
             ok = isinteger(idum,line,lp)
             if (ok) integ_prop(nprops)%lmax = idum
          elseif (equal(word,"deloc")) then
             integ_prop(nprops)%itype = itype_deloc
             str = trim(str) // "#deloca"
          elseif (equal(word,"deloc")) then
             integ_prop(nprops)%itype = itype_deloc
             str = trim(str) // "#deloca"
          elseif (equal(word,"source")) then
             integ_prop(nprops)%itype = itype_source
             str = trim(str) // "#source"
             integ_prop(nprops)%prop_name = str
             ok = eval_next(integ_prop(nprops)%x0(1),line,lp)
             ok = ok.and.eval_next(integ_prop(nprops)%x0(2),line,lp)
             ok = ok.and.eval_next(integ_prop(nprops)%x0(3),line,lp)
             if (.not.ok) then
                ! read an environment
                do while(.true.)
                   ok = getline(uin,oline,.true.,ucopy)
                   if (.not.ok) &
                      call ferror("fields_integrable","eof reading source function env.",faterr)
                   lp = 1
                   word = getword(oline,lp)
                   if (equal(word,"end").or.equal(word,"endintegrable")) then
                      nprops = nprops - 1
                      exit
                   else
                      lp = 1
                      ok = eval_next(x0(1),oline,lp)
                      ok = ok.and.eval_next(x0(2),oline,lp)
                      ok = ok.and.eval_next(x0(3),oline,lp)
                      if (.not.ok) &
                         call ferror("fields_integrable","error reading source function env.",faterr)
                   end if
                   integ_prop(nprops) = integ_prop(nprops-1)
                   integ_prop(nprops)%x0 = x0
                   if (cr%ismolecule) then
                      integ_prop(nprops)%x0 = integ_prop(nprops)%x0 / dunit - cr%molx0
                      integ_prop(nprops)%x0 = cr%c2x(integ_prop(nprops)%x0)
                   end if
                   nprops = nprops + 1
                   if (nprops > size(integ_prop)) &
                      call realloc(integ_prop,2*nprops)
                end do
                exit
             else
                ! convert the point read in a single line
                if (cr%ismolecule) then
                   integ_prop(nprops)%x0 = integ_prop(nprops)%x0 / dunit - cr%molx0
                   integ_prop(nprops)%x0 = cr%c2x(integ_prop(nprops)%x0)
                end if
             end if
          elseif (equal(word,"name")) then
             word = getword(line,lp)
             integ_prop(nprops)%prop_name = string(word)
          else if (len_trim(word) > 0) then
             call ferror("fields_integrable","Unknown extra keyword",faterr,line,syntax=.true.)
             nprops = nprops - 1
             return
          else
             exit
          end if
       end do
       if (trim(integ_prop(nprops)%prop_name) == "") then
          integ_prop(nprops)%prop_name = str
       end if
    end if

! [...]

       case(itype_source)
          sprop = "sf"

! [...]

       elseif (integ_prop(i)%itype == itype_source) then
          if (cr%ismolecule) then
             x0 = cr%x2c(integ_prop(i)%x0)
             x0 = (x0 + cr%molx0) * dunit
             sunit = "("// iunitname0(iunit) //")"
          else
             x0 = integ_prop(i)%x0
             sunit = "(cryst.)"
          end if
          write (uout,'(2X,4(A,2X),99(A))') &
             string(i,length=3,justify=ioj_right), string(sprop,length=4,justify=ioj_center), &
             string(integ_prop(i)%fid,length=5,justify=ioj_right), &
             string(integ_prop(i)%prop_name,10,ioj_left),&
             "at point: (", (trim(string(x0(j),'f',10,6))//", ",j=1,2),&
             trim(string(x0(3),'f',10,6)),") ", sunit

! in integration.f90:

       elseif (integ_prop(k)%itype == itype_source) then
          reason(k) = "Source function integrated separately (see table below)"
          cycle

! [...]

    ! source function
    call intgrid_source(fint,idprop,nattr,xgatt,idg,imtype,luw,sf)

! [...]

  !> Calculate the source function at given points using integration
  !> on a grid. The input is: the array containing the integrable
  !> field information on the basin field grid (fint), the index array
  !> relating 1->nprops to the fint array, the calculated number of
  !> attractors (natt), their positions in crystallographic
  !> coordinates (xgatt), the attractor assignment for each of the
  !> grid nodes (idg), the integration type (imtype: bader or yt), the
  !> weight file for YT, and the output multipolar moments.
  subroutine intgrid_source(fint,idprop,natt,xgatt,idg,imtype,luw,sf)
    use grid_tools, only: grid_laplacian
    use yt, only: yt_weights
    use fields, only: nprops, integ_prop, itype_source, f
    use struct_basic, only: cr
    use global, only: refden
    use tools_math, only: tosphere, genrlm_real
    use types, only: field
    use param, only: fourpi
    real*8, intent(in) :: fint(:,:,:,:)
    integer, intent(in) :: idprop(:)
    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: imtype
    integer, intent(in) :: luw
    real*8, allocatable, intent(inout) :: sf(:,:)

    integer :: i, j, k, l, n(3), ntot, np
    integer :: fid
    real*8, allocatable :: w(:,:,:)
    real*8 :: dv(3), p(3), r2
    type(field) :: faux

    ! allocate space for the source function
    np = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_source) cycle
       np = np + 1
    end do
    if (np == 0) return

    if (allocated(sf)) deallocate(sf)
    allocate(sf(natt,np))
    sf = 0d0

    ! size of the grid and YT weights
    n = f(refden)%n
    ntot = n(1)*n(2)*n(3)
    if (imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
    endif

    ! run over all properties for which the s.f. calculation is active
    np = 0
    do l = 1, nprops
       if (.not.integ_prop(l)%used) cycle
       if (.not.integ_prop(l)%itype == itype_source) cycle
       fid = integ_prop(l)%fid
       np = np + 1

       call grid_laplacian(f(fid),faux)

       do i = 1, n(1)
          do j = 1, n(2)
             do k = 1, n(3)
                p = real((/i-1,j-1,k-1/),8)
                dv = p/real(n,8) - integ_prop(l)%x0
                call cr%shortest(dv, r2)
                if (abs(r2) < 1d-8) cycle
                faux%f(i,j,k) = -faux%f(i,j,k) / (sqrt(r2)*fourpi)
             end do
          end do
       end do

       ! calcualate the source function contributions from each atom
       if (imtype == imtype_yt) then
          rewind(luw)
       endif
       do i = 1, natt
          if (imtype == imtype_yt) then
             call yt_weights(luw,i,w)
          end if
          if (imtype == imtype_bader) then
             w = 0d0
             where (idg == i)
                w = faux%f
             end where
             sf(i,np) = sum(w) * cr%omega / ntot
          else
             sf(i,np) = sum(w * faux%f) * cr%omega / ntot
          endif
       end do
    end do
    if (imtype == imtype_yt) then
       deallocate(w)
    endif

  end subroutine intgrid_source

! [...]

    ! Source function output
    if (present(sf)) then
       if (allocated(sf)) then
          write (uout,'("* Source function calculation ")')
          write (uout,'("  Please cite: ")')
          write (uout,'("  Carlo Gatti, Richard W. Bader, Chem. Phys. Lett. 287 (1998) 233.")')
          write (uout,'("  Christian Tantardini, Davide Ceresoli, Enrico Benassi, J. Comput. Chem. 37 (2016) 2133–2139.")')
          n = 0
          do l = 1, nprops
             if (.not.integ_prop(l)%used) cycle
             if (.not.integ_prop(l)%itype == itype_source) cycle
             n = n + 1
             write (uout,*)
             write (*,*) " property : ", l
             write (*,*) " point : ", integ_prop(l)%x0
             do i = 1, nattr
                write (*,*) " attractor : ", i, sf(i,n)
             end do
          end do
       end if
    end if


