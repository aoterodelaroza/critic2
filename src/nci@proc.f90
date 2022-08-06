! Copyright (c) 2007-2018 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

submodule (nci) proc
  implicit none

  !xx! private procedures
  ! subroutine write_cube_header(lu,l1,l2,periodic,nfrag,frag,x0x,x1x,x0,x1,nstep,xmat)
  ! subroutine write_cube_body(lu,n,c)
  ! function check_no_extra_word(line,lp,routine)
  ! function read_fragment(lu) result(fr)
  ! function readchk(x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag) result(lchk)
  ! subroutine writechk(x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag)

contains

  ! nci plots
  module subroutine nciplot()
    use systemmod, only: sy
    use fieldmod, only: field, type_grid, type_promol
    use struct_drivers, only: struct_write
    use grid1mod, only: agrid
    use global, only: fileroot, eval_next, dunit0, quiet, iunit, iunitname0
    use fragmentmod, only: fragment, realloc_fragment
    use tools_io, only: getline, lgetword, equal, uin, faterr, ferror, ucopy, &
       string, getword, uout, fopen_write, tictac, fclose
    use tools_math, only: eigsym, m_x2c_from_cellpar, matinv
    use types, only: scalar_value, realloc
    use param, only: pi, vsmall, bohrtoa, ifformat_as_grad, ifformat_as_hxx1,&
       ifformat_as_hxx2, ifformat_as_hxx3
    type(field) :: fgrho, fxx(3)
    type(scalar_value) :: res, resg
    character(len=:), allocatable :: line, word, oname, file
    logical :: ok, ok2, nstep_from_grid
    integer :: lp, istat
    integer :: i, j, k, iat, ifr, l, nmol
    integer :: lugc, ludc, luvmd, ludat, lupc
    real*8, allocatable, dimension(:,:,:) :: crho, cgrad, rhoat
    real*8, allocatable, dimension(:,:,:,:) :: rhofrag
    real*8 :: ehess(3), dimgrad, dist, rrho, rrho1, rrho2
    real*8 :: sumq, sumqp, sumv, rdum1, rdum2
    logical :: onlyneg, lchk, inter, isden
    integer :: findlimits, ithres, istep, iz
    logical :: periodic, usecore
    ! molmotif
    logical :: domolmotif
    logical, allocatable :: isdiscrete(:)
    ! chk
    logical :: usechk, dopromol
    ! fragments
    integer :: nfrag
    type(fragment), allocatable :: fr(:)
    type(fragment) :: fr0
    ! voids
    real*8 :: rho_void
    ! cube limits
    real*8 :: x(3), xd(3), x0(3), x1(3), x0x(3), x1x(3), xinc(3), xmat(3,3), xmati(3,3)
    real*8 :: xx0(3), xx1(3)
    integer :: nstep(3), imin(3), imax(3)
    ! cutoffs
    real*8 :: rhocut, dimcut, rhoplot, dimplot, rhoatl
    real*8, allocatable :: rhofragl(:)
    ! discarding rho parameter
    real*8 :: rhoparam, rhoparam2
    ! rthres
    real*8 :: rthres, srhorange(2)
    ! Cartesian matrix for the vmd coordinate system
    real*8 :: rchol(3,3), gg(3,3), aal(3), bbl(3), delta(3)
    logical :: isortho

    ! named constants
    real*8, parameter :: fthirds = 4d0/3d0
    real*8, parameter :: const = 2.d0*(3.d0*pi**2)**(1.d0/3.d0)
    real*8, parameter :: rthres_xyz = 1d0

    real*8, parameter :: xlist0(3,12) = reshape((/&
       0d0, 0d0, 0d0,&
       0d0, 0d0, 0d0,&
       0d0, 0d0, 0d0,&
       0d0, 0d0, 1d0,&
       0d0, 1d0, 1d0,&
       0d0, 1d0, 0d0,&
       1d0, 1d0, 0d0,&
       1d0, 0d0, 0d0,&
       1d0, 0d0, 1d0,&
       1d0, 1d0, 1d0,&
       1d0, 1d0, 1d0,&
       1d0, 1d0, 1d0/),shape(xlist0))
    real*8, parameter :: xlist1(3,12) = reshape((/&
       1d0, 0d0, 0d0,&
       0d0, 1d0, 0d0,&
       0d0, 0d0, 1d0,&
       0d0, 1d0, 1d0,&
       0d0, 1d0, 0d0,&
       1d0, 1d0, 0d0,&
       1d0, 0d0, 0d0,&
       1d0, 0d0, 1d0,&
       0d0, 0d0, 1d0,&
       0d0, 1d0, 1d0,&
       1d0, 0d0, 1d0,&
       1d0, 1d0, 0d0/),shape(xlist1))


    ! default values
    isden = .not.(sy%iref == 0)
    rhoparam = 0.95d0
    rhoparam2 = 0.75d0
    if (sy%c%ismolecule) then
       x0 = sy%c%molborder
       x1 = 1d0 - sy%c%molborder
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
       periodic = .false.
    else
       x0 = 0d0
       x1 = 1d0
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
       periodic = .true.
    end if
    xinc = 0.1d0
    oname = fileroot
    rhocut = 0.2d0
    nstep = -1
    findlimits = -1
    ithres = -1
    istep = -1

    ! initialize fragments
    nfrag = 0
    allocate(fr(1))

    onlyneg = .false.
    usechk = .true.
    rho_void = -1d0
    if (isden) then
       dimcut = 2.0d0
       dimplot=0.5d0
       rhoplot=0.1d0
    else
       dimcut = 1.0d0
       dimplot=0.3d0
       rhoplot=0.12d0
    end if
    rthres = 2.d0
    domolmotif = .false.
    usecore = sy%f(sy%iref)%usecore .and. any(sy%f(sy%iref)%zpsp /= -1)
    srhorange = (/ -1d30, 1d30 /)

    do while(.true.)
       ok = getline(uin,line,.true.,ucopy)
       lp = 1
       word = lgetword(line,lp)

       if (equal(word,'oname')) then
          oname = line(lp:)
       elseif (equal(word,'molmotif')) then
          domolmotif = .true.
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'cutoffs')) then
          ok = eval_next(rhocut,line,lp)
          ok = ok .and. eval_next(dimcut,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong cutoffs keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'rhoparam')) then
          ok = eval_next(rhoparam,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong rhoparam keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'rhoparam2')) then
          ok = eval_next(rhoparam2,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong rhoparam2 keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'cutplot')) then
          ok = eval_next(rhoplot,line,lp)
          ok = ok .and. eval_next(dimplot,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong cutplot keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'srhorange')) then
          ok = eval_next(rdum1,line,lp)
          ok = ok .and. eval_next(rdum2,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong cutplot keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return

          srhorange(1) = min(rdum1,rdum2)
          srhorange(2) = max(rdum1,rdum2)
       elseif (equal(word,'void')) then
          ok = eval_next(rho_void,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong void keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'rthres')) then
          ok = eval_next(rthres,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong cutplot keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          rthres = rthres / dunit0(iunit)
       elseif (equal(word,'increments')) then
          istep = 0
          ok = eval_next(xinc(1),line,lp)
          ok = ok .and. eval_next(xinc(2),line,lp)
          ok = ok .and. eval_next(xinc(3),line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong increments keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          xinc = xinc / dunit0(iunit)
       elseif (equal(word,'nstep')) then
          istep = 1
          ok = eval_next(nstep(1),line,lp)
          ok = ok .and. eval_next(nstep(2),line,lp)
          ok = ok .and. eval_next(nstep(3),line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong nstep keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          do i = 1, 3
             nstep(i) = max(nstep(i),2)
          end do
       elseif (equal(word,'onlyneg')) then
          onlyneg = .true.
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'nochk')) then
          usechk = .false.
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'cube')) then
          findlimits = 0
          periodic = .false.
          ok = eval_next(x0(1),line,lp)
          if (ok) then
             ok = ok .and. eval_next(x0(2),line,lp)
             ok = ok .and. eval_next(x0(3),line,lp)
             ok = ok .and. eval_next(x1(1),line,lp)
             ok = ok .and. eval_next(x1(2),line,lp)
             ok = ok .and. eval_next(x1(3),line,lp)
             if (.not.ok) then
                call ferror('nciplot','wrong cube syntax',faterr,line,syntax=.true.)
                return
             end if
             ok = check_no_extra_word(line,lp,'nciplot')
             if (.not.ok) return
             if (.not.sy%c%ismolecule) then
                x0 = sy%c%x2c(x0)
                x1 = sy%c%x2c(x1)
             else
                x0 = x0 / dunit0(iunit) - sy%c%molx0
                x1 = x1 / dunit0(iunit) - sy%c%molx0
             endif
             ithres = 0
          else
             x0 = 1d30
             x1 = -1d30
             ok = .false.
             do while (.true.)
                word = getword(line,lp)
                inquire(file=word,exist=ok2)
                if (.not.ok2) exit
                ok = .true.
                fr0 = sy%c%identify_fragment_from_xyz(word)
                do j = 1, fr0%nat
                   do k = 1, 3
                      x0(k) = min(x0(k),fr0%at(j)%r(k))
                      x1(k) = max(x1(k),fr0%at(j)%r(k))
                   end do
                end do
             end do
             ithres = 1
          endif
          if (.not.ok) then
             call ferror('nciplot','wrong cube keyword',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'fragment')) then
          if (findlimits == -1) findlimits = 1
          if (ithres == -1) ithres = 1
          nfrag = nfrag + 1
          if (nfrag > size(fr)) call realloc_fragment(fr,2*nfrag)
          word = getword(line,lp)
          if (equal(word,'')) then
             fr(nfrag) = read_fragment(uin)
          else
             fr(nfrag) = sy%c%identify_fragment_from_xyz(word)
          endif

       elseif (equal(word,'endnciplot').or.equal(word,'end')) then
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          exit
       else
          call ferror ('nciplot','unknown keyword',faterr,line,syntax=.true.)
          return
       end if
    end do
    if (nfrag > 0) call realloc_fragment(fr,nfrag)

    ! header
    write (uout,'("* NCIPLOT: non-covalent interactions")')
    write (uout,'("  Please cite:")')
    write (uout,'("  E. R. Johnson et al., J. Am. Chem. Soc., 132 (2010) 6498.")')
    write (uout,'("  A. Otero-de-la-Roza et al., Phys. Chem. Chem. Phys., 14 (2012) 12165.")')
    if (.not.quiet) call tictac("start nciplot")

    ! redefine the cube using the fragments
    if (findlimits == 1) then
       if (nfrag < 1) call ferror('nciplot','no fragments but findlimits=1',faterr)
       x0 = 1d30
       x1 = -1d30
       do i = 1, nfrag
          do j = 1, fr(i)%nat
             do k = 1, 3
                x0(k) = min(x0(k),fr(i)%at(j)%r(k))
                x1(k) = max(x1(k),fr(i)%at(j)%r(k))
             end do
          end do
       end do
       periodic = .false.
    endif

    ! apply the threshold
    if (ithres == 1) then
       x0 = x0 - rthres
       x1 = x1 + rthres
    end if

    ! convert to crystallographic coordinates
    x0 = sy%c%c2x(x0)
    x1 = sy%c%c2x(x1)
    x0x = x0
    x1x = x1

    ! determine the grid step
    nstep_from_grid = .false.
    if (periodic) then
       ! use xinc for periodic cells
       do i = 1, 3
          x1 = 0d0
          x1(i) = 1d0
          x1 = sy%c%x2c(x1)
          if (istep == -1 .and. sy%f(sy%iref)%type == type_grid) then
             nstep_from_grid = .true.
             nstep = sy%f(sy%iref)%grid%n
          elseif (istep == 0 .or. istep == -1) then
             nstep(i) = ceiling(norm2(x1) / xinc(i))
          end if
          xmat(:,i) = x1 / real(nstep(i),8)
          xinc(i) = norm2(xmat(:,i))
       end do
     else
       ! use the grid step
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
       xx0 = min(x0,x1)
       xx1 = max(x0,x1)
       x0 = xx0
       x1 = xx1

       if (istep == 0 .or. istep == -1) then
          nstep = ceiling(abs(x1 - x0) / xinc)
       else
          xinc = (x1-x0) / real(nstep,8)
       end if
       xmat = 0d0
       do i = 1, 3
          xmat(i,i) = xinc(i)
       end do
    end if

    ! use the promolecular density?
    dopromol = (nfrag > 0 .or. rho_void > 0d0)

    ! do stuff
    write(uout,'("+ File root: ",A)') trim(oname)
    write(uout,'("  Density cube file: ",A)') trim(oname) // "-dens.cube"
    write(uout,'("  RDG cube file: ",A)') trim(oname) // "-grad.cube"
    if (dopromol) then
       write(uout,'("  Promolecular density cube file: ",A)') trim(oname) // "-pdens.cube"
    end if
    write(uout,'("  RDG vs. density file: ",A)') trim(oname) // ".dat"
    write(uout,'("  VMD visualization script: ",A)') trim(oname) // ".vmd"
    write(uout,'("+ Density cutoff (for dat file): ",A)') string(rhocut,'f',decimal=6)
    write(uout,'("  RDG cutoff (for dat file): ",A)') string(dimcut,'f',decimal=6)
    write(uout,'("  Density plot cutoff: ",A)') string(rhoplot,'f',decimal=6)
    write(uout,'("  RDG plot cutoff: ",A)') string(dimplot,'f',decimal=6)
    write(uout,'("+ Cube origin (crystal coord.): ",3(A,X))') (string(x0x(j),'f',decimal=6),j=1,3)
    write(uout,'("  Cube end (crystal coord.): ",3(A,X))') (string(x1x(j),'f',decimal=6),j=1,3)
    write(uout,'("  Cube nodes in each direction: ",3(A,X))') (string(nstep(j)),j=1,3)
    write(uout,'("  Cube steps in each direction: ",3(A,X))') (string(xinc(j),'f',decimal=6),j=1,3)
    write(uout,'("  Cube side lengths in each direction (",A,"): ",3(A,X))') &
       iunitname0(iunit), (string(xinc(j)*nstep(j)*dunit0(iunit),'f',decimal=4),j=1,3)
    if (nstep_from_grid) &
       write(uout,'("  Values of NSTEP derived from reference field grid dimensions.")')
    if (sy%f(sy%iref)%type == type_grid.and..not.usecore) then
       write(uout,'("  Using fast Fourier transform for derivative calculation.")')
    else
       write(uout,'("  Using scalar field differentiation for derivative calculation.")')
    end if
    write(uout,*)

    ! allocate logical units and open files
    lugc = fopen_write(trim(oname)//"-grad.cube")
    ludc = fopen_write(trim(oname)//"-dens.cube")
    luvmd = fopen_write(trim(oname)//".vmd")
    ludat = fopen_write(trim(oname)//".dat")
    if (lugc < 0) call ferror('nciplot','error opening gradient cube',faterr)
    if (ludc < 0) call ferror('nciplot','error opening density cube',faterr)
    if (luvmd < 0) call ferror('nciplot','error opening vmd file',faterr)
    if (ludat < 0) call ferror('nciplot','error opening dat file',faterr)
    if (dopromol) then
       lupc = fopen_write(trim(oname)//"-pdens.cube")
       if (lupc < 0) call ferror('nciplot','error opening promolecular density cube',faterr)
    end if

    ! write cube headers
    call write_cube_header(lugc,'grad_cube','3d plot, reduced density gradient',periodic,nfrag,fr,x0x,x1x,x0,x1,nstep,xmat,isortho)
    call write_cube_header(ludc,'dens_cube','3d plot, density',periodic,nfrag,fr,x0x,x1x,x0,x1,nstep,xmat,isortho)
    if (dopromol) then
       call write_cube_header(lupc,'pdens_cube','3d plot, promolecular density',periodic,nfrag,fr,x0x,x1x,x0,x1,nstep,xmat,isortho)
    end if

    ! allocate memory for density and gradient
     allocate(crho(0:nstep(3)-1,0:nstep(2)-1,0:nstep(1)-1),stat=istat)
     if (istat /= 0) call ferror('nciplot','could not allocate memory for density cube',faterr)
     allocate(cgrad(0:nstep(3)-1,0:nstep(2)-1,0:nstep(1)-1),stat=istat)
     if (istat /= 0) call ferror('nciplot','could not allocate memory for grad',faterr)
     if (dopromol) then
        allocate(rhoat(0:nstep(3)-1,0:nstep(2)-1,0:nstep(1)-1),stat=istat)
        if (istat /= 0) call ferror('nciplot','could not allocate memory for rhoat',faterr)
     endif
     if (nfrag > 0) then
        allocate(rhofrag(0:nstep(3)-1,0:nstep(2)-1,0:nstep(1)-1,nfrag),stat=istat)
        if (istat /= 0) call ferror('nciplot','could not allocate memory for rhofrag',faterr)
     endif

     if (sy%f(0)%type /= type_promol) then
        call ferror('nciplot','promolecular density not found',faterr)
     end if

     ! calculate density, rdg,... read from chkpoint if available
     lchk = .false.
     if (usechk) then
        lchk = readchk(oname,x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag)
     end if

     if (.not.lchk) then
        ! Initialize gradrho out of the omp loop
        if (sy%f(sy%iref)%type == type_grid .and..not.usecore) then
           call fgrho%load_as_fftgrid(sy%c,-1,"",sy%f(sy%iref)%grid,ifformat_as_grad)
           call fxx(1)%load_as_fftgrid(sy%c,-1,"",sy%f(sy%iref)%grid,ifformat_as_hxx1)
           call fxx(2)%load_as_fftgrid(sy%c,-1,"",sy%f(sy%iref)%grid,ifformat_as_hxx2)
           call fxx(3)%load_as_fftgrid(sy%c,-1,"",sy%f(sy%iref)%grid,ifformat_as_hxx3)
           ! Set trilinear interpolation to prevent wonky negative values in the derivatives
           ! of the field in the low-density regions (which may have been noisy to begin with).
           call fgrho%grid%setmode("trilinear")
           call fxx(1)%grid%setmode("trilinear")
           call fxx(2)%grid%setmode("trilinear")
           call fxx(3)%grid%setmode("trilinear")
        endif

        ! calculate density, rdg, rhoat
        allocate(rhofragl(nfrag))
        rhofragl = 0d0
        !$omp parallel do private (x,res,resg,ehess,dimgrad,rhoatl,&
        !$omp xd,dist,rrho,rrho1,rrho2,iz) firstprivate(rhofragl)
        do i = 0, nstep(1)-1
           do j = 0, nstep(2)-1
              do k = 0, nstep(3)-1
                 x = x0 + i*xmat(:,1) + j*xmat(:,2) + k*xmat(:,3)

                 ! calculate properties at x: rho and rdg
                 if (sy%f(sy%iref)%type == type_grid.and..not.usecore) then
                    call sy%f(sy%iref)%grd(x,0,res)
                    call fgrho%grd(x,0,resg)
                    dimgrad = resg%f / (const*max(res%f,vsmall)**fthirds)
                    do l = 1, 3
                       call fxx(l)%grd(x,0,resg)
                       ehess(l) = resg%f
                    end do
                    if (count(ehess > 0d0) >= 2) then
                       ehess(2) = 1d0
                    else
                       ehess(2) = -1d0
                    endif
                 else
                    ! strangely enough, eig takes about the same time as counting the signs
                    ! and using sylvester's law of inertia.
                    call sy%f(sy%iref)%grd(x,2,res)
                    call eigsym(res%hf,3,ehess)
                    dimgrad = res%gfmod / (const*max(res%f,vsmall)**fthirds)
                 end if

                 ! promolecular density
                 if (dopromol) then
                    rhoatl = sy%f(0)%grd0(x)
                 endif

                 ! fragments
                 if (nfrag > 0) then
                    rhofragl = 0d0
                    do ifr = 1, nfrag
                       do iat = 1, fr(ifr)%nat
                          xd = x - fr(ifr)%at(iat)%r
                          dist = norm2(xd)
                          iz = fr(ifr)%spc(fr(ifr)%at(iat)%is)%z
                          if (.not.agrid(iz)%isinit) cycle
                          if (dist > agrid(iz)%rmax) cycle
                          dist = max(dist,agrid(iz)%r(1))
                          dist = max(dist,1d-14)
                          call agrid(iz)%interp(dist,rrho,rrho1,rrho2)
                          rrho = max(rrho,0d0)
                          rhofragl(ifr) = rhofragl(ifr) + rrho
                       end do
                    end do
                 end if

                 !$omp critical (cubewrite)
                 cgrad(k,j,i) = dimgrad
                 crho(k,j,i) = sign(res%f,ehess(2))*100.D0
                 if (dopromol) rhoat(k,j,i) = rhoatl
                 if (nfrag > 0) rhofrag(k,j,i,1:nfrag) = rhofragl(1:nfrag)
                 !$omp end critical (cubewrite)
              end do
           end do
        end do
        !$omp end parallel do
        deallocate(rhofragl)

        ! save the ncichk file
        if (usechk) call writechk(oname,x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag)
     endif

     ! calculate cutoff effect on grid
     sumq = 0d0
     sumqp = 0d0
     sumv = 0d0
     do i = 0, nstep(1)-1
        do j = 0, nstep(2)-1
           do k = 0, nstep(3)-1
              x = x0 + i*xmat(:,1) + j*xmat(:,2) + k*xmat(:,3)

              ! ! calculate properties at x
              inter = .true.
              if (nfrag > 0) then
                 inter = (sum(rhofrag(k,j,i,1:nfrag)) >= rhoparam2 * rhoat(k,j,i)) .and. &
                    all(rhofrag(k,j,i,1:nfrag) <= sum(rhofrag(k,j,i,1:nfrag))*rhoparam)
              end if
              if (rho_void > 0d0) then
                 inter = inter.and.(rhoat(k,j,i)<rho_void)
                 if (rhoat(k,j,i)<rho_void) then
                    sumq = sumq + abs(crho(k,j,i)) / 100d0
                    sumqp = sumqp + rhoat(k,j,i)
                    sumv = sumv + 1
                 endif
              endif

              ! apply cutoffs and write the dat file
              if ((abs(crho(k,j,i)) < rhocut*100d0) .and. (cgrad(k,j,i) < dimcut) .and. inter) then
                 write(ludat,'(1p,E15.7,E15.7)') crho(k,j,i)/100d0, cgrad(k,j,i)
              end if
              if (abs(crho(k,j,i)) > rhoplot*100d0 .or..not.inter) cgrad(k,j,i) = 100d0
              if (crho(k,j,i) < srhorange(1)*100d0 .or. crho(k,j,i) > srhorange(2)*100d0 .or..not.inter) cgrad(k,j,i) = 100d0
              if (crho(k,j,i) > 0 .and. onlyneg) cgrad(k,j,i) = 100d0
           end do
        end do
     end do
     sumq = sumq * sy%c%omega / (nstep(1)*nstep(2)*nstep(3))
     sumqp = sumqp * sy%c%omega / (nstep(1)*nstep(2)*nstep(3))
     sumv = sumv * sy%c%omega / (nstep(1)*nstep(2)*nstep(3))

     if (rho_void > 0d0) then
        write (uout,'("* Void charge (a.u.): ",A)') string(sumq,'f',decimal=6)
        write (uout,'("* Void promolecular charge (a.u.): ",A)') string(sumqp,'f',decimal=6)
        write (uout,'("* Void volume (bohr^3): ",A)') string(sumv,'f',decimal=6)
     endif

     ! write cubes
     call write_cube_body(ludc,nstep,crho)
     call write_cube_body(lugc,nstep,cgrad)
     if (dopromol) then
        call write_cube_body(lupc,nstep,rhoat)
     end if

     ! deallocate grid files
     if (allocated(crho)) deallocate(crho)
     if (allocated(cgrad)) deallocate(cgrad)
     if (allocated(rhoat)) deallocate(rhoat)
     if (allocated(rhofrag)) deallocate(rhofrag)

     ! compute the atomic positions for the xyz file
     file = trim(oname) // "_cell.xyz"
     if (nfrag == 0) then
        if (periodic) then
           fr0 = sy%c%listatoms_cells((/1,1,1/),.true.)
           if (domolmotif) then
              call sy%c%listmolecules(fr0,nmol,fr,isdiscrete)
              call fr0%merge_array(fr,.false.)
           end if
           if (allocated(isdiscrete)) deallocate(isdiscrete)
        else
           call fr0%init()
           fr0%nspc = sy%c%nspc
           call realloc(fr0%spc,fr0%nspc)
           fr0%spc = sy%c%spc

           xx0 = sy%c%c2x(x0 - rthres_xyz)
           xx1 = sy%c%c2x(x1 + rthres_xyz)
           imin = floor(min(xx0,xx1)) - 1
           imax = ceiling(max(xx0,xx1)) + 1
           do i = imin(1), imax(1)
              do j = imin(2), imax(2)
                 do k = imin(3), imax(3)
                    do l = 1, sy%c%ncel
                       x = sy%c%x2c(sy%c%atcel(l)%x + real((/i,j,k/),8))
                       if (x(1) > x0(1) .and. x(1) < x1(1) .and.&
                           x(2) > x0(2) .and. x(2) < x1(2) .and.&
                           x(3) > x0(3) .and. x(3) < x1(3)) then
                          fr0%nat = fr0%nat + 1
                          if (fr0%nat > size(fr0%at)) call realloc(fr0%at,2*fr0%nat)
                          fr0%at(fr0%nat)%r = x
                          fr0%at(fr0%nat)%x = sy%c%c2x(x)
                          fr0%at(fr0%nat)%cidx = l
                          fr0%at(fr0%nat)%idx = sy%c%atcel(l)%idx
                          fr0%at(fr0%nat)%lvec = (/i,j,k/)
                          fr0%at(fr0%nat)%is = sy%c%atcel(l)%is
                       end if
                    end do
                 end do
              end do
           end do
        end if
        call realloc(fr0%at,fr0%nat)
     else
        call fr0%merge_array(fr(1:nfrag),.false.)
     end if
     if (fr0%nat == 0) &
        call ferror('nciplot','no atoms in region',faterr)

     ! VMD screw-up log, entry #1:
     ! Transform the coordinates in fr0 to vmd coordinates. vmd
     ! internally converts the coordinate system to a Cartesian matrix
     ! that is lower triangular.
     gg = matmul(transpose(xmat),xmat)
     aal(1) = sqrt(gg(1,1))
     aal(2) = sqrt(gg(2,2))
     aal(3) = sqrt(gg(3,3))
     if (any(abs(aal) < 1d-2)) &
        call ferror('nciplot','region too small',faterr)
     bbl(1) = acos(gg(2,3) / aal(2) / aal(3)) * 180d0 / pi
     bbl(2) = acos(gg(1,3) / aal(1) / aal(3)) * 180d0 / pi
     bbl(3) = acos(gg(1,2) / aal(1) / aal(2)) * 180d0 / pi
     isortho = all(abs(bbl - 90d0) < 1d-4)
     rchol = m_x2c_from_cellpar(aal,bbl)
     xmati = xmat
     call matinv(xmati,3)
     rchol = transpose(matmul(rchol,xmati))

     ! VMD screw-up log, entry #3:
     ! The origin of the isosurface box is translated by -1/2 of
     ! a step in each direction. Translate the atoms and the cell
     ! by this delta
     delta = -0.5d0 * (xmat(:,1) + xmat(:,2) + xmat(:,3))

     ! Transform the atoms and the cell to agree with vmd
     do i = 1, fr0%nat
        fr0%at(i)%r = matmul(fr0%at(i)%r + delta,rchol)
     end do

     ! write the fragment to the file
     call fr0%writexyz(file,.false.)

     ! write vmd script
     write (luvmd,'("#!/usr/local/bin/vmd")')
     write (luvmd,'("")')
     write (luvmd,'("# Display settings            ")')
     write (luvmd,'("display projection   Orthographic            ")')
     write (luvmd,'("display nearclip set 0.000000            ")')
     write (luvmd,'("display depthcue   off")')
     write (luvmd,'("light 2 on")')
     write (luvmd,'("light 3 on")')
     write (luvmd,'("axes location off")')
     write (luvmd,'("stage location off")')
     write (luvmd,'("color Display {Background} white")')
     write (luvmd,'("")')
     write (luvmd,'("# load the xyz with the atoms")')
     write (luvmd,'("mol new ",A," type xyz first 0 last -1 step 1 filebonds 0 autobonds 1 waitfor all")') trim(oname)//"_cell.xyz"
     write (luvmd,'("")')
     write (luvmd,'("# Render cell atoms as CPK")')
     write (luvmd,'("mol top 0")')
     write (luvmd,'("mol modstyle 0 0 CPK 1.000000 0.300000 50.000000 50.000000")')
     write (luvmd,'("mol modcolor 0 0 Element")')
     write (luvmd,'("mol modmaterial 0 0 Opaque")')
     write (luvmd,'("")')
     if (nfrag == 0 .and. periodic) then
        write (luvmd,'("# Cell")')
        write (luvmd,'("draw color red")')
        do i = 1, 12
           xx0 = matmul(sy%c%x2c(xlist0(:,i)) + delta,rchol)
           xx1 = matmul(sy%c%x2c(xlist1(:,i)) + delta,rchol)
           write (luvmd,'("draw cylinder {",3(F14.4,X),"} {",3(F14.4,X),"} radius 0.05 resolution 20")') &
              xx0 * bohrtoa, xx1 * bohrtoa
        end do
        write (luvmd,'("")')
     end if
     write (luvmd,'("# Load the cubes")')
     write (luvmd,'("mol new ",A," type cube first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all")')&
        trim(oname)//"-dens.cube"
     write (luvmd,'("mol addfile ",A," type cube first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all")')&
        trim(oname)//"-grad.cube"
     write (luvmd,'("")')
     write (luvmd,'("# representation of the atoms")')
     write (luvmd,'("mol top 1")')
     write (luvmd,'("mol modstyle 0 top Isosurface ",F8.2," 1 0 0 1 1")') dimplot
     write (luvmd,'("mol modcolor 0 top Volume 0")')
     write (luvmd,'("mol modmaterial 0 top Opaque")')
     write (luvmd,'("mol scaleminmax top 0 ",F8.2,X,F8.2)') -rhoplot/2d0*100d0, rhoplot/2d0*100d0
     write (luvmd,'("mol selupdate 0 top 0")')
     write (luvmd,'("mol colupdate 0 top 0")')
     write (luvmd,'("color scale method BGR")')
     write (luvmd,'("")')
     write (luvmd,'("# center")')
     write (luvmd,'("mol top 0")')
     write (luvmd,'("display resetview")')
     write (luvmd,'("rotate x by 15")')
     write (luvmd,'("rotate y by 20")')
     write (luvmd,'("rotate z by 3")')

     ! close files
     call fclose(ludat)
     call fclose(luvmd)

     if (.not.quiet) call tictac("End NCIPLOT")

   end subroutine nciplot

  !xx! private procedures

   subroutine write_cube_header(lu,l1,l2,periodic,nfrag,frag,x0x,x1x,x0,x1,nstep,xmat,isortho)
    use systemmod, only: sy
    use fragmentmod, only: fragment
    integer, intent(in) :: lu
    character*(*), intent(in) :: l1, l2
    logical, intent(in) :: periodic
    integer, intent(in) :: nfrag
    type(fragment), intent(in) :: frag(nfrag)
    real*8, intent(in) :: x0x(3), x1x(3), x0(3), x1(3)
    integer, intent(in) :: nstep(3)
    real*8, intent(in) :: xmat(3,3)
    logical :: isortho

    integer :: i, j, nc
    integer :: i1, i2, i3, nmin(3), nmax(3)
    real*8 :: xx(3)
    logical :: iszero

    real*8, parameter :: eps = 0.1d0

    if (periodic) then
       nc = sy%c%ncel
    else
       if (nfrag == 0) then
          nc = 0
          nmin = floor(x0x) - 1
          nmax = ceiling(x1x) + 1
          do i = 1, sy%c%ncel
             do i1 = nmin(1), nmax(1)
                do i2 = nmin(2), nmax(2)
                   do i3 = nmin(3), nmax(3)
                      xx = sy%c%x2c(sy%c%atcel(i)%x + real((/i1,i2,i3/),8))
                      if (all(xx > x0-eps) .and. all(xx < x1+eps)) then
                         nc = nc + 1
                      end if
                   end do
                end do
             end do
          end do
       else
          nc = 0
          do i = 1, nfrag
             nc = nc + frag(i)%nat
          end do
       end if
    end if

    iszero = (nc == 0)
    if (iszero) nc = sy%c%ncel

    write(lu,*) trim(l1)
    write(lu,*) trim(l2)
    write(lu,'(I5,3(F12.6))') nc, x0
    write(lu,'(I5,3(F12.6))') nstep(1), xmat(:,1)

    ! VMD screw-up log, entry #2:
    ! vmd detects whether a cube is orthogonal or not by the (1,2) element
    ! of the Cartesian matrix only. If this value is zero, vmd will misalign
    ! the isosurface and its box, the latter being displayed correctly. To
    ! prevent this, nudge the (1,2) element a little bit if the cube is not
    ! orthogonal and that element happens to be zero.
    if (isortho .or. abs(xmat(1,2)) > 2.d-6) then
       write(lu,'(I5,3(F12.6))') nstep(2), xmat(:,2)
    else
       write(lu,'(I5,3(F12.6))') nstep(2), 1.01d-6, xmat(2:3,2)
    end if
    write(lu,'(I5,3(F12.6))') nstep(3), xmat(:,3)

    if (periodic .or. iszero) then
       do i = 1, sy%c%ncel
          write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') sy%c%spc(sy%c%atcel(i)%is)%z, 0d0, sy%c%atcel(i)%r
       end do
    else
       if (nfrag == 0) then
          do i = 1, sy%c%ncel
             do i1 = nmin(1), nmax(1)
                do i2 = nmin(2), nmax(2)
                   do i3 = nmin(3), nmax(3)
                      xx = sy%c%x2c(sy%c%atcel(i)%x + real((/i1,i2,i3/),8))
                      if (all(xx > x0-eps) .and. all(xx < x1+eps)) then
                         write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') &
                            sy%c%spc(sy%c%atcel(i)%is)%z, 0d0, xx
                      end if
                   end do
                end do
             end do
          end do
       else
          do i = 1, nfrag
             do j = 1, frag(i)%nat
                write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') sy%c%spc(frag(i)%at(j)%is)%z, 0d0, frag(i)%at(j)%r
             end do
          end do
       end if
    end if

  end subroutine write_cube_header

  subroutine write_cube_body(lu,n,c)
    use tools_io, only: fclose

    integer, intent(in) :: lu
    integer, intent(in) :: n(3)
    real*8, intent(in) :: c(0:n(3)-1,0:n(2)-1,0:n(1)-1)

    integer :: i, j, k

    do i = 0, n(1)-1
       do j = 0, n(2)-1
          write (lu,'(6(1x,1p,e13.5e3))') (c(k,j,i),k=0,n(3)-1)
       enddo
    enddo
    call fclose(lu)

  end subroutine write_cube_body

  function check_no_extra_word(line,lp,routine)
    use tools_io, only: getword, ferror, faterr
    character(len=:), allocatable :: aux2
    character*(*), intent(in) :: line, routine
    integer, intent(inout) :: lp
    logical :: check_no_extra_word

    check_no_extra_word = .true.
    aux2 = getword(line,lp)
    if (len_trim(aux2) > 0) then
       call ferror(trim(routine),'Unknown extra keyword',faterr,line,syntax=.true.)
       check_no_extra_word = .false.
    end if
  end function check_no_extra_word

  function read_fragment(lu) result(fr)
    use systemmod, only: sy
    use global, only: eval_next
    use tools_io, only: uin, getline, lgetword, ucopy, equal, ferror, faterr, usegui
    use fragmentmod, only: fragment
    use types, only: realloc
    use param, only: bohrtoa, icrd_cart

    type(fragment) :: fr
    integer, intent(in) :: lu

    character(len=:), allocatable :: line, word
    integer :: lp, id
    logical :: ok
    real*8 :: x(3)

    fr%nat = 0
    allocate(fr%at(10))
    fr%nspc = sy%c%nspc
    allocate(fr%spc(fr%nspc))
    fr%spc = sy%c%spc

    ! create a fragment from input
    do while (.true.)
       if (lu == uin.and..not.usegui) then
          ok = getline(lu,line,.true.,ucopy)
       else
          ok = getline(lu,line,.true.)
       endif
       lp = 1
       word = lgetword(line,lp)
       if (.not.equal(word,"end") .and. .not.equal(word,"endfragment")) then
          lp = 1
          ok = eval_next(x(1),line,lp)
          ok = ok .and. eval_next(x(2),line,lp)
          ok = ok .and. eval_next(x(3),line,lp)
          if (.not.ok) call ferror('nciplot','bad atom in fragment',faterr)
          x = x / bohrtoa - sy%c%molx0
          id = sy%c%identify_atom(x,icrd_cart)

          if (id > 0) then
             fr%nat = fr%nat + 1
             if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
             fr%at(fr%nat)%r = x
             fr%at(fr%nat)%x = sy%c%c2x(x)
             fr%at(fr%nat)%cidx = id
             fr%at(fr%nat)%idx = sy%c%atcel(id)%idx
             fr%at(fr%nat)%lvec = nint(fr%at(fr%nat)%x - sy%c%atcel(id)%x)
             fr%at(fr%nat)%is = sy%c%atcel(id)%is
          end if
       else
          exit
       end if
    end do
    call realloc(fr%at,fr%nat)

  end function read_fragment

  !> Write the NCI checkpoint file. x0, xmat, and nstep are checks to
  !> make sure that the checkpoint file corresponds to this structure.
  !> crho, cgrad, rhoat, and rhofrag contain the density, rdg,
  !> promolecular, and fragment promolecular densities on output, and
  !> they should already be allocated when passed to this function.
  function readchk(oname,x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag) result(lchk)
    use tools_io, only: uout, fopen_read, fclose
    character*(*), intent(in) :: oname
    real*8, intent(in) :: x0(3), xmat(3,3)
    integer, intent(in) :: nstep(3)
    real*8, allocatable, dimension(:,:,:), intent(inout) :: crho, cgrad, rhoat
    integer, intent(in) :: nfrag
    real*8, allocatable, dimension(:,:,:,:), intent(inout) :: rhofrag
    logical :: lchk

    integer :: luchk
    character(len=:), allocatable :: chkfile
    logical :: llat, llfrag
    real*8 :: xxx(3)
    integer :: ixx(3), nfrag0

    chkfile = trim(oname) // ".chk_nci"

    inquire(file=chkfile,exist=lchk)

    if (lchk) then
       ! open the file
       lchk = .false.
       write(uout,'("  Reading the checkpoint file: ",A/)') trim(chkfile)
       luchk = fopen_read(chkfile,"unformatted")

       ! check that the structure and fragments are the same
       read (luchk) llat, llfrag, nfrag0
       if (.not.llat .and. allocated(rhoat)) goto 99
       if (.not.llfrag .and. allocated(rhofrag)) goto 99
       if (nfrag0 /= nfrag) goto 99
       read (luchk) xxx
       if (any(abs(xxx-x0) > 1d-13)) goto 99
       read (luchk) xxx
       if (any(abs(xxx-xmat(:,1)) > 1d-13)) goto 99
       read (luchk) xxx
       if (any(abs(xxx-xmat(:,2)) > 1d-13)) goto 99
       read (luchk) xxx
       if (any(abs(xxx-xmat(:,3)) > 1d-13)) goto 99
       read (luchk) ixx
       if (any(abs(ixx-nstep) > 0)) goto 99

       ! read the actual densities
       read (luchk) crho, cgrad
       if (llat .and. allocated(rhoat)) then
          read (luchk) rhoat
          if (allocated(rhofrag) .and. llfrag) then
             read (luchk) rhofrag
          endif
       endif

       ! wrap up
       lchk = .true.
       call fclose(luchk)
    end if

    return
99  continue
    call fclose(luchk)
    write(uout,'("  Checkpoint file is not compatible -- restarting."/)')
  end function readchk

  !> Write the checkpoint file.
  subroutine writechk(oname,x0,xmat,nstep,crho,cgrad,rhoat,nfrag,rhofrag)
    use tools_io, only: uout, fopen_write, fclose
    character*(*), intent(in) :: oname
    real*8, intent(in) :: x0(3), xmat(3,3)
    integer, intent(in) :: nstep(3)
    real*8, allocatable, dimension(:,:,:), intent(in) :: crho, cgrad, rhoat
    integer, intent(in) :: nfrag
    real*8, allocatable, dimension(:,:,:,:), intent(in) :: rhofrag

    character(len=:), allocatable :: chkfile

    integer :: luchk

    ! open the file
    chkfile = trim(oname) // ".chk_nci"
    write(uout,'("  Writing the checkpoint file: ",A/)') trim(chkfile)
    luchk = fopen_write(chkfile,"unformatted")

    ! write it
    write (luchk) allocated(rhoat), allocated(rhofrag), nfrag
    write (luchk) x0
    write (luchk) xmat(:,1)
    write (luchk) xmat(:,2)
    write (luchk) xmat(:,3)
    write (luchk) nstep
    write (luchk) crho, cgrad
    if (allocated(rhoat)) write (luchk) rhoat
    if (allocated(rhofrag)) write (luchk) rhofrag

    ! close it
    call fclose(luchk)

  end subroutine writechk

end submodule proc
