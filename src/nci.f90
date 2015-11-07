! Copyright (c) 2015 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
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
!
module nci
  implicit none

  private
  public :: nciplot
  private :: write_cube_header
  private :: write_cube_body
  private :: check_no_extra_word
  private :: read_fragment

contains

  ! nci plots
  subroutine nciplot()
    use fields
    use struct
    use struct_basic
    use global
    use grid_tools
    use grd_atomic
    use grid1_tools
    use graphics
    use fragmentmod
    use tools_io
    use tools_math
    use param
    use types

    type(field) :: fgrho, fxx(3)
    type(scalar_value) :: res, resg
    character(len=:), allocatable :: line, word, oname, file
    logical :: ok, ok2
    integer :: lp, istat
    integer :: i, j, k, iat, ifr, l, nn
    integer :: lugc, ludc, luvmd, ludat, luchk, lupc
    real(kind=gk), allocatable, dimension(:,:,:) :: crho, cgrad, rhoat
    real(kind=gk), allocatable, dimension(:,:,:,:) :: rhofrag
    real*8 :: ehess(3), dimgrad, dist, rrho, rrho1, rrho2
    real*8 :: sumq, sumqp, sumv
    logical :: onlyneg, lchk, inter, llat, llfrag, isden
    integer :: findlimits, ithres, istep
    logical :: periodic, usecore
    integer :: nfrag0, ixx(3)
    real*8 :: xxx(3)
    ! molmotif
    logical :: domolmotif
    ! chk
    logical :: usechk, dopromol
    ! fragments
    integer :: nfrag
    type(fragment), allocatable :: fr(:)
    type(fragment) :: fr0
    ! voids
    real*8 :: rho_void
    ! cube limits
    real*8 :: x(3), xd(3), x0(3), x1(3), x0x(3), x1x(3), xinc(3), xmat(3,3)
    real*8 :: xx0(3), xx1(3)
    integer :: nstep(3)
    ! cutoffs
    real*8 :: rhocut, dimcut, rhoplot, dimplot, rhoatl
    real*8, allocatable :: rhofragl(:)
    ! discarding rho parameter
    real*8 :: rhoparam, rhoparam2
    ! rthres
    real*8 :: rthres

    ! named constants
    real*8, parameter :: fthirds = 4d0/3d0
    real*8, parameter :: const = 2.d0*(3.d0*pi**2)**(1.d0/3.d0)
    real*8, parameter :: vsmall = 1d-40
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
    fgrho%init = .false.
    fgrho%numerical = .false.
    fxx(:)%init = .false.
    fxx(:)%numerical = .false.
    isden = .not.(refden == 0)
    rhoparam = 0.95d0
    rhoparam2 = 0.75d0
    x0 = 0d0
    x1 = 1d0
    x0 = cr%x2c(x0)
    x1 = cr%x2c(x1)
    periodic = .true.
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
    usecore = f(refden)%usecore .and. any(cr%at(1:cr%nneq)%zpsp /= -1)

    ! header
    write (uout,'("* NCIPLOT: non-covalent interactions")')
    if (.not.quiet) call tictac("start nciplot")

    do while(.true.)
       ok = getline(uin,line,.true.,ucopy)
       lp = 1
       word = lgetword(line,lp)

       if (equal(word,'oname')) then
          oname = line(lp:)
       elseif (equal(word,'molmotif')) then
          domolmotif = .true.
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'cutoffs')) then
          ok = eval_next(rhocut,line,lp)
          ok = ok .and. eval_next(dimcut,line,lp)
          if (.not.ok) call ferror('nciplot','wrong cutoffs keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'rhoparam')) then
          ok = eval_next(rhoparam,line,lp)
          if (.not.ok) call ferror('nciplot','wrong rhoparam keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'rhoparam2')) then
          ok = eval_next(rhoparam2,line,lp)
          if (.not.ok) call ferror('nciplot','wrong rhoparam2 keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'cutplot')) then
          ok = eval_next(rhoplot,line,lp)
          ok = ok .and. eval_next(dimplot,line,lp)
          if (.not.ok) call ferror('nciplot','wrong cutplot keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'void')) then
          ok = eval_next(rho_void,line,lp)
          if (.not.ok) call ferror('nciplot','wrong void keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'rthres')) then
          ok = eval_next(rthres,line,lp)
          if (.not.ok) call ferror('nciplot','wrong cutplot keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'increments')) then
          istep = 0
          ok = eval_next(xinc(1),line,lp)
          ok = ok .and. eval_next(xinc(2),line,lp)
          ok = ok .and. eval_next(xinc(3),line,lp)
          if (.not.ok) call ferror('nciplot','wrong increments keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'nstep')) then
          istep = 1
          ok = eval_next(nstep(1),line,lp)
          ok = ok .and. eval_next(nstep(2),line,lp)
          ok = ok .and. eval_next(nstep(3),line,lp)
          if (.not.ok) call ferror('nciplot','wrong nstep keyword',faterr,line)
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'onlyneg')) then
          onlyneg = .true.
          call check_no_extra_word(line,lp,'nciplot')
       elseif (equal(word,'nochk')) then
          usechk = .false.
          call check_no_extra_word(line,lp,'nciplot')
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
             call check_no_extra_word(line,lp,'nciplot')
             x0 = cr%x2c(x0)
             x1 = cr%x2c(x1)
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
                fr0 = cr%identify_fragment_from_xyz(word)
                do j = 1, fr0%nat
                   do k = 1, 3
                      x0(k) = min(x0(k),fr0%at(j)%r(k))
                      x1(k) = max(x1(k),fr0%at(j)%r(k))
                   end do
                end do
             end do
             ithres = 1
          endif
          if (.not.ok) call ferror('nciplot','wrong cube keyword',faterr,line)
       elseif (equal(word,'fragment')) then
          if (findlimits == -1) findlimits = 1
          if (ithres == -1) ithres = 1
          nfrag = nfrag + 1
          if (nfrag > size(fr)) call realloc(fr,2*nfrag)
          word = getword(line,lp)
          if (equal(word,'')) then
             fr(nfrag) = read_fragment(uin)
          else
             fr(nfrag) = cr%identify_fragment_from_xyz(word)
          endif

       elseif (equal(word,'endnciplot').or.equal(word,'end')) then
          call check_no_extra_word(line,lp,'nciplot')
          exit
       else
          call ferror ('nciplot','unrecognized option',faterr,line)
       end if
    end do
    if (nfrag > 0) call realloc(fr,nfrag)

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
    x0 = cr%c2x(x0)
    x1 = cr%c2x(x1)
    x0x = x0
    x1x = x1

    ! determine the grid step
    if (periodic) then
       ! use xinc for periodic cells
       do i = 1, 3
          x1 = 0d0
          x1(i) = 1d0
          x1 = cr%x2c(x1)
          if (istep == -1 .and. f(refden)%type == type_grid) then
             nstep = f(refden)%n
          elseif (istep == 0 .or. istep == -1) then
             nstep(i) = ceiling(sqrt(dot_product(x1,x1)) / xinc(i))
          end if
          xmat(:,i) = x1 / real(nstep(i),8)
          xinc(i) = sqrt(dot_product(xmat(:,i),xmat(:,i)))
       end do
    else
       ! use the grid step
       x0 = cr%x2c(x0)
       x1 = cr%x2c(x1)
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
    write(uout,'("+ Density cube file: ",A)') trim(oname) // "-dens.cube"
    write(uout,'("+ RDG cube file: ",A)') trim(oname) // "-grad.cube"
    if (dopromol) then
       write(uout,'("+ Promolecular density cube file: ",A)') trim(oname) // "-pdens.cube"
    end if
    write(uout,'("+ RDG vs. density file: ",A)') trim(oname) // ".dat"
    write(uout,'("+ VMD visualization script: ",A)') trim(oname) // ".vmd"
    write(uout,'("+ Density cutoff (for dat file): ",A)') string(rhocut,'f',decimal=6)
    write(uout,'("+ RDG cutoff (for dat file): ",A)') string(dimcut,'f',decimal=6)
    write(uout,'("+ Density plot cutoff: ",A)') string(rhoplot,'f',decimal=6)
    write(uout,'("+ RDG plot cutoff: ",A)') string(dimplot,'f',decimal=6)
    write(uout,'("+ Cube origin (crystal coord.): ",3(A,X))') (string(x0x(j),'f',decimal=6),j=1,3)
    write(uout,'("+ Cube end (crystal coord.): ",3(A,X))') (string(x1x(j),'f',decimal=6),j=1,3)
    write(uout,'("+ Cube nodes in each direction: ",3(A,X))') (string(nstep(j)),j=1,3)
    write(uout,'("+ Cube steps in each direction: ",3(A,X))') (string(xinc(j),'f',decimal=6),j=1,3)
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
    call write_cube_header(lugc,'grad_cube','3d plot, reduced density gradient',periodic,nfrag,fr,x0x,x1x,x0,x1,nstep,xmat)
    call write_cube_header(ludc,'dens_cube','3d plot, density',periodic,nfrag,fr,x0x,x1x,x0,x1,nstep,xmat)
    if (dopromol) then
       call write_cube_header(lupc,'pdens_cube','3d plot, promolecular density',periodic,nfrag,fr,x0x,x1x,x0,x1,nstep,xmat)
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

     if (f(0)%type /= type_promol) then
        call ferror('nciplot','promolecular density not found',faterr)
     end if

     ! calculate density, rdg,... read from chkpoint if available
     if (usechk) then
        inquire(file=trim(oname)//".ncichk",exist=lchk)
     else
        lchk = .false.
     endif
99   continue
     if (lchk) then
        write(uout,'("  Reading the checkpoint file: ",A/)') trim(oname)//".ncichk"
        luchk = fopen_read(trim(oname)//".ncichk","unformatted")
        read (luchk) llat, llfrag, nfrag0
        if (.not.llat .and. allocated(rhoat)) goto 999
        if (.not.llfrag .and. allocated(rhofrag)) goto 999
        if (nfrag0 /= nfrag) goto 999
        read (luchk) xxx
        if (any(abs(xxx-x0) > 1d-13)) goto 999
        read (luchk) xxx
        if (any(abs(xxx-xmat(:,1)) > 1d-13)) goto 999
        read (luchk) xxx
        if (any(abs(xxx-xmat(:,2)) > 1d-13)) goto 999
        read (luchk) xxx
        if (any(abs(xxx-xmat(:,3)) > 1d-13)) goto 999
        read (luchk) ixx
        if (any(abs(ixx-nstep) > 0)) goto 999
        read (luchk) crho, cgrad
        if (llat .and. allocated(rhoat)) then
           read (luchk) rhoat
           if (allocated(rhofrag) .and. llfrag) then
              read (luchk) rhofrag
           endif
        endif
        lchk = .false.

999     continue
        call fclose(luchk)
        if (lchk) then
           write(uout,'("  Checkpoint file is not compatible -- restarting."/)')
           lchk = .false.
           goto 99
        endif
     else
        ! Initialize gradrho out of the omp loop
        if (f(refden)%type == type_grid .and..not.usecore) then
           call grid_gradrho(f(refden),fgrho)
           do i = 1, 3
              call grid_hxx(f(refden),fxx(i),i)
           end do
        endif

        ! calculate density, rdg, rhoat
        allocate(rhofragl(nfrag))
        rhofragl = 0d0
        !$omp parallel do private (x,res,resg,ehess,dimgrad,rhoatl,&
        !$omp xd,dist,rrho,rrho1,rrho2) firstprivate(rhofragl) schedule(dynamic)
        do i = 0, nstep(1)-1
           do j = 0, nstep(2)-1
              do k = 0, nstep(3)-1
                 x = x0 + i*xmat(:,1) + j*xmat(:,2) + k*xmat(:,3)

                 ! calculate properties at x: rho and rdg
                 if (f(refden)%type == type_grid.and..not.usecore) then
                    call grd(f(refden),x,0,res)
                    call grd(fgrho,x,0,resg)
                    dimgrad = resg%f / (const*max(res%f,vsmall)**fthirds)
                    do l = 1, 3
                       call grd(fxx(l),x,0,resg)
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
                    call grd(f(refden),x,2,res)
                    call eig(res%hf,ehess)
                    dimgrad = res%gfmod / (const*max(res%f,vsmall)**fthirds)           
                 end if

                 ! promolecular density
                 if (dopromol) then
                    rhoatl = grd0(f(0),x)
                 endif

                 ! fragments
                 if (nfrag > 0) then
                    rhofragl = 0d0
                    do ifr = 1, nfrag
                       do iat = 1, fr(ifr)%nat
                          xd = x - fr(ifr)%at(iat)%r
                          dist = sqrt(dot_product(xd,xd))
                          if (dist > agrid(fr(ifr)%at(iat)%z)%rmax) cycle
                          dist = max(dist,agrid(fr(ifr)%at(iat)%z)%r(1))
                          dist = max(dist,1d-14)
                          call grid1_interp(agrid(fr(ifr)%at(iat)%z),dist,rrho,rrho1,rrho2)
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
        if (usechk) then
           write(uout,'("  Writing the checkpoint file: ",A/)') trim(oname)//".ncichk"
           luchk = fopen_write(trim(oname)//".ncichk","unformatted")
           write (luchk) allocated(rhoat), allocated(rhofrag), nfrag
           write (luchk) x0
           write (luchk) xmat(:,1)
           write (luchk) xmat(:,2)
           write (luchk) xmat(:,3)
           write (luchk) nstep
           write (luchk) crho, cgrad
           if (allocated(rhoat)) write (luchk) rhoat
           if (allocated(rhofrag)) write (luchk) rhofrag
           call fclose(luchk)
        endif
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
              if (crho(k,j,i) > 0 .and. onlyneg) cgrad(k,j,i) = 100d0
           end do
        end do
     end do
     sumq = sumq * cr%omega / (nstep(1)*nstep(2)*nstep(3))
     sumqp = sumqp * cr%omega / (nstep(1)*nstep(2)*nstep(3))
     sumv = sumv * cr%omega / (nstep(1)*nstep(2)*nstep(3))

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

     ! write xyz file
     file = trim(oname) // "_cell.xyz"
     if (nfrag == 0) then
        if (periodic) then
           file = trim(file) // " border"
           if (domolmotif) file = trim(file) // " molmotif"
           call struct_write(file)
        else
           call fragment_init(fr0)
           xx0 = cr%c2x(x0 - rthres_xyz)
           xx1 = cr%c2x(x1 + rthres_xyz)
           do i = floor(xx0(1))-1, ceiling(xx1(1))+1
              do j = floor(xx0(2))-1, ceiling(xx1(2))+1
                 do k = floor(xx0(3))-1, ceiling(xx1(3))+1
                    do l = 1, cr%ncel
                       x = cr%atcel(l)%x + real((/i,j,k/),8)
                       if (x(1) > xx0(1) .and. x(1) < xx1(1) .and.&
                           x(2) > xx0(2) .and. x(2) < xx1(2) .and.&
                           x(3) > xx0(3) .and. x(3) < xx1(3)) then
                          fr0%nat = fr0%nat + 1
                          if (fr0%nat > size(fr0%at)) call realloc(fr0%at,2*fr0%nat)
                          fr0%at(fr0%nat)%x = x
                          fr0%at(fr0%nat)%r = cr%x2c(x)
                          fr0%at(fr0%nat)%cidx = l
                          fr0%at(fr0%nat)%idx = cr%atcel(l)%idx
                          fr0%at(fr0%nat)%lvec = (/i,j,k/)
                          fr0%at(fr0%nat)%z = cr%at(cr%atcel(l)%idx)%z
                       end if
                    end do
                 end do
              end do
           end do
           call writexyz(file,fr0)
        end if
     else
        call writexyz(file,fragment_merge_array(fr(1:nfrag)))
     end if

     ! ! write vmd script
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
           xx0 = cr%x2c(xlist0(:,i))
           xx1 = cr%x2c(xlist1(:,i))
           write (luvmd,'("draw cylinder {",3(F14.4,X),"} {",3(F14.4,X),"} radius 0.05 resolution 20")') &
              xx0 * bohrtoa, xx1 * bohrtoa
        end do
        write (luvmd,'("")')
     end if
     write (luvmd,'("# Load the cubes")')
     write (luvmd,'("mol new ",A," type cube first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all")') trim(oname)//"-dens.cube"
     write (luvmd,'("mol addfile ",A," type cube first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all")') trim(oname)//"-grad.cube"
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
     if (usechk) call fclose(luchk)
     call fclose(luvmd)

     if (.not.quiet) call tictac("End NCIPLOT")

  end subroutine nciplot

  subroutine write_cube_header(lu,l1,l2,periodic,nfrag,frag,x0x,x1x,x0,x1,nstep,xmat)
    use struct_basic
    use types

    integer, intent(in) :: lu
    character*(*), intent(in) :: l1, l2
    logical, intent(in) :: periodic
    integer, intent(in) :: nfrag
    type(fragment), intent(in) :: frag(nfrag)
    real*8, intent(in) :: x0x(3), x1x(3), x0(3), x1(3)
    integer, intent(in) :: nstep(3)
    real*8, intent(in) :: xmat(3,3)

    integer :: i, j, nc
    integer :: i1, i2, i3, nmin(3), nmax(3)
    real*8 :: xx(3)
    logical :: iszero

    real*8, parameter :: eps = 0.1d0

    if (periodic) then
       nc = cr%ncel
    else
       if (nfrag == 0) then
          nc = 0
          nmin = floor(x0x) - 1
          nmax = ceiling(x1x) + 1
          do i = 1, cr%ncel
             do i1 = nmin(1), nmax(1)
                do i2 = nmin(2), nmax(2)
                   do i3 = nmin(3), nmax(3)
                      xx = cr%x2c(cr%atcel(i)%x + real((/i1,i2,i3/),8))
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
    if (iszero) nc = cr%ncel

    write(lu,*) trim(l1)
    write(lu,*) trim(l2)
    write(lu,'(I5,3(F12.6))') nc, x0
    write(lu,'(I5,3(F12.6))') nstep(1), xmat(:,1)
    write(lu,'(I5,3(F12.6))') nstep(2), xmat(:,2)
    write(lu,'(I5,3(F12.6))') nstep(3), xmat(:,3)
    if (periodic .or. iszero) then
       do i = 1, cr%ncel
          write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') cr%at(cr%atcel(i)%idx)%z, 0d0, cr%atcel(i)%r
       end do
    else
       if (nfrag == 0) then
          do i = 1, cr%ncel
             do i1 = nmin(1), nmax(1)
                do i2 = nmin(2), nmax(2)
                   do i3 = nmin(3), nmax(3)
                      xx = cr%x2c(cr%atcel(i)%x + real((/i1,i2,i3/),8))
                      if (all(xx > x0-eps) .and. all(xx < x1+eps)) then
                         write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') &
                            cr%at(cr%atcel(i)%idx)%z, 0d0, xx
                      end if
                   end do
                end do
             end do
          end do
       else
          do i = 1, nfrag
             do j = 1, frag(i)%nat
                write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') cr%at(frag(i)%at(j)%idx)%z, 0d0, frag(i)%at(j)%r
             end do
          end do
       end if
    end if

  end subroutine write_cube_header

  subroutine write_cube_body(lu,n,c)
    use tools_io
    use types
    use param

    integer, intent(in) :: lu
    integer, intent(in) :: n(3)
    real(kind=gk), intent(in) :: c(0:n(3)-1,0:n(2)-1,0:n(1)-1)

    integer :: i, j, k

    do i = 0, n(1)-1
       do j = 0, n(2)-1
          write (lu,'(6(1x,e12.5))') (c(k,j,i),k=0,n(3)-1)
       enddo
    enddo
    call fclose(lu)

  end subroutine write_cube_body

  subroutine check_no_extra_word(line,lp,routine)
    use tools_io
    character(len=:), allocatable :: aux2
    character*(*), intent(in) :: line, routine
    integer, intent(inout) :: lp

    aux2 = getword(line,lp)
    if (len_trim(aux2) > 0) &
       call ferror(trim(routine),'Unknown extra keyword',faterr,line)

  end subroutine check_no_extra_word

  function read_fragment(lu) result(fr)
    use struct_basic
    use varbas
    use global
    use tools_io
    use types
    use param

    type(fragment) :: fr
    integer, intent(in) :: lu

    character(len=:), allocatable :: line, word
    integer :: lp, id
    logical :: ok
    real*8 :: x(3)

    fr%nat = 0
    allocate(fr%at(10))

    ! create a fragment from input
    do while (.true.)
       if (lu == uin) then
          ok = getline(lu,line,.true.,ucopy)
       else
          ok = getline(lu,line,.true.)
       endif
       lp = 1
       word = lgetword(line,lp)
       if (.not.equal(word,"end") .and. .not.equal(word,"endfragment")) then
          fr%nat = fr%nat + 1
          if (fr%nat > size(fr%at)) call realloc(fr%at,2*fr%nat)
          lp = 1
          ok = eval_next(x(1),line,lp)
          ok = ok .and. eval_next(x(2),line,lp)
          ok = ok .and. eval_next(x(3),line,lp)
          if (.not.ok) call ferror('nciplot','bad atom in fragment',faterr)
          x = x / bohrtoa
          id = cr%identify_atom(x,.true.)
          fr%at(fr%nat)%r = x
          fr%at(fr%nat)%x = cr%c2x(x)
          fr%at(fr%nat)%cidx = id
          fr%at(fr%nat)%idx = cr%atcel(id)%idx
          fr%at(fr%nat)%lvec = nint(fr%at(fr%nat)%x - cr%atcel(id)%x)
          fr%at(fr%nat)%z = cr%at(fr%at(fr%nat)%idx)%z
       else
          exit
       end if
    end do
    word = getword(line,lp)
    if (len_trim(word) > 0) &
       call ferror('read_fragment','Unknown extra keyword',faterr,line)
    call realloc(fr%at,fr%nat)
    
  end function read_fragment

end module nci
