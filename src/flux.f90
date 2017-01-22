! Copyright (c) 2015 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Code for the FLUXPRINT environment: 3d plotting.
module flux
  implicit none

  private

  !xx! list of routines
  public  :: fluxprint
  private :: flx_initialize
  public  :: flx_end
  private :: flx_printpath
  private :: flx_point
  private :: flx_ncp
  private :: flx_bcp
  private :: flx_graph
  private :: flx_findthetagrid
  private :: ball_radius

  logical :: flx_init = .false. !< is fluxprint initialized?

  ! General info
  integer :: flx_n
  integer :: flx_iup
  integer, parameter :: flx_m = 4000

  ! output format
  character*3 :: outfmt
  integer :: luout, lumtl
  character(len=:), allocatable :: outfile

  ! Gradient path info (cryst.)
  real*8, dimension(:,:), allocatable :: flx_x
  real*8, dimension(:), allocatable :: flx_rho
  real*8, dimension(:,:), allocatable :: flx_grad
  real*8, dimension(:,:,:), allocatable :: flx_h

contains

  subroutine fluxprint()
    use varbas, only: ncpcel, cpcel, cp
    use struct_basic, only: cr
    use global, only: eval_next, dunit
    use tools_io, only: uout, getline, lgetword, equal, ferror, faterr, uin, &
       ucopy, getword, lower
    use param, only: jmlcol
    
    type flxorder
       character*3 :: id, method
       real*8 :: x(3)
       integer :: iup, cpid, ntheta, nphi, lvec(3), num
       integer :: rgb(3)
    end type flxorder
    integer :: norder
    type(flxorder), allocatable :: order(:), orderaux(:)
    character*3 :: method
    character(len=:), allocatable :: line, word
    integer :: local_flx_nsym
    integer :: cpid, lvecx(3), iup, lp, nn, nphi, ntheta
    integer :: num_points, i, rgb(3)
    real*8 :: xpoint(3)
    logical :: ok, nosym

    integer, parameter :: irgb(3) = (/255,179,77/)

    ! set initial values for flx
    flx_init = .false.
    local_flx_nsym = -1
    nosym = .false.
    outfmt = "cml"
    allocate(order(10))
    norder = 0
    rgb = -1

    ! header
    write (uout,'(A)') "* FLUXPRINT: 3d representation of cell, CPs and gradient paths"
    write (uout,*)
    do while(.true.)
       ok = getline(uin,line,eofstop=.true.,ucopy=ucopy)
       lp = 1
       word = lgetword(line,lp)

       if (equal(word,'nosym')) then
          nosym = .true.
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'color')) then
          ok=eval_next(rgb(1),line,lp)
          ok= ok .and. eval_next(rgb(2),line,lp)
          ok= ok .and. eval_next(rgb(3),line,lp)
          if (.not. ok) then
             call ferror ('fluxprint','wrong COLOR syntax',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'text')) then
          outfmt = "txt"
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'tessel').or.equal(word,'tess')) then
          outfmt = "tss"
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'obj')) then
          outfmt = "obj"
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'ply')) then
          outfmt = "ply"
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'off')) then
          outfmt = "off"
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'cml')) then
          outfmt = "cml"
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'shells')) then
          ok = eval_next(local_flx_nsym,line,lp)
          if (.not. ok) then
             call ferror ('fluxprint','wrong SHELLS syntax',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return
       else if (equal(word,'point')) then
          ok=eval_next(iup ,line,lp)
          ok= ok .and. eval_next(xpoint(1),line,lp)
          ok= ok .and. eval_next(xpoint(2),line,lp)
          ok= ok .and. eval_next(xpoint(3),line,lp)
          if (.not. ok) then
             call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return
          if (cr%ismolecule) &
             xpoint = cr%c2x(xpoint / dunit - cr%molx0)

          if (iup /= 1 .and. iup /= -1 .and. iup /= 0) then
             call ferror('fluxprint','iup must be +1, 0 or -1',faterr,line,syntax=.true.)
             return
          end if

          norder = norder + 1
          if (norder > size(order)) call realloc_flxorder()
          order(norder)%id = "poi"
          order(norder)%iup = iup
          order(norder)%x = xpoint
          if (any(rgb < 0)) then
             order(norder)%rgb = irgb
          else
             order(norder)%rgb = rgb
          end if
       elseif (equal(word,'ncp')) then
          ok=eval_next(cpid,line,lp)
          if (cpid <= 0 .or. cpid > ncpcel .or. (cpcel(cpid)%typ /= -3)) then
             call ferror ('fluxprint','ncp identifier not known.',faterr,line,syntax=.true.)
             return
          end if
          ok=ok.and.eval_next(ntheta,line,lp)
          ok=ok.and.eval_next(nphi,line,lp)
          if (.not. ok) then
             call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
             return
          end if

          lvecx = 0
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'lvec')) then
                ok=ok.and.eval_next(lvecx(1),line,lp)
                ok=ok.and.eval_next(lvecx(2),line,lp)
                ok=ok.and.eval_next(lvecx(3),line,lp)
                if (.not. ok) then
                   call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
                   return
                end if
             elseif (len_trim(word) > 0) then
                call ferror('fluxprint','Unknown extra keyword',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do

          norder = norder + 1
          if (norder > size(order)) call realloc_flxorder()
          order(norder)%id = "ncp"
          order(norder)%cpid = cpid
          order(norder)%ntheta = ntheta
          order(norder)%nphi = nphi
          order(norder)%lvec = lvecx
          if (any(rgb < 0).and.cp(cpid)%isnuc) then
             order(norder)%rgb = jmlcol(:,cr%at(cr%atcel(cpid)%idx)%z)
          else
             order(norder)%rgb = rgb
          end if

       elseif (equal(word,'bcp')) then
          ok=eval_next(cpid,line,lp)
          if (cpid < 0 .or. cpid > ncpcel .or. (cpcel(cpid)%typ /= -1) ) then
             call ferror ('fluxprint','bcp identifier not recognized.',faterr,line,syntax=.true.)
             return
          end if

          ok=ok.and.eval_next(iup,line,lp)
          if (iup /= 1 .and. iup /= -1 .and. iup /= 0) then
             call ferror('fluxprint','iup must be +1, 0 or -1',faterr,line,syntax=.true.)
             return
          end if
          if (iup /= 1) then
             ok=ok.and.eval_next(num_points,line,lp)
          else
             num_points = 0
          end if
          if (.not. ok) then
             call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
             return
          end if

          method = "bra"
          lvecx = 0
          do while (.true.)
             word = getword(line,lp)
             word = lower(word)
             if (equal(word,'lvec')) then
                ok=ok.and.eval_next(lvecx(1),line,lp)
                ok=ok.and.eval_next(lvecx(2),line,lp)
                ok=ok.and.eval_next(lvecx(3),line,lp)
                if (.not. ok) then
                   call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
                   return
                end if
             else if (equal(word,'braindead')) then
                method = "bra"
             else if (equal(word,'quotient')) then
                method = "quo"
             else if (equal(word,'dynamical')) then
                method = "dyn"
             else if (equal(word,'h1')) then
                method = "h1 "
             else if (len_trim(word) > 0) then
                call ferror('fluxprint','Unknown extra keyword',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do

          norder = norder + 1
          if (norder > size(order)) call realloc_flxorder()
          order(norder)%id = "bcp"
          order(norder)%iup = iup
          order(norder)%cpid = cpid
          order(norder)%num = num_points
          order(norder)%method = method
          order(norder)%lvec = lvecx
          if (any(rgb < 0)) then
             order(norder)%rgb = irgb
          else
             order(norder)%rgb = rgb
          end if

       elseif (equal(word,'rcp')) then
          ok=eval_next(cpid,line,lp)
          if (cpid < 0 .or. cpid > ncpcel .or. (cpcel(cpid)%typ /= 1) ) then
             call ferror ('fluxprint','rcp identifier not recognized.',faterr,line,syntax=.true.)
             return
          end if

          ok=ok.and.eval_next(iup,line,lp)
          if (iup /= 1 .and. iup /= -1 .and. iup /= 0) then
             call ferror('fluxprint','iup must be +1, 0 or -1',faterr,line,syntax=.true.)
             return
          end if
          if (iup /= -1) then
             ok=ok.and.eval_next(num_points,line,lp)
          else
             num_points = 0
          end if
          if (.not. ok) then
             call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
             return
          end if

          method = "bra"
          lvecx = 0
          do while (.true.)
             word = getword(line,lp)
             word = lower(word)
             if (equal(word,'lvec')) then
                ok=ok.and.eval_next(lvecx(1),line,lp)
                ok=ok.and.eval_next(lvecx(2),line,lp)
                ok=ok.and.eval_next(lvecx(3),line,lp)
                if (.not. ok) then
                   call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
                   return
                end if
             else if (equal(word,'braindead')) then
                method = "bra"
             else if (equal(word,'quotient')) then
                method = "quo"
             else if (equal(word,'dynamical')) then
                method = "dyn"
             else if (len_trim(word) > 0) then
                call ferror('fluxprint','Unknown extra keyword',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do

          norder = norder + 1
          if (norder > size(order)) call realloc_flxorder()
          order(norder)%id = "rcp"
          order(norder)%iup = iup
          order(norder)%cpid = cpid
          order(norder)%num = num_points
          order(norder)%method = method
          order(norder)%lvec = lvecx
          if (any(rgb < 0)) then
             order(norder)%rgb = irgb
          else
             order(norder)%rgb = rgb
          end if

       elseif (equal(word,'ccp')) then
          ok = eval_next(cpid,line,lp)
          if (cpid < 0 .or. cpid > ncpcel .or. (cpcel(cpid)%typ /= 3)) then
             call ferror ('fluxprint','ccp identifier not recognized.',faterr,line,syntax=.true.)
             return
          end if
          ok=ok.and.eval_next(ntheta,line,lp)
          ok=ok.and.eval_next(nphi,line,lp)
          if (.not. ok) then
             call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
             return
          end if

          lvecx = 0
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,'lvec')) then
                ok=ok.and.eval_next(lvecx(1),line,lp)
                ok=ok.and.eval_next(lvecx(2),line,lp)
                ok=ok.and.eval_next(lvecx(3),line,lp)
                if (.not. ok) then
                   call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
                   return
                end if
             elseif (len_trim(word) > 0) then
                call ferror('fluxprint','Unknown extra keyword',faterr,line,syntax=.true.)
                return
             else
                exit
             end if
          end do

          norder = norder + 1
          if (norder > size(order)) call realloc_flxorder()
          order(norder)%id = "ccp"
          order(norder)%cpid = cpid
          order(norder)%ntheta = ntheta
          order(norder)%nphi = nphi
          order(norder)%lvec = lvecx
          if (any(rgb < 0)) then
             order(norder)%rgb = irgb
          else
             order(norder)%rgb = rgb
          end if

       elseif (equal(word,'graph')) then
          ok = eval_next(nn,line,lp)
          if (.not. ok) then
             call ferror ('fluxprint','Wrong syntax',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word()
          if (.not.ok) return

          norder = norder + 1
          if (norder > size(order)) call realloc_flxorder()
          order(norder)%id = "gra"
          order(norder)%num = nn
          if (any(rgb < 0)) then
             order(norder)%rgb = irgb
          else
             order(norder)%rgb = rgb
          end if

       elseif (equal(word,'end').or.equal(word,'endfluxprint')) then
          ok = check_no_extra_word()
          if (.not.ok) return
          exit
       else
          call ferror('fluxprint','Unknown keyword',faterr,line,syntax=.true.)
          return
       end if
    end do

    ! run the commands, in order
    call flx_initialize(nosym)
    do i = 1, norder
       if (order(i)%id == "poi") then
          call flx_point(order(i)%x,order(i)%iup,local_flx_nsym,order(i)%rgb)
       elseif (order(i)%id == "ncp") then
          call flx_ncp(order(i)%cpid,order(i)%ntheta,order(i)%nphi,local_flx_nsym,order(i)%lvec,order(i)%rgb)
       elseif (order(i)%id == "bcp") then
          call flx_bcp(order(i)%cpid,order(i)%iup,order(i)%num,local_flx_nsym,order(i)%method,order(i)%lvec,order(i)%rgb)
       elseif (order(i)%id == "rcp") then
          call flx_bcp(order(i)%cpid,order(i)%iup,order(i)%num,local_flx_nsym,order(i)%method,order(i)%lvec,order(i)%rgb)
       elseif (order(i)%id == "ccp") then
          call flx_ncp(order(i)%cpid,order(i)%ntheta,order(i)%nphi,local_flx_nsym,order(i)%lvec,order(i)%rgb)
       elseif (order(i)%id == "gra") then
          call flx_graph(local_flx_nsym,order(i)%num,rgb=order(i)%rgb)
       else
          call ferror('fluxprint','Unknown order id',faterr)
       end if
    end do
    call flx_end(.false.)
    deallocate(order)

  contains
    
    function check_no_extra_word()
      character(len=:), allocatable :: aux2
      logical :: check_no_extra_word
      aux2 = getword(line,lp)
      check_no_extra_word = .true.
      if (len_trim(aux2) > 0) then
         call ferror('fluxprint','Unknown extra keyword',faterr,line)
         check_no_extra_word = .false.
      end if
    end function check_no_extra_word

    subroutine realloc_flxorder()
      allocate(orderaux(2*norder))
      orderaux(1:norder-1) = order(1:norder-1)
      call move_alloc(orderaux,order)
    end subroutine realloc_flxorder

  end subroutine fluxprint

  !> Initialize the module, allocate memory for a gradient
  !> path. prune, plot gradient path, one point each prune distance
  !> step along the path. If noballs, do not write the BALL commands
  !> to the tessel output (optional). If nocell, do not write the CELL
  !> command to the tessel output (optional). If nosym, do not write
  !> the symmetry elements.
  subroutine flx_initialize(nosym,noballs,nocell)
    use varbas, only: ncpcel, cpcel, ncp, cp
    use global, only: fileroot
    use struct_basic, only: cr
    use struct_writers, only: struct_write_3dmodel, struct_write_mol
    use tools_io, only: faterr, ferror, uout, fopen_write

    logical, intent(in) :: nosym
    logical, intent(in), optional :: noballs, nocell

    integer :: i
    character(2) :: label
    logical :: doballs, docell

    if (flx_init) then
       call ferror('flx_initialize','The flux module is already initialized.',faterr)
    end if

    ! Initialize variables
    flx_n = 0
    flx_init = .true.

    doballs = .true.
    if (present(noballs)) then
       if (noballs) doballs = .false.
    end if
    docell = .true.
    if (present(nocell)) then
       if (nocell) docell = .false.
    end if

    ! Allocate memory
    allocate(flx_x(3,flx_m))
    allocate(flx_rho(flx_m))
    allocate(flx_grad(3,flx_m))
    allocate(flx_h(3,3,flx_m))

    ! Connect and initialize tessel file
    if (outfmt == "tss") then
       outfile = trim(fileroot) // "_flux.tess"
       write (uout,'(A,A/)') "* Writing tessel file : ", trim(outfile)
       luout = fopen_write(outfile)

       write (luout,'(A)') "set camangle 75 -10 45"
       write (luout,'(A)') "set background background {color rgb <1,1,1>}"
       write (luout,'(A)') "set use_planes .false."
       write (luout,'(A)') "set ball_texture finish{specular 0.2 roughness 0.1 reflection 0.1}"
       write (luout,'(A)') "set equalscale noscale"
       ! molecule
       write (luout,'(A)') "molecule"
       ! crystal
       write (luout,'(X,A)') "crystal"
       !! general info
       write (luout,'(2X,A)') "title gradient path representation."
       write (luout,'(2X,A)') "symmatrix seitz"
       if (nosym) then
          write (luout,'(3X,3(F5.2,X),F15.12)') 1d0, 0d0, 0d0, 0d0
          write (luout,'(3X,3(F5.2,X),F15.12)') 0d0, 1d0, 0d0, 0d0
          write (luout,'(3X,3(F5.2,X),F15.12)') 0d0, 0d0, 1d0, 0d0
       else
          do i = 1, cr%ncv
             write (luout,'(3X,A,3(F15.12,X))') "cen ", cr%cen(:,i)
          end do
          write (luout,'(3X,A)') "#"
          do i = 1, cr%neqv
             write (luout,'(3X,3(F5.2,X),F15.12)') cr%rotm(1,:,i)
             write (luout,'(3X,3(F5.2,X),F15.12)') cr%rotm(2,:,i)
             write (luout,'(3X,3(F5.2,X),F15.12)') cr%rotm(3,:,i)
             write (luout,'(3X,A)') "#"
          end do
       end if
       write (luout,'(2X,A)') "endsymmatrix"
       write (luout,'(2X,A,6(F10.6" "))') "cell", cr%aa, cr%bb
       write (luout,'(2X,A)') "crystalbox  -2.30 -2.30 -2.30 2.30 2.30 2.30"
       write (luout,'(2X,A,6(F6.3,X))') "clippingbox ",-0.02,-0.02,-0.02,+1.02,+1.02,+1.02
       if (nosym) then
          do i = 1,ncpcel
             if (cpcel(i)%idx <= cr%nneq) then
                label = trim(cr%at(cpcel(i)%idx)%name)
                if (label(2:2) == " ") label(2:2) = "_"
             else if (cpcel(i)%typ == -3) then
                label = "XX"
             else if (cpcel(i)%typ == -1)  then
                label = "YY"
             else if (cpcel(i)%typ == 1) then
                label = "ZZ"
             else
                label = "XZ"
             end if
             write (luout,'(2X,A,3(F10.6," "),A2,I4.4,a)') "neq ",cpcel(i)%x,label,i," 0"
          end do
       else
          do i = 1, ncp
             if (i <= cr%nneq) then
                label = trim(cr%at(i)%name)
                if (label(2:2) == " ") label(2:2) = "_"
             else if (cp(i)%typ == -3) then
                label = "XX"
             else if (cp(i)%typ == -1)  then
                label = "YY"
             else if (cp(i)%typ == 1) then
                label = "ZZ"
             else
                label = "XZ"
             end if
             write (luout,'(2X,A,3(F10.6," "),A2,I4.4,a)') "neq ",cp(i)%x,label,i," 0"
          end do
       end if
       write (luout,'(X,A)') "endcrystal"
       if (docell) then
          write (luout,'(X,A)') "unitcell radius 0.01 rgb 1.0 0.5 0.5 many"
       end if
       if (doballs) then
          if (nosym) then
             do i = 1,ncpcel
                if (cpcel(i)%idx <= cr%nneq) then
                   label = trim(cr%at(cpcel(i)%idx)%name)
                   if (label(2:2) == " ") label(2:2) = "_"
                else if (cpcel(i)%typ == -3) then
                   label = "XX"
                else if (cpcel(i)%typ == -1)  then
                   label = "YY"
                else if (cpcel(i)%typ == 1) then
                   label = "ZZ"
                else
                   label = "XZ"
                end if
                write (luout,'(X,A5,A,I4.4,A,F6.2)') "ball ",label,i," jmol radius", ball_radius(i,.true.)
             end do
          else
             do i = 1,ncp
                if (i <= cr%nneq) then
                   label = trim(cr%at(i)%name)
                   if (label(2:2) == " ") label(2:2) = "_"
                else if (cp(i)%typ == -3) then
                   label = "XX"
                else if (cp(i)%typ == -1)  then
                   label = "YY"
                else if (cp(i)%typ == 1) then
                   label = "ZZ"
                else
                   label = "XZ"
                end if
                write (luout,'(X,A5,A,I4.4,A,F6.2)') "ball ",label,i," jmol radius", ball_radius(i,.false.)
             end do
          end if
       end if
       ! freehand
       ! 1 - cores | 2 - ncps | 3 - bcps | 4 - rcps | 5 - ccps | 6 - paths | 7 - misc.
       write (luout,'(X,A)') "freehand"
       write (luout,'(A)') "# 1 core, 2 ncp, 3 bcp, 4 rcp, 5 ccp, 6 path, 7 ???"
       write (luout,'(A)') "# 1 dark gray, 2 red, 3 green, 4 light blue, 5 yellow, 6 gold, 7 orange"
       write (luout,'(2X,A)') "type 1 pointrad 0.20 pointrgb 0.2 0.2 0.2"
       write (luout,'(2X,A)') "type 2 pointrad 0.10 pointrgb 0.6274 0.0000 0.2588"
       write (luout,'(2X,A)') "type 3 pointrad 0.06 pointrgb 0.0588 0.5098 0.0588"
       write (luout,'(2X,A)') "type 4 pointrad 0.06 pointrgb 0.7 1.0 1.0"
       write (luout,'(2X,A)') "type 5 pointrad 0.10 pointrgb 1.0000 0.7059 0.2745"
       write (luout,'(2X,A)') "type 6 pointrad 0.02 pointrgb 1.0 0.7 0.3"
       write (luout,'(2X,A)') "type 7 pointrad 0.03 pointrgb 1.0 0.5 0.1"
    elseif (outfmt == "txt") then
       ! connect and initalize flux file
       outfile = trim(fileroot) // "_flux.txt"
       write (uout,'(A,A/)') "* Writing paths to file: ", trim(outfile)
       luout = fopen_write(outfile)
    elseif (outfmt=="obj" .or. outfmt=="ply" .or. outfmt=="off") then
       ! connect and initialize obj/ply/off file 
       outfile = trim(fileroot) // "_flux." // outfmt
       write (uout,'(A,A/)') "* Writing paths to file: ", trim(outfile)
       call struct_write_3dmodel(cr,outfile,outfmt,(/1,1,1/),.true.,.false.,.false.,&
          .true.,.true.,-1d0,(/0d0,0d0,0d0/),-1d0,(/0d0,0d0,0d0/),luout,lumtl)
    elseif (outfmt=="cml") then
       ! connect and initialize cml file 
       outfile = trim(fileroot) // "_flux." // outfmt
       write (uout,'(A,A/)') "* Writing paths to file: ", trim(outfile)
       call struct_write_mol(cr,outfile,outfmt,(/1,1,1/),.true.,.false.,.false.,&
          .false.,0d0,.false.,1,-1d0,(/0d0,0d0,0d0/),-1d0,(/0d0,0d0,0d0/),luout)
    endif

  end subroutine flx_initialize

  !> Deallocate memory, end the tessel output and close files.
  subroutine flx_end(molmotif)
    use tools_io, only: fclose
    use graphics, only: graphics_close
    use global, only: fileroot

    logical, intent(in) :: molmotif

    if (.not.flx_init) return

    ! Nullify variables
    flx_init = .false.

    ! close it
    if (outfmt == "obj" .or. outfmt == "ply" .or. outfmt == "off") then
       call graphics_close(outfmt,luout,lumtl)
    elseif (outfmt == "cml") then
       write (luout,'(" </atomArray>")')
       write (luout,'("</molecule>")')
       call fclose(luout)
    elseif (outfmt == "tss") then
       ! End tessel output
       write (luout,'(X,A)') "endfreehand"
       ! print formats
       if (molmotif) then
          write (luout,'(X,"molmotif allmaincell")') 
       end if
       write (luout,'(X,"vrml ",A,".wrl")') trim(adjustl(fileroot))
       write (luout,'(X,"# povray ",A,".pov")') trim(adjustl(fileroot))

       ! endmolecule
       write (luout,'(A)') "endmolecule"

       ! end tasks
       write (luout,'("# run povray -d +ft +I",A,".pov +O",A,".tga +W2000 +H2000 +A")') &
          trim(adjustl(fileroot)), trim(adjustl(fileroot))
       write (luout,'("# run convert ",A,".tga -bordercolor white -border 1x1 -trim +repage ",A,".png")') &
          trim(adjustl(fileroot)), trim(adjustl(fileroot))
       write (luout,'("# run rm -f ",A,".tga")') trim(adjustl(fileroot))
       write (luout,'("reset")')
       write (luout,'("end")')
       call fclose(luout)
    elseif (outfmt == "txt") then
       call fclose(luout)
    end if

    ! Deallocate memory
    deallocate(flx_x)
    deallocate(flx_rho)
    deallocate(flx_grad)
    deallocate(flx_h)

  end subroutine flx_end

  !> Print gradient path info to standard output.
  subroutine flx_printpath(rgb0)
    use struct_basic, only: cr
    use graphics, only: graphics_ball
    use global, only: dunit
    use tools_io, only: string
    use param, only: bohrtoa
    integer, intent(in) :: rgb0(3)

    integer :: i, j
    real*8 :: x(3)
    real*8, parameter :: rrad = 0.15d0

    real*8 :: maux(4,4)

    if (outfmt == "txt") then
       write (luout,'(A,A50)') "# Gradient path "
       write (luout,'(A,I2)') "# Iup  : ", flx_iup
       write (luout,'(A,I6)') "# Number of points  : ", flx_n

       if (.not.cr%ismolecule) then
          maux = 0d0
          maux(1:3,1:3) = cr%crys2car
          maux(4,4) = 1d0
          write (luout,'(A/1p,4("#",4E22.14/))') "# Crys2Car : ", ((maux(i,j),j=1,4),i=1,4)
          maux(1:3,1:3) = cr%car2crys
          write (luout,'(A/1p,4("#",4E22.14/))') "# Car2Crys : ", ((maux(i,j),j=1,4),i=1,4)
          write (luout,'(A,X,A10,X,12(A20,X))') "#","x","y","z","rho","rhox","rhoy","rhoz",&
             "rhoxx","rhoxy","rhoxz","rhoyy","rhoyz","rhozz"
          do i = 1,flx_n
             write (luout,'(13(E20.12,X))') flx_x(:,i), flx_rho(i),&
                flx_grad(:,i), flx_h(1,:,i), flx_h(2,:,i), flx_h(3,:,i)
          end do
       else
          write (luout,'(A,X,A10,X,12(A20,X))') "#","x","y","z","rho","rhox","rhoy","rhoz",&
             "rhoxx","rhoxy","rhoxz","rhoyy","rhoyz","rhozz"
          do i = 1,flx_n
             x = flx_x(:,i)
             if (cr%ismolecule) &
                x = (cr%x2c(x) + cr%molx0) * dunit
             write (luout,'(13(E20.12,X))') x, flx_rho(i),&
                flx_grad(:,i), flx_h(1,:,i), flx_h(2,:,i), flx_h(3,:,i)
          end do
       end if
       write (luout,'(A/)') "# End gradient path"

    elseif (outfmt == "tss") then
       write (luout,'("# ")') 
       write (luout,'(2X,A)') "curve balls type 6"
       do i=1,flx_n
          write (luout,'(3X,3(E20.12,X))') flx_x(:,i)
       end do
       write (luout,'(2X,A)') "endcurve"
    elseif (outfmt=="obj" .or. outfmt=="off" .or. outfmt=="ply" .or. outfmt=="cml") then
       do i=1,flx_n
          x = cr%x2c(flx_x(:,i))
          if (outfmt == "obj" .or. outfmt == "off" .or. outfmt == "ply") then
             call graphics_ball(outfmt,luout,x,rgb0,rrad)
          elseif (outfmt == "cml") then
             if (.not.cr%ismolecule) then
                write (luout,'("<atom id=""a",A,""" elementType=""Xz"" xFract=""",A,""" yFract=""",A,""" zFract=""",A,"""/>")')&
                   string(i), (trim(string(flx_x(j,i),'f',18,10)),j=1,3)
             else
                x = x + cr%molx0
                write (luout,'("<atom id=""a",A,""" elementType=""Xz"" x3=""",A,""" y3=""",A,""" z3=""",A,"""/>")') &
                   string(i), (trim(string(x(j) * bohrtoa,'f',18,10)),j=1,3)
             end if
          end if
       end do
    end if

  end subroutine flx_printpath

  !> Use flx_printpath to print the path info. Also, obtain the
  !> symmetry equivalent positions of x (cryst. coords.), transform
  !> the path in flx_x with the symmetry operations and print the
  !> results. The behaviour depends on the value of flxsym. If -1,
  !> print only the path without symmetry replication. If 0, print
  !> every copy in the main unit cell.  If > 0, print every copy
  !> contained in flxsym shells of unit cells.  flx_symprintpath is a
  !> wrapper for flx_printpath() so that this routine is not used
  !> directly.
  subroutine flx_symprintpath(x,flxsym,rgb)
    use struct_basic, only: cr

    integer, intent(in) :: flxsym
    real*8, dimension(3), intent(in) :: x
    integer, intent(in) :: rgb(3)

    real*8, parameter :: flx_epsx = 1d-4

    real*8, dimension(3) :: tempx
    real*8, allocatable :: sympos(:,:)
    integer, allocatable :: symrotm(:), symcenv(:)
    integer :: mmult
    integer :: shedge, shsize

    integer :: j, k, l, ln
    real*8, dimension(3) :: templvec

    real*8, parameter :: eps = 1d-5

    ! print symmetry equivalent?
    if (flxsym >= 0) then
       shedge = 2*flxsym + 3
       shsize = shedge**3
       ! reproduce the point using the space group symmetry operations
       tempx = x
       call cr%symeqv(tempx,mmult,sympos,symrotm,symcenv,eps)
       do j = 1,mmult
          ! apply symmetry operations to the path
          do l = 0, shsize-1
             ln = l
             templvec(1) = mod(ln,shedge) - 1
             ln = ln / shedge
             templvec(2) = mod(ln,shedge) - 1
             ln = ln / shedge
             templvec(3) = mod(ln,shedge) - 1

             ! only bcps in the aggregation of cells or its boundary
             if (any(sympos(:,j)+templvec < -flxsym-flx_epsx) .or. any(sympos(:,j)+templvec > 1.d0+flxsym+flx_epsx)) then
                cycle
             end if

             do k = 1, flx_n
                flx_x(:,flx_n+k) = flx_x(:,k)
                flx_x(:,k) = matmul(cr%rotm(:,1:3,symrotm(j)),flx_x(:,k))
                flx_x(:,k) = flx_x(:,k) + cr%rotm(:,4,symrotm(j))
                flx_x(:,k) = flx_x(:,k) + cr%cen(:,symcenv(j))
                flx_x(:,k) = flx_x(:,k) + templvec

                ! return the point to the big cell
                where (flx_x(:,k) > 1.d0+flxsym+flx_epsx) flx_x(:,k) = flx_x(:,k) - 2*flxsym - 1.d0
                where (flx_x(:,k) < -flxsym-flx_epsx)     flx_x(:,k) = flx_x(:,k) + 2*flxsym + 1.d0

             end do

             ! print the new bond path
             call flx_printpath(rgb)

             ! retrieve the original bond path
             flx_x(:,1:flx_n) = flx_x(:,flx_n+1:2*flx_n)
          end do
       end do
       deallocate(sympos,symrotm,symcenv)
    else
       call flx_printpath(rgb)
    end if

  end subroutine flx_symprintpath

  !> Print gradient paths from point x (cryst. coords.), with step
  !> step, gradient norm criterion eps and in the iup
  !> direction. flxsym = -1, ignore symmetry ; 0, apply x symmetry
  !> to the unit cell and its boundary; >0 to a shell of flxsym unit
  !> cells.
  subroutine flx_point(x,iup,flxsym,rgb)
    use navigation, only: gradient, prunepath
    use fields, only: f
    use struct_basic, only: cr
    use global, only: refden, prunedist

    integer, intent(in) :: iup
    real*8, dimension(3), intent(in) :: x
    integer, intent(in) :: flxsym
    integer, intent(in) :: rgb(3)

    real*8, dimension(3) :: tempx, xini
    integer :: ier

    xini = x
    tempx = cr%x2c(x)

    if (iup /= 0) then
       call gradient(f(refden),tempx,iup,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
       call prunepath(cr,flx_n,flx_x,prunedist)
       call flx_symprintpath(xini,flxsym,rgb)
    else
       call gradient(f(refden),tempx,1,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
       call prunepath(cr,flx_n,flx_x,prunedist)
       call flx_symprintpath(xini,flxsym,rgb)
       tempx = cr%x2c(x)
       call gradient(f(refden),tempx,-1,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
       call prunepath(cr,flx_n,flx_x,prunedist)
       call flx_symprintpath(xini,flxsym,rgb)
    end if

  end subroutine flx_point

  !> Print gradient paths starting at the ncp or ccp id (all the ids
  !> are from the complete cp list), with ntheta and nphi points. Add
  !> lvec to the position of the cp. flxsym = -1, ignore symmetry ; 0,
  !> apply x symmetry to the unit cell and its boundary; >0 to a shell
  !> of flxsym unit cells.
  subroutine flx_ncp(id,ntheta,nphi,flxsym,lvec,rgb)
    use navigation, only: gradient, prunepath
    use fields, only: f
    use varbas, only: ncpcel, cpcel
    use global, only: refden, prunedist
    use struct_basic, only: cr
    use tools_io, only: ferror, faterr
    use param, only: pi

    integer, intent(in) :: id, ntheta, nphi
    integer, intent(in) :: flxsym
    integer, intent(in), dimension(3), optional :: lvec
    integer, intent(in) :: rgb(3)

    real*8, parameter :: change = 0.1d0

    real*8, dimension(3) :: xncp, xpoint
    real*8 :: phi, theta
    integer :: i, j
    integer :: ier, iup
    real*8, dimension(3) :: xini

    if (id <= 0 .or. id > ncpcel) call ferror('flx_ncp','CP identifier < 0 or > # of CPs',faterr)

    if (cpcel(id)%typ == -3) then
       iup = -1
    else if (cpcel(id)%typ == 3) then
       iup = +1
    else
       call ferror('flx_ncp','CP is not CCP or NCP',faterr)
    end if
    xncp = cpcel(id)%x
    if (present(lvec)) then
       xncp = xncp + real(lvec,8)
    end if
    xini = xncp
    xncp = cr%x2c(xncp)

    do i = 1,nphi
       do j = 1,ntheta
          phi = pi * real(i,8) / real(nphi+1,8)
          theta = 2.d0 * pi * real(j,8) / real(ntheta,8)
          xpoint = xncp + change * (/ cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi) /)
          call gradient(f(refden),xpoint,iup,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
          call prunepath(cr,flx_n,flx_x,prunedist)
          call flx_symprintpath(xini,flxsym,rgb)
       end do
    end do

    ! special points on the unit sphere
    phi = 0d0
    theta = 0d0
    xpoint = xncp + change * (/ cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi) /)
    call gradient(f(refden),xpoint,iup,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
    call prunepath(cr,flx_n,flx_x,prunedist)
    call flx_symprintpath(xini,flxsym,rgb)
    phi = pi
    theta = 0d0
    xpoint = xncp + change * (/ cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi) /)
    call gradient(f(refden),xpoint,iup,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
    call prunepath(cr,flx_n,flx_x,prunedist)
    call flx_symprintpath(xini,flxsym,rgb)

  end subroutine flx_ncp

  !> Print gradient paths starting at the bcp id (all the ids are from
  !> the complete cp list), in the direction iup. iup = 0, both
  !> direction. npoints only applies to the -1 direction. Add lvec to
  !> the position of the cp. flxsym = -1, ignore symmetry ; 0, apply x
  !> symmetry to the unit cell and its boundary; >0 to a shell of
  !> flxsym unit cells. Bcpmethod determines the method for finding
  !> starting angular grid around the bcp. "bra" (braindead),
  !> uniform distribution; "quo" (quotient), remap the uniform
  !> distribution using the quotient of eigenvalues, "dyn"
  !> (dynamical), use the linearized dynamical system to predict a
  !> good approximation to a grid. "h1 " (heuristic 1), experimental
  !> heuristic method, using information of the partially fluxed
  !> initial points. (h1 and dyn are experimental)
  subroutine flx_bcp(id,iup,npoints,flxsym,bcpmethod,lvec,rgb)
    use navigation, only: gradient, prunepath
    use fields, only: f, grd
    use varbas, only: ncpcel, cpcel
    use struct_basic, only: cr
    use global, only: refden, prunedist
    use tools_math, only: eig
    use tools_io, only: ferror, faterr
    use param, only: pi
    use types, only: scalar_value

    integer, intent(in) :: id, iup, npoints
    integer, intent(in) :: flxsym
    character*3, intent(in) :: bcpmethod
    integer, intent(in), dimension(3), optional :: lvec
    integer, intent(in) :: rgb(3)

    real*8, parameter :: change = 0.1d0

    real*8, dimension (3) :: xbcp, xpoint, xini
    real*8 :: angle, cangle, sangle
    integer :: i
    real*8, dimension(3) :: reval
    real*8, dimension(3,3) :: evec
    real*8, dimension(3) :: v1, v2, vup, v
    real*8 :: ev1, ev2
    integer :: n
    real*8 :: R
    real*8, dimension(npoints) :: thetavec
    integer :: ier, ircp
    type(scalar_value) :: res
    
    if (id <= 0 .or. id > ncpcel) call ferror('flx_bcp','CP identifier < 0 or > # of CPs',faterr)

    if (cpcel(id)%typ == -1) then
       ircp = 1
    else if (cpcel(id)%typ == 1) then
       ircp = -1
    else
       call ferror('flx_bcp','CP is not BCP or RCP',faterr)
    end if

    xbcp = cpcel(id)%x
    if (present(lvec)) then
       xbcp = xbcp + real(lvec,8)
    end if
    xini = xbcp
    xbcp = cr%x2c(xbcp)
    call grd(f(refden),xbcp,2,res)
    evec = res%hf
    call eig(evec,reval)

    if (ircp * reval(1) > 0) then
       vup = evec(:,1)
       if (ircp * reval(2) < ircp * reval(3)) then
          v1 = evec(:,2)
          ev1 = reval(2)
          v2 = evec(:,3)
          ev2 = reval(3)
       else
          v1 = evec(:,3)
          ev1 = reval(3)
          v2 = evec(:,2)
          ev2 = reval(2)
       end if
    else if (ircp * reval(2) > 0) then
       vup = evec(:,2)
       if (ircp * reval(1) < ircp * reval(3)) then
          v1 = evec(:,1)
          ev1 = reval(1)
          v2 = evec(:,3)
          ev2 = reval(3)
       else
          v1 = evec(:,3)
          ev1 = reval(3)
          v2 = evec(:,1)
          ev2 = reval(1)
       end if
    else
       vup = evec(:,3)
       if (ircp * reval(1) < ircp * reval(2)) then
          v1 = evec(:,1)
          ev1 = reval(1)
          v2 = evec(:,2)
          ev2 = reval(2)
       else
          v1 = evec(:,2)
          ev1 = reval(2)
          v2 = evec(:,1)
          ev2 = reval(1)
       end if
    end if

    v1 = v1 / sqrt(dot_product(v1,v1))
    v2 = v2 / sqrt(dot_product(v2,v2))
    vup = vup / sqrt(dot_product(vup,vup))

    if (iup == 0 .or. iup == ircp) then
       xpoint = xbcp + change * vup
       call gradient(f(refden),xpoint,ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
       call prunepath(cr,flx_n,flx_x,prunedist)
       call flx_symprintpath(xini,flxsym,rgb)
       xpoint = xbcp - change * vup
       call gradient(f(refden),xpoint,ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
       call prunepath(cr,flx_n,flx_x,prunedist)
       call flx_symprintpath(xini,flxsym,rgb)
    end if

    if (iup == 0 .or. iup == -ircp) then
       if (bcpmethod == "bra") then
          ! braindead: use an uniform distribution
          ! write (uout,'("+ bra : fluxing the uniform distribution.")')
          do i = 1,npoints
             angle = 2d0 * pi * real(i-1,8)/real(npoints,8)
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * sangle + v2 * cangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)
          end do
       else if (bcpmethod == "quo") then
          ! use eigenvalue quotient to remap the uniform point distribution
          ! write (uout,'("+ quo : fluxing a uniform distribution remapped with the eigen-quotient.")')
          n = npoints / 2
          do i = 1,n
             angle = pi * (real(i,8)-1.d0-(real(n,8)-1.d0)*.5d0)/(real(n,8)-1.d0)
             angle = sign(abs(angle)**(ev1/ev2) / (pi/2.d0)**(ev1/ev2-1),angle)
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * sangle + v2 * cangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)
             angle = angle + pi
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * sangle + v2 * cangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)
          end do
       else
          ! write (uout,'("+ dyn: coarse exploration to determine bcp-ccp distance")')
          ! coarse exploration to determine the distance from bcp to ccps
          r = 1d20
          do i = 1,50
             angle = 2d0*pi * (i-1) / 49d0
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * cangle + v2 * sangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             ! last point before newton -> not converted to the main cell
             xpoint = cr%x2c(flx_x(:,flx_n-1))
             r = min(r,sqrt((xpoint(1)-xbcp(1))**2+(xpoint(2)-xbcp(2))**2+(xpoint(3)-xbcp(3))**2))
          end do

          ! write (uout,'("+ dyn: bcp-ccp minimum distance: ",F20.10)') r

          ! print bcp-ias using an exponentially adapted grid
          n = max(npoints/4-1,2)
          call flx_findthetagrid(ev1,ev2,change,r,n,thetavec)

          ! write (uout,'("+ dyn: fluxing the remapped angles")')
          do i = 1,n
             angle = thetavec(i)
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * cangle + v2 * sangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)

             angle = angle + pi
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * cangle + v2 * sangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)

             angle = -thetavec(i) + pi
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * cangle + v2 * sangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)

             angle = angle + pi
             cangle = cos(angle)
             sangle = sin(angle)
             v = v1 * cangle + v2 * sangle
             xpoint = xbcp + change * v
             call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
             call prunepath(cr,flx_n,flx_x,prunedist)
             call flx_symprintpath(xini,flxsym,rgb)

          end do

          ! Manually do 0, pi/2, pi, 3pi/2
          ! write (uout,'("+ dyn: fluxing special angles.")')
          xpoint = xbcp + change * v1
          call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
          call prunepath(cr,flx_n,flx_x,prunedist)
          call flx_symprintpath(xini,flxsym,rgb)
          xpoint = xbcp + change * v2
          call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
          call prunepath(cr,flx_n,flx_x,prunedist)
          call flx_symprintpath(xini,flxsym,rgb)
          xpoint = xbcp - change * v1
          call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
          call prunepath(cr,flx_n,flx_x,prunedist)
          call flx_symprintpath(xini,flxsym,rgb)
          xpoint = xbcp - change * v2
          call gradient(f(refden),xpoint,-ircp,flx_n,flx_m,ier,2,flx_x,flx_rho,flx_grad,flx_h,up2beta=.false.)
          call prunepath(cr,flx_n,flx_x,prunedist)
          call flx_symprintpath(xini,flxsym,rgb)
       end if
       ! write (uout,*)
    end if

  end subroutine flx_bcp

  !> Represents the graph. flxsym is the number of cell shells
  !> that are going to be represented. With flxsym = -1, just the
  !> indicated cp is represented.  flxsym = 0, fill the first
  !> cell. flxsym > 0, fill flxsym shells of unit cells.  igraph is
  !> an integer containing 2 flags: 1, represent ring paths; 2,
  !> represent bond paths. The value of igraph is the sum of
  !> the selected integers. cpid and lvec are optional
  !> parameters. If cpid is not given, the whole graph in the first
  !> unit cell is drawn, with max(0,flxsym) additional shells. If
  !> cpid is present, the lines and paths associated to that cp are
  !> represented, with a lattice displacement given by the optional
  !> integer lvec vector.
  subroutine flx_graph(flxsym,igraph,cpid,lvec,rgb)
    use varbas, only: ncp, cp, cpcel, ncpcel

    integer, intent(in) :: flxsym
    integer, intent(in) :: igraph
    integer, intent(in), optional :: cpid
    integer, intent(in), dimension(3), optional :: lvec
    integer, intent(in) :: rgb(3)

    integer :: i, l, ln
    integer :: m, temp
    logical :: dorcp, dobcp
    integer, dimension(3) :: templvec, locallvec
    integer :: localsym

    ! parse flags
    dorcp = .false.
    dobcp = .false.
    temp = igraph
    if (mod(temp,2) == 1) dorcp = .true.
    temp = temp / 2
    if (mod(temp,2) == 1) dobcp = .true.
    temp = temp / 2

    ! lattice vector
    if (present(lvec)) then
       locallvec = lvec
    else
       locallvec = (/ 0, 0, 0/)
    end if

    ! graphcp case
    if (present(cpid)) then
       ! ncp
       if (cpcel(cpid)%typ == -3) then
          ! all the bond paths of the bcp bonded to the ncp
          do i = 1, ncpcel
             if (cpcel(i)%typ /= -1) cycle
             do l = 0, 26
                ln = l
                templvec(1) = mod(ln,3) - 1
                ln = ln / 3
                templvec(2) = mod(ln,3) - 1
                ln = ln / 3
                templvec(3) = mod(ln,3) - 1
                if (cpcel(i)%ipath(1) == cpid .and. all(locallvec == cpcel(i)%ilvec(:,1) + templvec) .or. &
                   cpcel(i)%ipath(2) == cpid .and. all(locallvec == cpcel(i)%ilvec(:,2) + templvec)) then
                   call flx_bcp(i,1,1,flxsym,"dyn",templvec,rgb)
                end if
             end do
          end do

       ! ccp
       else if (cpcel(cpid)%typ == 3) then
          ! all the ring paths of the rcp bonded to the ccp
          do i = 1, ncpcel
             if (cpcel(i)%typ /= 1) cycle
             do l = 0, 26
                ln = l
                templvec(1) = mod(ln,3) - 1
                ln = ln / 3
                templvec(2) = mod(ln,3) - 1
                ln = ln / 3
                templvec(3) = mod(ln,3) - 1
                if (cpcel(i)%ipath(1) == cpid .and. all(locallvec == cpcel(i)%ilvec(:,1) + templvec) .or. &
                   cpcel(i)%ipath(2) == cpid .and. all(locallvec == cpcel(i)%ilvec(:,2) + templvec)) then
                   call flx_bcp(i,-1,1,flxsym,"dyn",templvec,rgb)
                end if
             end do
          end do
          
       ! bcp
       ! bond paths for the bcp
       else if (cpcel(cpid)%typ == -1) then
          call flx_bcp(cpid,1,1,flxsym,"dyn",locallvec,rgb)

       ! rcp
       ! ring paths for the rcp
       else if (cpcel(cpid)%typ == 1) then
          call flx_bcp(cpid,-1,1,flxsym,"dyn",locallvec,rgb)
       end if

    else
       ! graph case
       localsym = max(flxsym,0)
       m = 0
       do i = 1, ncp
          if (dobcp) then
             ! compute bcp bond paths
             if (cp(i)%typ == -1) then
                call flx_bcp(m+1,1,1,localsym,"dyn",rgb=rgb)
             end if
          end if

          ! compute rcp ring paths
          if (dorcp) then
             if (cp(i)%typ == 1) then
                call flx_bcp(m+1,-1,1,localsym,"dyn",rgb=rgb)
             end if
          end if

          m = m + cp(i)%mult
       end do
    end if

  end subroutine flx_graph

  !> xxxx experimental xxxx
  !> Find a vector thetavec of angles (radian) such that an uniform
  !> angular distribution of points on the boundary circle of radius
  !> R around a bcp with (negative) eigenvalues lx, ly transform to
  !> the points (r0,thetavec) using the electron density gradient
  !> flux. Always R > r0. This routine is used to determine a set of
  !> starting angles that samples uniformly a given bcp-IAS.
  subroutine flx_findthetagrid(lx,ly,r0,R,n,thetavec)
    use param, only: pi

    real*8, intent(in) :: lx, ly, r0, R
    integer, intent(in) :: n
    real*8, dimension(n), intent(out) :: thetavec

    integer :: i
    real*8 :: tmin, tmax, ta, tb, told, localf, localfp
    real*8, dimension(n) :: t
    logical :: critical

    tmin = (log(r0/R)/min(lx,ly))
    tmax = (log(r0/R)/max(lx,ly))

    ! write (uout,*) "+ dyn: findthetagrid quantities."
    do i = 1, n
       ! Fill the initial theta vector
       thetavec(i) = pi / 2d0 * real(i,8) / real(n+1,8)

       ! Numerically find the integration time in the linearized system (Newton)
       t(i) = (tmin+tmax) / 2
       localf = cos(thetavec(i))**2 * exp(2d0*lx*t(i)) + sin(thetavec(i))**2 * exp(2d0*ly*t(i)) - (r0/R)**2
       critical = .false.
       do while (abs(localf) > 1.d-10)
          localf = cos(thetavec(i))**2 * exp(2d0*lx*t(i)) + sin(thetavec(i))**2 * exp(2d0*ly*t(i)) - (r0/R)**2
          localfp = cos(thetavec(i))**2 * exp(2d0*lx*t(i)) *2d0*lx + sin(thetavec(i))**2 * exp(2d0*ly*t(i)) *2d0*ly
          told = t(i)
          t(i) = t(i) - localf / localfp
          if (t(i) > tmax .or. t(i) < tmin) then
             critical = .true.
             exit
          end if
       end do

       ! If newton did not converge, fall back to bisection
       if (critical) then
          ta = tmin
          tb = tmax
          t(i) = (ta+tb)/2d0
          do while (abs(ta-tb) > 0.5d-10)
             localf = cos(thetavec(i))**2 * exp(2d0*lx*t(i)) + sin(thetavec(i))**2 * exp(2d0*ly*t(i)) - (r0/R)**2
             if (localf < 0) then
                tb = t(i)
             else
                ta = t(i)
             end if
             t(i) = (ta+tb)/2d0
          end do
       end if

       thetavec(i) = atan2(sin(thetavec(i))*exp((ly-lx)*t(i)),cos(thetavec(i)))

    end do

  end subroutine flx_findthetagrid

  !> Given the complete or reduced list index of a CP, 
  !> output the appropriate ball radius.
  function ball_radius(i,cel)
    use varbas, only: cpcel, cp
    use struct_basic, only: cr
    use param, only: maxzat
    
    integer, intent(in) :: i
    logical, intent(in) :: cel
    real*8 :: ball_radius

    integer :: typ, z
    
    z = 0
    if (cel) then
       if (cpcel(i)%isnuc) then
          typ = 0
          z = cr%at(cr%atcel(i)%idx)%z
       else
          typ = cpcel(i)%typ
       end if
    else
       if (cp(i)%isnuc) then
          typ = 0
          z = cr%at(i)%z
       else
          typ = cp(i)%typ
       end if
    end if

    ball_radius = 0.5d0
    select case(typ)
    case(3)
       ball_radius = 0.05d0
    case(1)
       ball_radius = 0.05d0
    case(-1)
       ball_radius = 0.1d0
    case(-3)
       ball_radius = 0.2d0
    case(0)
       if (z <= 2) then ! H-He
          ball_radius = 0.25d0
       else if (z <= 10) then ! Li-Ne
          ball_radius = 0.35d0
       else if (z <= 18) then ! Na-Ar
          ball_radius = 0.45d0
       else if (z <= 36) then ! K-Kr
          ball_radius = 0.55d0
       else if (z <= 54) then ! Rb-Xe
          ball_radius = 0.60d0
       else if (z <= 86) then ! Cs-Rn
          ball_radius = 0.650
       else if (z <= maxzat) then ! Fr-...
          ball_radius = 0.70d0
       else
          ball_radius = 0.20d0 ! cps
       end if
    end select
  end function ball_radius

end module flux
