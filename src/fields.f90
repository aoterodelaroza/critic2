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

! Evaluation of scalar fields
module fields
  use types
  implicit none

  private

  integer, parameter, public :: type_promol = 0 !< promolecular density
  integer, parameter, public :: type_grid = 1 !< grid format
  integer, parameter, public :: type_wien = 2 !< wien2k format
  integer, parameter, public :: type_elk  = 3 !< elk format
  integer, parameter, public :: type_pi   = 4 !< pi format
  integer, parameter, public :: type_wfn  = 6 !< molecular wavefunction format
  integer, parameter, public :: type_promol_frag = 7 !< promolecular density from a fragment
  integer, parameter, public :: type_ghost = 8 !< a ghost field

  public :: fields_load
  public :: fields_unload
  public :: fields_init
  public :: fields_end
  public :: fields_propty
  public :: fields_integrable
  public :: fields_integrable_report
  public :: fields_pointprop
  public :: fields_pointprop_report
  public :: goodfield
  public :: getfieldnum
  public :: setfield
  public :: grd
  public :: grd0
  public :: grdall
  public :: fields_typestring
  public :: benchmark
  public :: taufromelf
  public :: testrmt
  public :: fields_fcheck
  public :: fields_feval

  ! integrable properties
  integer, parameter, public :: mprops = 101
  integer, public :: nprops
  type(integrable), public :: integ_prop(mprops)
  integer, parameter, public :: itype_v = 1
  integer, parameter, public :: itype_f = 2
  integer, parameter, public :: itype_fval = 3
  integer, parameter, public :: itype_gmod = 4
  integer, parameter, public :: itype_lap = 5
  integer, parameter, public :: itype_lapval = 6
  integer, parameter, public :: itype_expr = 7
  integer, parameter, public :: itype_mpoles = 8
  integer, parameter, public :: itype_deloc = 9
  character*10, parameter, public :: itype_names(9) = (/&
     "Volume    ","Field     ","Field (v) ","Gradnt mod","Laplacian ",&
     "Laplcn (v)","Expression","Multipoles","Deloc indx"/)
  
  ! pointprop properties
  integer, public :: nptprops
  type(pointpropable), public :: point_prop(mprops)

  ! scalar fields
  integer, public :: mf
  type(field), allocatable, public :: f(:)
  logical, allocatable, public :: fused(:)

  ! eps to move to the main cell
  real*8, parameter :: flooreps = 1d-4 ! border around unit cell

contains

  subroutine fields_init()

    integer :: i

    ! allocate space for all fields
    mf = 5
    if (.not.allocated(fused)) allocate(fused(0:mf))
    fused = .false.
    if (.not.allocated(f)) allocate(f(0:mf))
    do i = 0, mf
       f(i)%init = .false.
       f(i)%numerical = .false.
    end do

    ! initialize the promolecular density
    f(0)%init = .true.
    f(0)%type = type_promol
    f(0)%usecore = .false.
    f(0)%typnuc = -3
    fused(0) = .true.

    ! integrable properties, initialize
    nprops = 1
    integ_prop(:)%used = .false.

    ! point properties, initialize
    nptprops = 0

    ! integrable: volume
    integ_prop(1)%used = .true.
    integ_prop(1)%itype = itype_v
    integ_prop(1)%fid = 0
    integ_prop(1)%prop_name = "Volume"

  end subroutine fields_init

  subroutine fields_end()

    if (allocated(fused)) deallocate(fused)
    if (allocated(f)) deallocate(f)

  end subroutine fields_end

  subroutine fields_load(line,id)
    use elk_private
    use wien_private
    use wfn_private
    use pi_private
    use struct_basic
    use grid_tools
    use global
    use arithmetic
    use tools_io
    use types
    use param

    character*(*), intent(in) :: line
    integer, intent(out) :: id

    integer :: lp, lp2, i, j, id1, id2
    character(len=:), allocatable :: file, wext1, wext2, word, word2, file2, file3, expr
    integer :: zz, n(3)
    logical :: ok
    real*8 :: renv0(3,cr%nenv), xp(3), rhopt
    integer :: idx0(cr%nenv), zenv0(cr%nenv), ix, iy, iz, oid, oid2
    real*8 :: xd(3,3)
    integer :: nid, nwan
    integer, allocatable :: idlist(:)
    type(fragment) :: fr
    logical :: isfrag

    ! read and parse
    lp=1
    file = getword(line,lp)
    word = file(index(file,dirsep,.true.)+1:)
    wext1 = word(index(word,'.',.true.)+1:)
    word = file(index(file,dirsep,.true.)+1:)
    wext2 = word(index(word,'_',.true.)+1:)

    ! if it is one of the special keywords, change the extension
    ! and read the next
    if (equal(file,"wien")) then
       file = getword(line,lp)
       wext1 = "clmsum"
       wext2 = wext1
    elseif (equal(file,"elk")) then
       file = getword(line,lp)
       wext1 = "OUT"
       wext2 = wext1
    elseif (equal(file,"pi")) then
       file = getword(line,lp)
       wext1 = "ion"
       wext2 = wext1
    elseif (equal(file,"cube")) then
       file = getword(line,lp)
       wext1 = "cube"
       wext2 = wext1
    elseif (equal(file,"abinit")) then
       file = getword(line,lp)
       wext1 = "DEN"
       wext2 = wext1
    elseif (equal(file,"vasp")) then
       file = getword(line,lp)
       wext1 = "CHGCAR"
       wext2 = wext1
    elseif (equal(file,"vaspchg")) then
       file = getword(line,lp)
       wext1 = "CHG"
       wext2 = wext1
    elseif (equal(file,"qub")) then
       file = getword(line,lp)
       wext1 = "qub"
       wext2 = wext1
    elseif (equal(file,"xsf")) then
       file = getword(line,lp)
       wext1 = "xsf"
       wext2 = wext1
    elseif (equal(file,"elfgrid")) then
       file = getword(line,lp)
       wext1 = "grid"
       wext2 = wext1
    elseif (equal(file,"siesta")) then
       file = getword(line,lp)
       wext1 = "RHO"
       wext2 = wext1
    elseif (equal(file,"wannier")) then
       file = ""
       wext1 = "wannier"
       wext2 = wext1
    elseif (equal(file,"as")) then
       file = ""
       wext1 = "as"
       wext2 = wext1
    elseif (equal(file,"copy")) then
       ok = eval_next(oid,line,lp)
       if (.not.ok) call ferror("fields_load","wrong LOAD COPY syntax",faterr,line)
       word = lgetword(line,lp)
       if (equal(word,'to')) then
          ok = eval_next(id,line,lp)
          if (id < 0 .or. id > mf) &
             call ferror("fields_load","erroneous field id in LOAD COPY",faterr,line)
          if (.not.ok) call ferror("fields_load","wrong LOAD COPY syntax",faterr,line)
          if (fused(id)) call fields_unload(id)
       else
          id = getfieldnum()
       end if
       f(id) = f(oid)
       fused(id) = .true.
       write (uout,'("* COPIED scalar field from slot ",A," to ",A/)') string(oid), string(id)
       return
    end if

    ! allocate slot
    id = getfieldnum()

    write (uout,'("* LOAD scalar field in slot number: ",A)') string(id)
    if (equal(wext1,'cube')) then
       call grid_read_cube(file,f(id))
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'DEN') .or. equal(wext2,'DEN')) then
       call grid_read_abinit(file,f(id))
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'RHO') .or. equal(wext1,'BADER') .or.&
       equal(wext1,'DRHO') .or. equal(wext1,'LDOS') .or.&
       equal(wext1,'VT') .or. equal(wext1,'VH')) then
       call grid_read_siesta(file,f(id))
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'CHGCAR')) then
       call grid_read_vasp(file,f(id),cr%omega)
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'CHG') .or. equal(wext1,'ELFCAR')) then
       call grid_read_vasp(file,f(id),1d0)
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'qub')) then
       call grid_read_qub(file,f(id))
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'xsf')) then
       call grid_read_xsf(file,f(id))
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'grid')) then
       call grid_read_elk(file,f(id))
       f(id)%type = type_grid
       f(id)%file = file
    else if (equal(wext1,'wfn')) then
       call wfn_read_wfn(file,f(id))
       f(id)%type = type_wfn
       call wfn_register_struct(cr%ncel,cr%atcel)
       f(id)%file = file
    else if (equal(wext1,'wfx')) then
       call wfn_read_wfx(file,f(id))
       f(id)%type = type_wfn
       call wfn_register_struct(cr%ncel,cr%atcel)
       f(id)%file = file
    else if (equal(wext1,'clmsum')) then
       file2 = getword(line,lp)
       call wien_read_clmsum(file,file2,f(id))
       f(id)%type = type_wien
       f(id)%file = file
    else if (equal(wext1,'OUT')) then
       file2 = getword(line,lp)
       file3 = getword(line,lp)
       if (file3 == "") then
          call elk_read_out(f(id),file,file2)
       else
          call elk_read_out(f(id),file,file2,file3)
       end if
       f(id)%type = type_elk
       f(id)%file = file
    else if (equal(wext1,'promolecular')) then
       word = getword(line,lp)
       if (len_trim(word) > 0) then
          fr = cr%identify_fragment_from_xyz(word)
          f(id)%type = type_promol_frag
          f(id)%fr = fr
       else
          f(id)%type = type_promol
       end if
       f(id)%init = .true.
    else if (equal(wext1,'ion')) then
       ! read all the ions
       f(id)%file = ""
       do while(.true.)
          lp2 = lp
          ok = eval_next(zz,line,lp)
          if (ok) then
             call pi_read_ion(file,f(id),zz)
          else
             word = lgetword(line,lp2)
             zz = zatguess(word)
             if (zz == -1) call ferror("fields_load","syntax: load file.ion {atidx/atsym} ...",faterr,line)
             do i = 1, cr%nneq
                if (cr%at(i)%z == zz) call pi_read_ion(file,f(id),i)
             end do
             f(id)%file = adjustl(trim(f(id)%file) // " " // file)
             lp = lp2
          endif
          lp2 = lp
          file = getword(line,lp)
          wext1 = file(index(file,'.',.true.)+1:)
          if (.not.equal(wext1,'ion')) then
             lp = lp2
             exit
          end if
       end do
       f(id)%type = type_pi
       f(id)%init = .true.
       f(id)%exact = .false.
       f(id)%file = trim(f(id)%file)

       ! register structural info in the pi module
       do i = 1, cr%nenv
          renv0(:,i) = cr%atenv(i)%r
          idx0(i) = cr%atenv(i)%idx
          zenv0(i) = cr%at(cr%atenv(i)%idx)%z
       end do
       call pi_register_struct(cr%nenv,renv0,idx0,zenv0)

       ! fill the interpolation tables of the field
       call fillinterpol(f(id))

       ! output info
       write (uout,'("* aiPI input, from ion files")')
       write (uout,'(2X,A,2X,A)') string("Atom"), string("File")
       do i = 1, cr%nneq
          if (.not.f(id)%pi_used(i)) then
             write (uout,'("Non-equivalent atom missing in the load : ",I3)') i
             call ferror("field_load","missing atoms in pi field description",faterr)
          end if
          write (uout,'(4X,A,2X,A)') string(i), string(f(id)%piname(i))
       end do
       write (uout,*)
    else if (equal(wext1,'wannier')) then
       ok = isinteger(n(1),line,lp)
       ok = ok .and. isinteger(n(2),line,lp)
       ok = ok .and. isinteger(n(3),line,lp)
       if (.not.ok) &
          call ferror("fields_load","Error reading wannier supercell",faterr,line)

       nwan = 0
       f(id)%file = ""
       do while(.true.)
          nwan = nwan + 1
          file = getword(line,lp)
          if (len_trim(file) == 0) exit
          ! read another wannier function and accumulate
          call grid_read_xsf(file,f(id),nwan,n,cr%omega)
          f(id)%file = f(id)%file // file // " "
       end do
       f(id)%type = type_grid
       f(id)%file = trim(f(id)%file)

    else if (equal(wext1,'as')) then
       lp2 = lp
       word = getword(line,lp)
       if (equal(word,"promolecular").or.equal(word,"core").or.&
           equal(word,"grad").or.equal(word,"lap").or.equal(word,"taufromelf")) then
          ! load as {promolecular,core,grad,lap,taufromelf}

          if (equal(word,"promolecular").or.equal(word,"core")) then
             ok = eval_next(n(1),line,lp)
             if (ok) then
                ok = ok .and. eval_next(n(2),line,lp)
                ok = ok .and. eval_next(n(3),line,lp)
                if (.not.ok) &
                   call ferror("fields_load","wrong grid size for promolecular/core",faterr,line)
             else
                word2 = getword(line,lp)
                if (equal(word2,"sizeof")) then
                   ok = eval_next(zz,line,lp)
                   if (.not.ok) call ferror("fields_load","wrong sizeof",faterr,line)
                   if (.not.goodfield(zz,type_grid)) call ferror("fields_load","field not allocated",faterr,line)
                   n = f(zz)%n
                else
                   call ferror("fields_load","wrong grid size",faterr,line)
                end if
             end if
          else
             ok = eval_next(oid,line,lp)
             if (.not.ok) &
                call ferror("fields_load","wrong grid for lap/grad/etc.",faterr,line)
             if (.not.goodfield(oid)) call ferror("fields_load","field not allocated",faterr,line)
             n = f(oid)%n
          end if

          if (equal(word,"promolecular").or.equal(word,"core").or.goodfield(oid,type_grid)) then 
             if (equal(word,"promolecular").or.equal(word,"core")) then
                ! maybe we are given a fragment?
                word2 = getword(line,lp)
                isfrag = (len_trim(word2) > 0)
                if (isfrag) fr = cr%identify_fragment_from_xyz(word2)
             end if
             ! build a grid 
             f(id)%init = .true.
             f(id)%n = n
             f(id)%type = type_grid
             f(id)%mode = mode_trispline
             allocate(f(id)%f(n(1),n(2),n(3)))
             if (equal(word,"promolecular")) then
                if (isfrag) then
                   call grid_rhoat(f(id),f(id),2,fr)
                else
                   call grid_rhoat(f(id),f(id),2)
                end if
             elseif (equal(word,"core")) then
                if (isfrag) then
                   call grid_rhoat(f(id),f(id),3,fr)
                else
                   call grid_rhoat(f(id),f(id),3)
                end if
             elseif (equal(word,"lap")) then
                call grid_laplacian(f(oid),f(id))
             elseif (equal(word,"grad")) then
                call grid_gradrho(f(oid),f(id))
             end if
          elseif ((goodfield(oid,type_wien).or.goodfield(oid,type_elk)).and.equal(word,"lap")) then
             ! calculate the laplacian of a lapw scalar field
             f(id) = f(oid)
             if (f(oid)%type == type_wien) then
                call wien_tolap(f(id))
             else
                call elk_tolap(f(id))
             end if
          else
             call ferror("fields_load","Incorrect scalar field type",faterr,line)
          endif

       elseif (equal(word,"clm")) then
          ! load as clm
          word = lgetword(line,lp)
          if (equal(word,'add').or.equal(word,'sub')) then
             ok = eval_next(id1,line,lp)
             ok = ok .and. eval_next(id2,line,lp)
             if (.not.ok) call ferror("fields_load","wrong syntax in LOAD AS CLM",faterr,line)
             if (.not.goodfield(id1).or..not.goodfield(id2)) call ferror("fields_load","field not allocated",faterr,line)
             if (f(id1)%type/=f(id2)%type) call ferror("fields_load","fields not the same type",faterr,line)
             if (f(id1)%type/=type_wien.and.f(id1)%type/=type_elk) call ferror("fields_load","incorrect type",faterr,line) 
             if (f(id1)%type == type_wien) then
                ! wien
                f(id) = f(id1)
                if (any(shape(f(id1)%slm) /= shape(f(id2)%slm))) &
                   call ferror("fields_load","nonconformant slm in wien fields",faterr,line) 
                if (any(shape(f(id1)%sk) /= shape(f(id2)%sk))) &
                   call ferror("fields_load","nonconformant sk in wien fields",faterr,line) 
                if (equal(word,'add')) then
                   f(id)%slm = f(id)%slm + f(id2)%slm
                   f(id)%sk = f(id)%sk + f(id2)%sk
                   if (allocated(f(id)%ski)) &
                      f(id)%ski = f(id)%ski + f(id2)%ski
                else
                   f(id)%slm = f(id)%slm - f(id2)%slm
                   f(id)%sk = f(id)%sk - f(id2)%sk
                   if (allocated(f(id)%ski)) &
                      f(id)%ski = f(id)%ski - f(id2)%ski
                end if
             else
                ! elk
                f(id) = f(id1)
                if (any(shape(f(id1)%rhomt) /= shape(f(id2)%rhomt))) &
                   call ferror("fields_load","nonconformant rhomt in elk fields",faterr,line) 
                if (any(shape(f(id1)%rhok) /= shape(f(id2)%rhok))) &
                   call ferror("fields_load","nonconformant rhok in elk fields",faterr,line) 
                if (equal(word,'add')) then
                   f(id)%rhomt = f(id)%rhomt + f(id2)%rhomt
                   f(id)%rhok = f(id)%rhok + f(id2)%rhok
                else
                   f(id)%rhomt = f(id)%rhomt - f(id2)%rhomt
                   f(id)%rhok = f(id)%rhok - f(id2)%rhok
                end if
             end if
          else
             call ferror("fields_load","invalid CLM syntax",faterr,line)
          end if
       elseif (isexpression_or_word(expr,line,lp2)) then
          ! load as "blehblah"
          lp = lp2
          ok = eval_next(n(1),line,lp)
          if (ok) then
             ! load as "blehblah" n1 n2 n3 
             ok = ok .and. eval_next(n(2),line,lp)
             ok = ok .and. eval_next(n(3),line,lp)
          else
             lp2 = lp
             word2 = getword(line,lp)
             if (equal(word2,"sizeof")) then
                ! load as "blehblah" sizeof
                ok = eval_next(zz,line,lp)
                if (.not.ok) call ferror("fields_load","wrong sizeof",faterr,line)
                if (.not.goodfield(zz,type_grid)) call ferror("fields_load","field not allocated",faterr,line)
                n = f(zz)%n
             elseif (equal(word2,"ghost")) then
                ! load as "blehblah" ghost
                ! ghost field -> no grid size
                n = 0
                call fields_in_eval(expr,nid,idlist)
                if (nid == 0) call ferror("fields_load","No fields referenced in LOAD AS expression",faterr,line)
             else
                ! load as "blehblah"
                lp = lp2
                call fields_in_eval(expr,nid,idlist)
                n = 0
                do i = 1, nid
                   if (.not.goodfield(idlist(i),type_grid)) cycle
                   do j = 1, 3
                      n(j) = max(n(j),f(idlist(i))%n(j))
                   end do
                end do
                deallocate(idlist)
                if (any(n < 1)) then
                   ! nothing read, must be a ghost
                   call fields_in_eval(expr,nid,idlist)
                   if (nid == 0) call ferror("fields_load","No fields referenced in LOAD AS expression",faterr,line)
                endif
             end if
          end if

          if (any(n < 1)) then
             ! Ghost field
             do i = 1, nid
                if (.not.fused(idlist(i))) &
                   call ferror('fields_load','Tried to define ghost field using and undefined field',faterr)
             end do

             f(id)%init = .true.
             f(id)%type = type_ghost
             f(id)%expr = string(expr)
             f(id)%nf = nid
             if (allocated(f(id)%fused)) deallocate(f(id)%fused)
             allocate(f(id)%fused(nid))
             f(id)%fused = idlist
             f(id)%usecore = .false.
             f(id)%numerical = .true.
             write (uout,'("* GHOST field with expression: ",A)') expr
             write (uout,'("  Using NUMERICAL derivatives")') 
          else
             ! Grid field
             ! prepare the cube
             f(id)%init = .true.
             f(id)%n = n
             f(id)%type = type_grid
             f(id)%mode = mode_trispline
             allocate(f(id)%f(n(1),n(2),n(3)))

             ! dimensions
             xd = eye
             do i = 1, 3
                xd(:,i) = cr%x2c(xd(:,i))
                xd(:,i) = xd(:,i) / real(n(i),8)
             end do

             write (uout,'("* GRID from expression: ",A)') expr
             !$omp parallel do private (xp,rhopt) schedule(dynamic)
             do ix = 0, n(1)-1
                do iy = 0, n(2)-1
                   do iz = 0, n(3)-1
                      xp = real(ix,8) * xd(:,1) + real(iy,8) * xd(:,2) &
                         + real(iz,8) * xd(:,3)
                      rhopt = eval_hard_fail(expr,xp,fields_fcheck,fields_feval)
                      !$omp critical (fieldwrite)
                      f(id)%f(ix+1,iy+1,iz+1) = rhopt
                      !$omp end critical (fieldwrite)
                   end do
                end do
             end do
             !$omp end parallel do
             write (uout,'("  Grid dimensions : ",3(A,2X))') (string(f(id)%n(j)),j=1,3)
             write (uout,'("  First elements... ",3(A,2X))') (string(f(id)%f(1,1,j),'e',decimal=12),j=1,3)
             write (uout,'("  Last elements... ",3(A,2X))') (string(f(id)%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=0,2)
             write (uout,'("  Sum of elements... ",A)') string(sum(f(id)%f(:,:,:)),'e',decimal=12)
             write (uout,'("  Sum of squares of elements... ",A)') string(sum(f(id)%f(:,:,:)**2),'e',decimal=12)
             write (uout,*)
          end if
       else
          call ferror("fields_load","invalid AS syntax",faterr,line)
       end if
    else
       call ferror('fields_load','unrecognized file format',faterr,line)
    end if

    ! fill misc values
    f(id)%c2x = cr%car2crys
    f(id)%x2c = cr%crys2car
    f(id)%typnuc = -3

    ! use cores?
    if (f(id)%type == type_grid) then
       f(id)%usecore = .true.
    else
       f(id)%usecore = .false.
    end if

    ! parse the rest of the line
    call setfield(id,line(lp:))

    ! test the muffin tin discontinuity
    call testrmt(id,0)

  end subroutine fields_load

  subroutine fields_unload(id)

    integer, intent(in) :: id

    fused(id) = .false.
    f(id)%init = .false.
    f(id)%name = ""
    f(id)%file = ""
    if (allocated(f(id)%f)) deallocate(f(id)%f)
    if (allocated(f(id)%fwan)) deallocate(f(id)%fwan)
    if (allocated(f(id)%c2)) deallocate(f(id)%c2)
    if (allocated(f(id)%lm)) deallocate(f(id)%lm)
    if (allocated(f(id)%lmmax)) deallocate(f(id)%lmmax)
    if (allocated(f(id)%slm)) deallocate(f(id)%slm)
    if (allocated(f(id)%sk)) deallocate(f(id)%sk)
    if (allocated(f(id)%ski)) deallocate(f(id)%ski)
    if (allocated(f(id)%tauk)) deallocate(f(id)%tauk)
    if (allocated(f(id)%tauki)) deallocate(f(id)%tauki)
    if (allocated(f(id)%rotloc)) deallocate(f(id)%rotloc)
    if (allocated(f(id)%rnot)) deallocate(f(id)%rnot)
    if (allocated(f(id)%rmt)) deallocate(f(id)%rmt)
    if (allocated(f(id)%dx)) deallocate(f(id)%dx)
    if (allocated(f(id)%jri)) deallocate(f(id)%jri)
    if (allocated(f(id)%multw)) deallocate(f(id)%multw)
    if (allocated(f(id)%iatnr)) deallocate(f(id)%iatnr)
    if (allocated(f(id)%pos)) deallocate(f(id)%pos)
    if (allocated(f(id)%iop)) deallocate(f(id)%iop)
    if (allocated(f(id)%iz)) deallocate(f(id)%iz)
    if (allocated(f(id)%tau)) deallocate(f(id)%tau)
    if (allocated(f(id)%krec)) deallocate(f(id)%krec)
    if (allocated(f(id)%inst)) deallocate(f(id)%inst)
    if (allocated(f(id)%spr)) deallocate(f(id)%spr)
    if (allocated(f(id)%spr_a)) deallocate(f(id)%spr_a)
    if (allocated(f(id)%spr_b)) deallocate(f(id)%spr_b)
    if (allocated(f(id)%nrmt)) deallocate(f(id)%nrmt)
    if (allocated(f(id)%vgc)) deallocate(f(id)%vgc)
    if (allocated(f(id)%igfft)) deallocate(f(id)%igfft)
    if (allocated(f(id)%xcel)) deallocate(f(id)%xcel)
    if (allocated(f(id)%iesp)) deallocate(f(id)%iesp)
    if (allocated(f(id)%rhomt)) deallocate(f(id)%rhomt)
    if (allocated(f(id)%rhok)) deallocate(f(id)%rhok)
    if (allocated(f(id)%piname)) deallocate(f(id)%piname)
    if (allocated(f(id)%pi_used)) deallocate(f(id)%pi_used)
    if (allocated(f(id)%nsym)) deallocate(f(id)%nsym)
    if (allocated(f(id)%naos)) deallocate(f(id)%naos)
    if (allocated(f(id)%naaos)) deallocate(f(id)%naaos)
    if (allocated(f(id)%nsto)) deallocate(f(id)%nsto)
    if (allocated(f(id)%nasto)) deallocate(f(id)%nasto)
    if (allocated(f(id)%nn)) deallocate(f(id)%nn)
    if (allocated(f(id)%z)) deallocate(f(id)%z)
    if (allocated(f(id)%xnsto)) deallocate(f(id)%xnsto)
    if (allocated(f(id)%c)) deallocate(f(id)%c)
    if (allocated(f(id)%nelec)) deallocate(f(id)%nelec)
    if (allocated(f(id)%pgrid)) deallocate(f(id)%pgrid)
    if (allocated(f(id)%icenter)) deallocate(f(id)%icenter)
    if (allocated(f(id)%itype)) deallocate(f(id)%itype)
    if (allocated(f(id)%e)) deallocate(f(id)%e)
    if (allocated(f(id)%occ)) deallocate(f(id)%occ)
    if (allocated(f(id)%cmo)) deallocate(f(id)%cmo)
    if (allocated(f(id)%fused)) deallocate(f(id)%fused)

  end subroutine fields_unload

  !> Calculates the properties of the scalar field f at point x0 and fills the
  !> blanks.
  subroutine fields_propty(id,x0,res,verbose,allfields)
    use struct_basic
    use tools_math
    use tools_io
    use arithmetic
    use types
    implicit none

    integer, intent(in) :: id
    real*8, dimension(:), intent(in)  :: x0
    type(scalar_value), intent(out) :: res
    logical, intent(in) :: verbose, allfields

    real*8 :: xp(3), fres
    integer :: i, j, k

    ! get the scalar field properties
    xp = cr%x2c(x0)
    call grd(f(id),xp,2,res)

    ! r and s
    res%hfevec = res%hf
    call rsindex(res%hfevec,res%hfeval,res%r,res%s)

    if (verbose) then
       if (res%isnuc) then
          write (uout,'("  Type : nucleus")') 
       else
          write (uout,'("  Type : (",a,",",a,")")') string(res%r), string(res%s)
       end if
       write (uout,'("  Field value (f): ",A)') string(res%f,'e',decimal=9)
       write (uout,'("  Field value, valence (fval): ",A)') string(res%fval,'e',decimal=9)
       write (uout,'("  Gradient (grad f): ",3(A,2X))') (string(res%gfort(j),'e',decimal=9),j=1,3)
       write (uout,'("  Gradient norm (|grad f|): ",A)') string(res%gfmodort,'e',decimal=9)
       write (uout,'("  Laplacian (del2 f): ",A)') string(res%del2fort,'e',decimal=9)
       write (uout,'("  Laplacian, valence (del2 fval): ",A)') string(res%del2fval,'e',decimal=9)
       write (uout,'("  Hessian eigenvalues: ",3(A,2X))') (string(res%hfeval(j),'e',decimal=9),j=1,3)
       write (uout,'("  Hessian:")')
       do j = 1, 3
          write (uout,'(4x,1p,3(A,2X))') (string(res%hfort(j,k),'e',decimal=9,length=16,justify=4), k = 1, 3)
       end do
       ! Write ellipticity, if it is a candidate for bond critical point
       if (res%r == 3 .and. res%s == 1 .and..not.res%isnuc) then
          write (uout,'("  Ellipticity (l_1/l_2 - 1): ",A)') string(res%hfeval(1)/res%hfeval(2)-1.d0,'e',decimal=9)
       endif
       if (allfields) then
          do i = 0, mf
             if (fused(i) .and. i /= id) then
                call grd(f(i),xp,2,res)
                write (uout,'("  Field ",A," (f,|grad|,lap): ",3(A,2X))') string(i),&
                   string(res%f,'e',decimal=9), string(res%gfmod,'e',decimal=9), &
                   string(res%del2f,'e',decimal=9)
             end if
          end do
       end if
       ! properties at points defined by the user
       do i = 1, nptprops
          fres = eval_hard_fail(point_prop(i)%expr,xp,fields_fcheck,fields_feval)
          write (uout,'(2X,A," (",A,"): ",A)') string(point_prop(i)%name),&
             string(point_prop(i)%expr), string(fres,'e',decimal=9)
       end do
    end if

  end subroutine fields_propty

  !> Define fields as integrable. Atomic integrals for these fields
  !> will be calculated in the basins of the reference field.
  subroutine fields_integrable(line)
    use global
    use tools_io

    character*(*), intent(in) :: line

    logical :: ok
    integer :: id, lp, lpold, idum
    character(len=:), allocatable :: word, expr, str
    logical :: useexpr

    ! read input
    lp=1
    useexpr = .false.
    ok = eval_next(id,line,lp)
    if (.not.ok) then
       lpold = lp
       ! clear the integrable properties list
       word = lgetword(line,lp)
       if (equal(word,'clear')) then
          nprops = 3
          integ_prop(:)%used = .false.
          integ_prop(1)%used = .true.
          integ_prop(1)%itype = itype_v
          integ_prop(1)%fid = 0
          integ_prop(1)%prop_name = "Volume"
          integ_prop(2)%used = .true.
          integ_prop(2)%itype = itype_fval
          integ_prop(2)%fid = refden
          integ_prop(2)%prop_name = "Pop"
          integ_prop(3)%used = .true.
          integ_prop(3)%itype = itype_lapval
          integ_prop(3)%fid = refden
          integ_prop(3)%prop_name = "Lap"
          call fields_integrable_report()
          return
       else 
          lp = lpold
          if (isexpression_or_word(expr,line,lp)) then
             useexpr = .true.
          else
             call ferror("fields_integrable","Unknown integrable syntax",faterr,line)
          endif
       end if
    end if

    ! add property
    nprops = nprops + 1
    if (nprops > mprops) &
       call ferror("fields_integrable","too many props",faterr)
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
       str = "$" // string(id,2,pad0=.true.) // "_f"
       do while (.true.)
          word = lgetword(line,lp)
          if (equal(word,"f")) then
             integ_prop(nprops)%itype = itype_f
          elseif (equal(word,"fval")) then
             integ_prop(nprops)%itype = itype_fval
             str = "$" // string(id,2,pad0=.true.) // "_fval"
          elseif (equal(word,"gmod")) then
             integ_prop(nprops)%itype = itype_gmod
             str = "$" // string(id,2,pad0=.true.) // "_gmod"
          elseif (equal(word,"lap")) then
             integ_prop(nprops)%itype = itype_lap
             str = "$" // string(id,2,pad0=.true.) // "_lap"
          elseif (equal(word,"lapval")) then
             integ_prop(nprops)%itype = itype_lapval
             str = "$" // string(id,2,pad0=.true.) // "_lval"
          elseif (equal(word,"multipoles") .or. equal(word,"multipole")) then
             integ_prop(nprops)%itype = itype_mpoles
             str = "$" // string(id,2,pad0=.true.) // "_mpol"
             ok = isinteger(idum,line,lp)
             if (ok) integ_prop(nprops)%lmax = idum
          elseif (equal(word,"deloc")) then
             integ_prop(nprops)%itype = itype_deloc
             str = "$" // string(id,2,pad0=.true.) // "_deloc"
          elseif (equal(word,"name")) then
             word = getword(line,lp)
             integ_prop(nprops)%prop_name = string(word)
          else if (len_trim(word) > 0) then
             call ferror("fields_integrable","Unknown extra keyword",faterr,line)
          else
             exit
          end if
       end do
       if (trim(integ_prop(nprops)%prop_name) == "") then
          integ_prop(nprops)%prop_name = str
       end if
    end if

    ! report
    call fields_integrable_report()

  end subroutine fields_integrable

  !> Report for the integrable fields.
  subroutine fields_integrable_report()
    use tools_io

    integer :: i, id
    character*4 :: sprop

    write (uout,'("* List of integrable properties")')
    write (uout,'("# ",4(A,2X))') &
       string("Id",length=3,justify=ioj_right), string("Type",length=4,justify=ioj_center), &
       string("Field",length=5,justify=ioj_right), string("Name")
    do i = 1, nprops
       if (.not.integ_prop(i)%used) cycle
       id = integ_prop(i)%fid
       if (.not.goodfield(id)) cycle

       select case(integ_prop(i)%itype)
       case(itype_v)
          sprop = "v"
       case(itype_f)
          sprop = "f"
       case(itype_fval)
          sprop = "fval"
       case(itype_gmod)
          sprop = "gmod"
       case(itype_lap)
          sprop = "lap"
       case(itype_lapval)
          sprop = "lval"
       case(itype_expr)
          sprop = "expr"
       case(itype_mpoles)
          sprop = "mpol"
       case(itype_deloc)
          sprop = "dloc"
       case default
          call ferror('grdall','unknown property',faterr)
       end select

       if (integ_prop(i)%itype == itype_expr) then
          write (uout,'(2X,2(A,2X),2X,"""",A,"""")') &
             string(i,length=3,justify=ioj_right), string(sprop,length=4,justify=ioj_center), &
             string(integ_prop(i)%expr)
       else
          write (uout,'(2X,4(A,2X))') &
             string(i,length=3,justify=ioj_right), string(sprop,length=4,justify=ioj_center), &
             string(integ_prop(i)%fid,length=5,justify=ioj_right), string(integ_prop(i)%prop_name)
       end if
    end do
    write (uout,*)

  end subroutine fields_integrable_report

  subroutine fields_pointprop(line0)
    use arithmetic
    use global
    use tools_io

    character*(*), intent(in) :: line0

    logical :: isblank
    integer :: lp, n, i
    character(len=:), allocatable :: expr, word, lword, line
    integer, allocatable :: idlist(:)
    
    ! Get the name
    lp = 1
    line = line0
    word = getword(line,lp)
    lword = lower(word)
    if (equal(lword,'clear')) then
       do i = 1, nptprops
          point_prop(i)%name = ""
          point_prop(i)%expr = ""
          point_prop(i)%nf = 0
          if (allocated(point_prop(i)%fused)) deallocate(point_prop(i)%fused)
       end do
       nptprops = 0
       write(uout,'("* No properties to be calculated at points (POINTPROP)"/)')
       return
    elseif (equal(lword,'gtf')) then
       lp = 1 
       line = "gtf(" // string(refden) // ")"
    elseif (equal(lword,'vtf')) then
       lp = 1 
       line = "vtf(" // string(refden) // ")"
    elseif (equal(lword,'htf')) then
       lp = 1 
       line = "htf(" // string(refden) // ")"
    elseif (equal(lword,'gtf_kir')) then
       lp = 1 
       line = "gtf_kir(" // string(refden) // ")"
    elseif (equal(lword,'vtf_kir')) then
       lp = 1 
       line = "vtf_kir(" // string(refden) // ")"
    elseif (equal(lword,'htf_kir')) then
       lp = 1 
       line = "htf_kir(" // string(refden) // ")"
    elseif (equal(lword,'gkin')) then
       lp = 1 
       line = "gkin(" // string(refden) // ")"
    elseif (equal(lword,'kkin')) then
       lp = 1 
       line = "kkin(" // string(refden) // ")"
    elseif (equal(lword,'lag')) then
       lp = 1 
       line = "lag(" // string(refden) // ")"
    elseif (equal(lword,'elf')) then
       lp = 1 
       line = "elf(" // string(refden) // ")"
    elseif (equal(lword,'vir')) then
       lp = 1 
       line = "vir(" // string(refden) // ")"
    elseif (equal(lword,'he')) then
       lp = 1 
       line = "he(" // string(refden) // ")"
    elseif (equal(lword,'lol')) then
       lp = 1 
       line = "lol(" // string(refden) // ")"
    elseif (equal(lword,'lol_kir')) then
       lp = 1 
       line = "lol_kir(" // string(refden) // ")"
    elseif (len_trim(lword) == 0) then
       call fields_pointprop_report()
       return
    endif

    ! Clean up the string
    expr = ""
    isblank = .false.
    do i = lp, len_trim(line)
       if (line(i:i) == " ") then
          if (.not.isblank) then
             expr = expr // line(i:i)
             isblank = .true.
          end if
       elseif (line(i:i) /= """" .and. line(i:i) /= "'") then
          expr = expr // line(i:i)
          isblank = .false.
       endif
    end do
    expr = trim(adjustl(expr))
    if (len_trim(expr) == 0) &
       call ferror("fields_pointprop","Wrong arithmetic expression",faterr,line)

    ! Determine the fields in the expression, check that they are defined
    call fields_in_eval(expr,n,idlist)
    do i = 1, n
       if (.not.goodfield(idlist(i))) &
          call ferror("fields_pointprop","Unknown field in arithmetic expression",faterr,expr)
    end do

    ! Add this pointprop to the list
    nptprops = nptprops + 1
    point_prop(nptprops)%name = word
    point_prop(nptprops)%expr = expr
    point_prop(nptprops)%nf = n
    if (allocated(point_prop(nptprops)%fused)) deallocate(point_prop(nptprops)%fused)
    allocate(point_prop(nptprops)%fused(n))
    point_prop(nptprops)%fused = idlist

    ! Print report
    call fields_pointprop_report()

    ! Clean up
    if (allocated(idlist)) deallocate(idlist)

  end subroutine fields_pointprop

  !> Report for the pointprop fields.
  subroutine fields_pointprop_report()
    use tools_io
  
    integer :: i, id
    character*4 :: sprop
  
    write(uout,'("* List of properties to be calculated at (critical) points (POINTPROP)")')
    write(uout,'("# Id      Name       Expression")')
    do i = 1, nptprops
       write (uout,'(2X,2(A,2X),2X,"""",A,"""")') &
          string(i,length=3,justify=ioj_right), &
          string(point_prop(i)%name,length=10,justify=ioj_center), &
          string(point_prop(i)%expr)
    end do
    write (uout,*)
  
  end subroutine fields_pointprop_report


  !> Returns .true. if the field is initialized and usable.
  logical function goodfield(id,type) result(good)

    integer, intent(in) :: id
    integer, intent(in), optional :: type

    good = (id >= 0 .and. id <= mf)
    if (.not.good) return
    good = fused(id)
    if (.not.good) return
    good = f(id)%init
    if (.not.good) return
    if (present(type)) then
       good = (f(id)%type == type)
       if (.not.good) return
    end if

  end function goodfield

  !> Gives an unused field number
  integer function getfieldnum() result(id)
    use tools_io

    integer :: nmf
    type(field), allocatable :: nf(:)
    logical, allocatable :: nfused(:)

    logical :: found

    found = .false.
    do id = 1, mf
       if (.not.fused(id)) then
          found = .true.
          exit
       end if
    end do
    if (.not.found) then
       nmf = mf + 5
       ! reallocate the fused
       allocate(nfused(0:nmf))
       nfused(0:mf) = fused(0:mf)
       call move_alloc(nfused,fused)
       fused(mf+1:nmf) = .false.
       ! reallocate the fields
       allocate(nf(0:nmf))
       nf(0:mf) = f(0:mf)
       call move_alloc(nf,f)
       ! update the mf and get the first one
       id = mf + 1
       mf = nmf
    end if
    fused(id) = .true.
    
  end function getfieldnum

  subroutine setfield(id,line)
    use struct_basic
    use grid_tools
    use global
    use tools_io
    use param

    integer, intent(in) :: id
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word, aux
    integer :: lp, j
    logical :: ok
    real*8 :: norm

    ! parse the rest of the line
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'trispline')) then
          f(id)%mode = mode_trispline
          if (allocated(f(id)%c2)) deallocate(f(id)%c2)
       else if (equal(word,'trilinear')) then
          f(id)%mode = mode_trilinear
       else if (equal(word,'nearest')) then
          f(id)%mode = mode_nearest
       else if (equal(word,'exact')) then
          f(id)%exact = .true.
       else if (equal(word,'approximate')) then
          f(id)%exact = .false.
       else if (equal(word,'rhonorm')) then
          if (.not.f(id)%cnorm) then
             f(id)%cnorm = .true.
             if (allocated(f(id)%slm)) &
                f(id)%slm(:,1,:) = f(id)%slm(:,1,:) / sqfp
          end if
       else if (equal(word,'vnorm')) then
          if (f(id)%cnorm) then
             f(id)%cnorm = .false.
             if (allocated(f(id)%slm)) &
                f(id)%slm(:,1,:) = f(id)%slm(:,1,:) * sqfp
          end if
       else if (equal(word,'core')) then
          f(id)%usecore = .true.
       else if (equal(word,'nocore')) then
          f(id)%usecore = .false.
       else if (equal(word,'numerical')) then
          f(id)%numerical = .true.
       else if (equal(word,'analytical')) then
          f(id)%numerical = .false.
       else if (equal(word,'typnuc')) then
          ok = eval_next(f(id)%typnuc,line,lp)
          if (.not.ok) call ferror("setfield","wrong typnuc",faterr,line)
          if (f(id)%typnuc /= -3 .and. f(id)%typnuc /= -1 .and. &
              f(id)%typnuc /= +1 .and. f(id)%typnuc /= +3) then
             call ferror("setfield","wrong typnuc (-3,-1,+1,+3)",faterr,line)
          end if
       else if (equal(word,'normalize')) then
          ok = eval_next(norm,line,lp)
          if (.not. ok) call ferror('setfield','normalize real number missing',faterr,line)
          f(id)%f = f(id)%f / (sum(f(id)%f) * cr%omega / real(product(f(id)%n),8)) * norm
          if (allocated(f(id)%c2)) deallocate(f(id)%c2)
       else if (equal(word,'name')) then
          ok = isexpression_or_word(aux,line,lp)
          if (.not. ok) call ferror('setfield','wrong name keyword',faterr,line)
          f(id)%name = aux
       else if (len_trim(word) > 0) then
          call ferror('setfield','Unknown extra keyword',faterr,line)
       else
          exit
       end if
    end do

    write (uout,'("+ Flags for field number: ",A)') string(id)
    if (len_trim(f(id)%name) > 0) then
       write (uout,'("  Name: ",A)') string(f(id)%name)
    else
       write (uout,'("  Name: <field_",A,">")') string(id)
    endif
    if (len_trim(f(id)%file) > 0) then
       write (uout,'("  From: ",A)') string(f(id)%file)
    else
       write (uout,'("  From: <generated>")')
    end if
    write (uout,'("  Type: ",A)') fields_typestring(f(id)%type)
    select case (f(id)%type)
    case(type_grid)
       write (uout,'("  Grid dimensions: ",3(A,2X))') (string(f(id)%n(j)),j=1,3)
       write (uout,'("  Interpolation mode (1=nearest,2=linear,3=spline): ",A)') string(f(id)%mode)
       write (uout,'("  Cell integral (grid SUM) (",A,") = ",A)') &
          string(id), string(sum(f(id)%f) * cr%omega / real(product(f(id)%n),8),'f',decimal=8)
       write (uout,'("  Min: ",A)') string(minval(f(id)%f),'e',decimal=8)
       write (uout,'("  Average: ",A)') string(sum(f(id)%f) / real(product(f(id)%n),8),'e',decimal=8)
       write (uout,'("  Max: ",A)') string(maxval(f(id)%f),'e',decimal=8)
    case(type_wien)
       write (uout,'("  Density-style normalization? ",L)') f(id)%cnorm
    case(type_pi)
       write (uout,'("  Exact calculation? ",L)') f(id)%exact
    end select
    write (uout,'("  Use core densities? ",L)') f(id)%usecore
    write (uout,'("  Numerical derivatives? ",L)') f(id)%numerical
    write (uout,'("  Nuclear CP signature: ",A)') string(f(id)%typnuc)
    write (uout,*)

  end subroutine setfield

  ! Calculate the scalar field f at point v (cartesian) and its derivatives
  ! up to nder. Return the results in res
  subroutine grd(f,v,nder,res)
    use grd_atomic
    use grid_tools
    use struct_basic
    use wfn_private
    use pi_private
    use elk_private
    use wien_private
    use arithmetic
    use types
    use tools_io
    use tools_math

    type(field), intent(inout) :: f
    real*8, intent(in) :: v(3) !< Target point in cartesian or spherical coordinates
    integer, intent(in) :: nder !< Number of derivatives to calculate
    type(scalar_value), intent(out) :: res

    real*8 :: wx(3), wc(3), dist, x(3)
    integer :: i, nid, lvec(3), idx(3)
    real*8 :: rho, grad(3), h(3,3)
    real*8 :: fval(3,-ndif_jmax:ndif_jmax), fzero
    logical :: isgrid

    real*8, parameter :: hini = 1d-3, errcnv = 1d-8
    real*8, parameter :: neargrideps = 1d-12

    if (.not.f%init) call ferror("grd","field not initialized",faterr)

    ! initialize misc quantities 
    res%gkin = 0d0
    res%vir = 0d0

    ! numerical derivatives
    res%isnuc = .false.
    if (f%numerical) then
       fzero = grd0(f,v)
       res%f = fzero
       res%gf = 0d0
       res%hf = 0d0
       if (nder > 0) then
          ! x
          fval(1,:) = 0d0
          fval(1,0) = fzero
          res%gf(1) = der1i((/1d0,0d0,0d0/),v,hini,errcnv,fval(1,:),f,grd0)
          ! y
          fval(2,:) = 0d0
          fval(2,0) = fzero
          res%gf(2) = der1i((/0d0,1d0,0d0/),v,hini,errcnv,fval(2,:),f,grd0)
          ! z
          fval(3,:) = 0d0
          fval(3,0) = fzero
          res%gf(3) = der1i((/0d0,0d0,1d0/),v,hini,errcnv,fval(3,:),f,grd0)
          if (nder > 1) then
             ! xx, yy, zz
             res%hf(1,1) = der2ii((/1d0,0d0,0d0/),v,0.5d0*hini,errcnv,fval(1,:),f,grd0)
             res%hf(2,2) = der2ii((/0d0,1d0,0d0/),v,0.5d0*hini,errcnv,fval(2,:),f,grd0)
             res%hf(3,3) = der2ii((/0d0,0d0,1d0/),v,0.5d0*hini,errcnv,fval(3,:),f,grd0)
             ! xy, xz, yz
             res%hf(1,2) = der2ij((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),v,hini,hini,errcnv,f,grd0)
             res%hf(1,3) = der2ij((/1d0,0d0,0d0/),(/0d0,0d0,1d0/),v,hini,hini,errcnv,f,grd0)
             res%hf(2,3) = der2ij((/0d0,1d0,0d0/),(/0d0,0d0,1d0/),v,hini,hini,errcnv,f,grd0)
             ! final
             res%hf(2,1) = res%hf(1,2)
             res%hf(3,1) = res%hf(1,3)
             res%hf(3,2) = res%hf(2,3)
          end if
       end if
       res%gfmod = norm(res%gf)
       res%del2f = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)
       ! valence quantities
       res%fval = res%f
       res%del2fval = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)
       ! fill orthogonal qtys.
       res%gfort = res%gf
       res%hfort = res%hf
       res%gfmodort = res%gfmod
       res%del2fort = res%del2f
       res%gfsph = 0d0
       res%hfsph = 0d0
       res%gfmodsph = 0d0
       res%del2fsph = 0d0
       return
    end if

    ! To the main cell. Add a small safe zone around the limits of the unit cell
    ! to prevent precision problems.
    wx = cr%c2x(v)
    do i = 1, 3
       if (wx(i) < -flooreps .or. wx(i) > 1d0+flooreps) &
          wx(i) = wx(i) - real(floor(wx(i)),8)
    end do
    wc = cr%x2c(wx)

    ! type selector
    select case(f%type)
    case(type_grid)
       isgrid = .false.
       if (nder == 0) then
          ! maybe we can get the grid point directly
          x = modulo(wx,1d0) * f%n
          idx = nint(x)
          isgrid = all(abs(x-idx) < neargrideps)
       end if
       if (isgrid) then
          idx = modulo(idx,f%n)+1
          res%f = f%f(idx(1),idx(2),idx(3))
          res%gf = 0d0
          res%hf = 0d0
       else
          call grinterp(f,wx,res%f,res%gf,res%hf)
          res%gf = matmul(transpose(cr%car2crys),res%gf)
          res%hf = matmul(matmul(transpose(cr%car2crys),res%hf),cr%car2crys)
       endif

    case(type_wien)
       call wien_rho2(f,wx,res%f,res%gf,res%hf)
       res%gf = matmul(transpose(cr%car2crys),res%gf)
       res%hf = matmul(matmul(transpose(cr%car2crys),res%hf),cr%car2crys)

    case(type_elk)
       call elk_rho2(f,wx,nder,res%f,res%gf,res%hf)
       res%gf = matmul(transpose(cr%car2crys),res%gf)
       res%hf = matmul(matmul(transpose(cr%car2crys),res%hf),cr%car2crys)

    case(type_pi)
       call pi_rho2(f,wc,res%f,res%gf,res%hf)
       ! transformation not needed because of pi_register_struct:
       ! all work done in cartesians in a finite environment.

    case(type_wfn)
       call wfn_rho2(f,wc,nder,res%f,res%gf,res%hf,res%gkin,res%vir)
       ! transformation not needed because all work done in cartesians
       ! in a finite environment. wfn assumes the crystal structure
       ! resulting from load xyz/wfn/wfx (molecule at the center of 
       ! a big cube).

    case(type_promol)
       call grda_promolecular(wx,res%f,res%gf,res%hf,nder,.false.)
       ! not needed because grd_atomic uses struct.

    case(type_promol_frag)
       call grda_promolecular(wx,res%f,res%gf,res%hf,nder,.false.,f%fr)
       ! not needed because grd_atomic uses struct.

    case(type_ghost)
       res%f = eval_hard_fail(f%expr,wc,fields_fcheck,fields_feval)
       res%gf = 0d0
       res%hf = 0d0

    case default
       call ferror("grd","unknown scalar field type",faterr)
    end select

    ! save the valence-only value
    res%fval = res%f
    res%del2fval = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)

    ! augment with the core if applicable
    if (f%usecore .and. any(cr%at(1:cr%nneq)%zpsp /= -1)) then
       call grda_promolecular(wx,rho,grad,h,nder,.true.)
       res%f = res%f + rho
       res%gf  = res%gf + grad
       res%hf = res%hf + h
    end if

    ! If it's on a nucleus, nullify the gradient (may not be zero in
    ! grid fields, for instance)
    nid = 0
    call cr%nearest_atom(wx,nid,dist,lvec)
    res%isnuc = (dist < 1d-5)
    if (res%isnuc) res%gf = 0d0
    res%gfmod = norm(res%gf)
    res%del2f = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)

    ! fill orthogonal qtys.
    res%gfort = res%gf
    res%hfort = res%hf
    res%gfmodort = res%gfmod
    res%del2fort = res%del2f
    res%gfsph = 0d0
    res%hfsph = 0d0
    res%gfmodsph = 0d0
    res%del2fsph = 0d0

  end subroutine grd

  subroutine grdall(xpos,lprop)
    use arithmetic
    use tools_io

    real*8, intent(in) :: xpos(3) !< Point (cartesian).
    real*8, intent(out) :: lprop(nprops) !< the properties vector.

    type(scalar_value) :: res(0:mf)
    logical :: fdone(0:mf)
    integer :: i, id

    lprop = 0d0
    fdone = .false.
    do i = 1, nprops
       if (.not.integ_prop(i)%used) cycle
       if (integ_prop(i)%itype == itype_expr) then
          lprop(i) = eval_hard_fail(integ_prop(i)%expr,xpos,fields_fcheck,fields_feval)
       else
          id = integ_prop(i)%fid
          if (.not.fused(id)) cycle
          if (.not.fdone(id).and.integ_prop(i)%itype /= itype_v) &
             call grd(f(id),xpos,2,res(id))

          select case(integ_prop(i)%itype)
          case(itype_v)
             lprop(i) = 1
          case(itype_f)
             lprop(i) = res(id)%f
          case(itype_fval)
             lprop(i) = res(id)%fval
          case(itype_gmod)
             lprop(i) = res(id)%gfmod
          case(itype_lap)
             lprop(i) = res(id)%del2f
          case(itype_lapval)
             lprop(i) = res(id)%del2fval
          case(itype_mpoles)
             lprop(i) = res(id)%f
          case(itype_deloc)
             lprop(i) = res(id)%f
          case default
             call ferror('grdall','unknown property',faterr)
          end select
       end if
    end do

  end subroutine grdall

  !> Calculate only the value of the scalar field at the given point
  !> (v in cartesian).
  function grd0(f,v)
    use grd_atomic
    use grid_tools
    use struct_basic
    use wfn_private
    use pi_private
    use elk_private
    use wien_private
    use arithmetic
    use types
    use tools_io
    use tools_math
    use param

    type(field), intent(inout) :: f
    real*8, dimension(3), intent(in) :: v !< Target point in cartesian or spherical coordinates.
    real*8 :: grd0

    real*8 :: wx(3), wc(3)
    integer :: i
    real*8 :: h(3,3), grad(3), rho, rhoaux, gkin, vir

    ! To the main cell. Add a small safe zone around the limits of the unit cell
    ! to prevent precision problems.
    wx = cr%c2x(v)
    do i = 1, 3
       if (wx(i) < -flooreps .or. wx(i) > 1d0+flooreps) &
          wx(i) = wx(i) - real(floor(wx(i)),8)
    end do
    wc = cr%x2c(wx)

    ! type selector
    select case(f%type)
    case(type_grid)
       call grinterp(f,wx,rho,grad,h)
    case(type_wien)
       call wien_rho2(f,wx,rho,grad,h)
    case(type_elk)
       call elk_rho2(f,wx,0,rho,grad,h)
    case(type_pi)
       call pi_rho2(f,wc,rho,grad,h)
    case(type_wfn)
       call wfn_rho2(f,wc,0,rho,grad,h,gkin,vir)
    case(type_promol)
       call grda_promolecular(wx,rho,grad,h,0,.false.)
    case(type_promol_frag)
       call grda_promolecular(wx,rho,grad,h,0,.false.,f%fr)
    case(type_ghost)
       rho = eval_hard_fail(f%expr,wc,fields_fcheck,fields_feval)
    case default
       call ferror("grd","unknown scalar field type",faterr)
    end select

    if (f%usecore .and. any(cr%at(1:cr%nneq)%zpsp /= -1)) then
       call grda_promolecular(wx,rhoaux,grad,h,0,.true.)
       rho = rho + rhoaux
    end if
    grd0 = rho

  end function grd0

  !> Return a string description of the field type
  function fields_typestring(itype) result(s)
    use tools_io
    integer, intent(in) :: itype
    character(len=:), allocatable :: s
    
    select case (itype)
    case (type_promol)
       s = "promolecular"
    case (type_grid)
       s = "grid"
    case (type_wien)
       s = "wien2k"
    case (type_elk)
       s = "elk"
    case (type_pi)
       s = "pi"
    case (type_wfn)
       s = "wfn/wfx"
    case (type_promol_frag)
       s = "promolecular fragment"
    case (type_ghost)
       s = "ghost field"
    case default
       call ferror('fields_typestring','unknown field type',faterr)
    end select

  endfunction fields_typestring

  !> Do a benchmark of the speed of the external module in calculating grd and grdall
  !> by using npts random points in the unit cell.
  subroutine benchmark(npts)
    use struct_basic
    use global
    use tools_io
    use types
    use wien_private
    use elk_private
    implicit none

    integer, intent(in) :: npts

    integer :: wpts(cr%nneq+1)
    integer :: i, j
    real*8 :: x(3), aux(3), dist
    integer :: c1, c2, rate
    real*8, allocatable :: randn(:,:,:)
    type(scalar_value) :: res
    logical :: inrmt

    write (uout,'("* Benchmark of the field ")')
    write (uout,'("* Field : ",I2)') refden

    if (f(refden)%type == type_wien .or. f(refden)%type == type_elk) then
       wpts = 0
       allocate(randn(3,npts,cr%nneq+1))

       out: do while (any(wpts(1:cr%nneq+1) < npts))
          call random_number(x)
          do i = 1, cr%ncel
             if (cr%atcel(i)%idx > cr%nneq) cycle
             aux = x - cr%atcel(i)%x
             aux = cr%x2c(aux - nint(aux))
             dist = dot_product(aux,aux)
             if (f(refden)%type == type_wien) then
                inrmt = (dist < wien_rmt_atom(f(refden),cr%at(cr%atcel(i)%idx)%x)**2)
             else
                inrmt = (dist < elk_rmt_atom(f(refden),cr%at(cr%atcel(i)%idx)%x)**2)
             end if
             if (inrmt) then
                if (wpts(cr%atcel(i)%idx) >= npts) cycle out
                wpts(cr%atcel(i)%idx) = wpts(cr%atcel(i)%idx) + 1
                randn(:,wpts(cr%atcel(i)%idx),cr%atcel(i)%idx) = cr%x2c(x)
                cycle out
             end if
          end do
          if (wpts(cr%nneq+1) >= npts) cycle out
          wpts(cr%nneq+1) = wpts(cr%nneq+1) + 1
          randn(:,wpts(cr%nneq+1),cr%nneq+1) = cr%x2c(x)
       end do out

       write (uout,'("* Benchmark of muffin / interstitial grd ")')
       write (uout,'("* Number of points per zone : ",I8)') npts
       do i = 1, cr%nneq+1
          call system_clock(count=c1,count_rate=rate)
          do j = 1, npts
             call grd(f(refden),randn(:,j,i),0,res)
          end do
          call system_clock(count=c2)
          if (i <= cr%nneq) then
             write (uout,'("* Atom : ",I3)') i
          else
             write (uout,'("* Interstitial ")')
          end if
          write (uout,'("* Total wall time : ",F20.6," s ")') real(c2-c1,8) / rate
          write (uout,'("* Avg. wall time per call : ",F20.8," us ")') &
             real(c2-c1,8) / rate / real(npts,8) * 1d6
       end do
       write (uout,*)
       deallocate(randn)

    else
       allocate(randn(3,npts,1))

       call random_number(randn)
       do i = 1, npts
          randn(1:3,i,1) = cr%x2c(randn(1:3,i,1))
       end do

       ! grd
       write (uout,'("  Benchmark of the grd call ")')
       write (uout,'("  Number of points : ",I8)') npts
       call system_clock(count=c1,count_rate=rate)
       do i = 1, npts
          call grd(f(refden),randn(:,i,1),0,res)
       end do
       call system_clock(count=c2)
       write (uout,'("  Total wall time : ",F20.6," s ")') real(c2-c1,8) / rate
       write (uout,'("  Avg. wall time per call : ",F20.8," us ")') &
          real(c2-c1,8) / rate / real(npts,8) * 1d6
       write (uout,*)

       ! ! grdall
       ! write (uout,'("* Benchmark of the grdall call ")')
       ! write (uout,'("* Number of points : ",I8)') npts
       ! call system_clock(count=c1,count_rate=rate)
       ! do i = 1, npts
       !    call grdall(randn(:,i,1),dumprop,dumatom)
       ! end do
       ! call system_clock(count=c2)
       ! write (uout,'("* Total wall time : ",F20.6," s ")') real(c2-c1,8) / rate
       ! write (uout,'("* Avg. wall time per call : ",F20.8," us ")') &
       !    real(c2-c1,8) / rate / real(npts,8) * 1d6
       ! write (uout,*)

       deallocate(randn)

    end if

  end subroutine benchmark

  ! test the muffin tin discontinuity
  subroutine testrmt(id,ilvl)
    use struct_basic
    use global
    use tools_io
    use wien_private
    use elk_private
    use types
    use param

    integer, intent(in) :: id, ilvl

    integer :: n, i, j
    integer :: ntheta, nphi
    integer :: nt, np
    real*8 :: phi, theta, dir(3), xnuc(3), xp(3)
    real*8 :: fin, fout, gfin, gfout
    real*8 :: r, rmt
    character*4 :: label
    character(len=:), allocatable :: linefile
    integer :: luline, luplane
    integer  :: npass(cr%nneq), nfail(cr%nneq)
    logical :: ok
    real*8 :: epsm, epsp, mepsm, mepsp, dif, dosum, mindif, maxdif
    type(scalar_value) :: res

    real*8, parameter :: eps = 1d-3

    if (f(id)%type /= type_wien .and. f(id)%type /= type_elk) return

    write (uout,'("* Muffin-tin discontinuity test")')

    ntheta = 10
    nphi = 10
    do n = 1, cr%nneq
       if (f(id)%type == type_wien) then
          rmt = wien_rmt_atom(f(id),cr%at(n)%x)
       else
          rmt = elk_rmt_atom(f(id),cr%at(n)%x)
       end if
       mepsm = 0d0
       mepsp = 0d0
       if (ilvl > 1) then
          write (uout,'("+ Analysis of the muffin tin discontinuity for atom ",A)') string(n)
          write (uout,'("  ntheta = ",A)') string(ntheta)
          write (uout,'("  nphi = ",A)') string(nphi)
       end if

       xnuc = cr%at(n)%x
       if (ilvl > 1) write (uout,'("  coords = ",3(A,X))') (string(xnuc(j),'f',decimal=9),j=1,3)
       xnuc = cr%x2c(xnuc)

       if (ilvl > 1) then
          write (uout,'("  rmt = ",A)') string(rmt,'f',decimal=7)
          write (uout,'(2(A8,X),6(A12,X),A4)') "Azim.", "Polar", "f_in",&
             "f_out", "f_in-f_out", "gf_in", "gf_out", "gf_in-gf_out", "ok?"
          write (uout,'(100("-"))')
       end if
       npass(n) = 0
       nfail(n) = 0
       dosum = 0d0
       mindif = 1d10
       maxdif = -1d10
       if (ilvl > 1) then
          ! write line
          linefile = "plane_" // string(n,2,pad0=.true.) // ".dbg"
          luplane = fopen_write(linefile)
          write (luplane,'("#",A,I3)') " atom: ", n
          write (luplane,'("#",A,1p,3(E20.13,X))') " at: ", cr%at(n)%x
          write (luplane,'("#",A,1p,3(E20.13,X))') " atc: ", xnuc
          write (luplane,'("#",A,1p,E20.13)') " rmt: ", rmt
          write (luplane,'("#  theta phi in out")')
       end if
       do nt = 1, ntheta
          do np = 0, nphi
             if ((np == 0 .or. np == nphi) .and. nt /= 1) cycle
             phi = real(np,8) * pi / nphi
             theta = real(nt,8) * 2d0 * pi / ntheta
             dir(1) = 1d0 * cos(theta) * sin(phi)
             dir(2) = 1d0 * sin(theta) * sin(phi)
             dir(3) = 1d0 * cos(phi)

             xp = xnuc + (rmt+eps) * dir
             call grd(f(id),xp,1,res)
             fout = res%f
             gfout = dot_product(res%gf,xp-xnuc) / (rmt+eps)
             xp = xnuc + (rmt-eps) * dir
             call grd(f(id),xp,1,res)
             fin = res%f
             gfin = dot_product(res%gf,xp-xnuc) / (rmt-eps)

             dif = fout - fin
             dosum = dosum + dif * dif
             if (dif < mindif) mindif = dif
             if (dif > maxdif) maxdif = dif

             if (ilvl > 1) then
                ! write line
                write (luplane,'(1p,4(E20.13,X))') theta, phi, fin, fout
             end if

             if (gfin*gfout > 0d0) then
                label = "pass"
                npass(n) = npass(n) + 1
             else
                label = "fail"
                if (ilvl > 2) then
                   ! write line
                   linefile = "line_" // string(n,2,pad0=.true.) // "_" // string(nt,3,pad0=.true.) //&
                      "_" // string(np,3,pad0=.true.) // ".dbg"
                   luline = fopen_write(linefile)
                   write (luline,'("#",A,I3)') " atom: ", n
                   write (luline,'("#",A,1p,3(E20.13,X))') " at: ", cr%at(n)%x
                   write (luline,'("#",A,1p,3(E20.13,X))') " atc: ", xnuc
                   write (luline,'("#",A,1p,E20.13)') " rmt: ", rmt
                   write (luline,'("#",A,1p,3(E20.13,X))') " dir: ", dir
                   write (luline,'("#",A,1p,E20.13)') " r_ini: ", 0.50d0 * rmt
                   write (luline,'("#",A,1p,E20.13)') " r_end: ", 4.50d0 * rmt
                   do i = 0, 1000
                      r = 0.50d0 * rmt + (real(i,8) / 1000) * 4d0 * rmt
                      xp = xnuc + r * dir
                      call grd(f(id),xp,1,res)
                      write (luline,'(1p,3(E20.13,X))') r, res%f, dot_product(res%gf,xp-xnuc) / r
                   end do
                   call fclose(luline)
                end if
                epsm = 0d0
                epsp = 0d0
                mepsm = min(epsm,mepsm)
                mepsp = max(epsp,mepsp)
                nfail(n) = nfail(n) + 1
             end if
             if (ilvl > 1) write (uout,'(2(F8.4,X),1p,6(E12.4,X),0p,A4)') &
                theta, phi, fin, fout, fin-fout, gfin, gfout, gfin-gfout, label
          end do
       end do
       dosum = sqrt(dosum / (npass(n)+nfail(n)))
       if (ilvl > 1) then
          write (uout,'(100("-"))')
          write (uout,*)
       end if
       if (ilvl > 1) then
          call fclose(luplane)
       end if
       if (nfail(n) > 0) then
          write (uout,'("  Atom: ",A," delta_m = ",A," delta_p = ",A)') &
             string(n), string(mepsm+1d-3,'f',decimal=6), string(mepsp+1d-3,'f',decimal=6)
          write (uout,*)
       end if
       write (uout,'("  Atom: ",A," RMS/max/min(fout-fin) = ",3(A,2X))') &
          string(n), string(dosum,'f',decimal=6), string(maxdif,'f',decimal=6), &
          string(mindif,'f',decimal=6)
    end do

    ok = .true.
    if (ilvl > 0) then
       write (uout,'("+ Summary ")')
       write (uout,'(A4,3(X,A7))') "Atom", "Pass", "Fail", "Total"
    end if
    do n = 1, cr%nneq
       if (ilvl > 0) write (uout,'(I4,3(X,I7))') n, npass(n), nfail(n), npass(n)+nfail(n)
       ok = ok .and. (nfail(n) == 0)
    end do
    if (ilvl > 0) write (uout,*)
    write (uout,'("+ Assert - no spurious CPs on the muffin tin surface: ",L1/)'), ok
    if (.not.ok) call ferror('testrmt','Spurious CPs on the muffin tin surface!',warning)

  end subroutine testrmt

  !> Calculate the kinetic energy density from the elf
  subroutine taufromelf(ielf,irho,itau)
    use grid_tools
    use tools_io
    use types, only: gk
    use param

    integer, intent(in) :: ielf, irho, itau

    integer :: igrad
    real(kind=gk), allocatable :: g(:,:,:)

    ! check that the fields are good
    if (.not.goodfield(ielf,type_grid)) &
       call ferror("taufromelf","wrong elf field",faterr)
    if (.not.goodfield(irho,type_grid)) &
       call ferror("taufromelf","wrong rho field",faterr)
    if (any(f(ielf)%n /= f(irho)%n)) &
       call ferror("taufromelf","incongruent sizes of rho and elf",faterr)

    ! copy the elf field for now
    if (allocated(f(itau)%c2)) deallocate(f(itau)%c2)
    f(itau) = f(ielf)
    fused(itau) = .true.
    if (allocated(f(itau)%f)) deallocate(f(itau)%f)

    ! allocate a temporary field for the gradient
    igrad = getfieldnum()
    call grid_gradrho(f(irho),f(igrad))
    
    allocate(g(f(ielf)%n(1),f(ielf)%n(2),f(ielf)%n(3)))
    g = sqrt(1d0 / max(min(f(ielf)%f,1d0),1d-14) - 1d0)
    g = g * (3d0/10d0 * (3d0*pi**2)**(2d0/3d0) * max(f(irho)%f,0d0)**(5d0/3d0))
    g = g + 0.125d0 * f(igrad)%f**2 / f(irho)%f
    call move_alloc(g,f(itau)%f)

    ! unload the temporary field
    call fields_unload(igrad)

  end subroutine taufromelf

  !> Check that the id is a grid and is a sane field. Wrapper
  !> around goodfield() to pass it to the arithmetic module.
  function fields_fcheck(id)
    logical :: fields_fcheck
    integer, intent(in) :: id

    fields_fcheck = goodfield(id)

  end function fields_fcheck

  !> Evaluate the field at a point. Wrapper around grd() to pass
  !> it to the arithmetic module. 
  function fields_feval(id,nder,x0)
    use types, only: scalar_value
    type(scalar_value) :: fields_feval
    integer, intent(in) :: id, nder
    real*8, intent(in) :: x0(3)

    call grd(f(id),x0,nder,fields_feval)

  end function fields_feval

end module fields
