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

submodule (systems) proc
  use types, only: thread_info
  implicit none

  !xx! private procedures
  ! function initialization_thread_worker(arg)

contains

  !> Launch the initialization threads, which will go over all systems
  !> trying to initialize them.
  module subroutine launch_initialization_thread()
    use interfaces_threads, only: thrd_create
    use gui_main, only: force_quit_threads
    integer :: i, idum

    force_quit_threads = .false.
#ifdef _THREADS
    do i = 1, nthread
       idum = thrd_create(c_loc(thread(i)), c_funloc(initialization_thread_worker), c_loc(thread_ti(i)))
    end do
#else
    idum = initialization_thread_worker(c_loc(thread_ti(1)))
#endif

  end subroutine launch_initialization_thread

  !> Force the initialization threads to quit and wait for them to
  !> finish before returning.
  module subroutine kill_initialization_thread()
    use interfaces_threads, only: wrap_thrd_join
    use gui_main, only: force_quit_threads
    integer :: i
    integer(c_int) :: res, idum

#ifdef _THREADS
    force_quit_threads = .true.
    do i = 1, nthread
       idum = wrap_thrd_join(c_loc(thread(i)),res)
    end do
    force_quit_threads = .false.
#endif

  end subroutine kill_initialization_thread

  ! Returns .true. if there are any initialization threads running
  module function are_threads_running()
    logical :: are_threads_running

    are_threads_running = any(thread_ti(1:nthread)%active)

  end function are_threads_running

  ! Re-write the seed names from the full-path names to remove as much
  ! of the prefixes as possible
  module subroutine system_shorten_names()
    use param, only: dirsep

    integer :: idx, i, i1
    character(len=:), allocatable :: str, removed

    ! reset all names to the full-path name
    do i = 1, nsys
       if (sysc(i)%status /= sys_empty .and..not.sysc(i)%renamed) &
          sysc(i)%seed%name = sysc(i)%fullname
    end do

    ! remove as much as possible from the beginning
    if (nsys > 0) then
       i1 = 0
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty.and..not.sysc(i)%renamed) then
             i1 = i
             exit
          end if
       end do
       if (i1 > 0) then
          idx = len_trim(sysc(i1)%seed%name)+1
          removed = ""
          main: do while (.true.)
             ! grab string up to the first dirsep
             idx = index(sysc(i1)%seed%name(1:idx-1),dirsep,.true.)
             if (idx == 0) exit main
             str = sysc(i1)%seed%name(1:idx)

             ! check all names start with the same string
             do i = i1+1, nsys
                if (sysc(i)%status == sys_empty.or.sysc(i)%renamed) cycle
                if (len_trim(sysc(i)%seed%name) < idx) exit
                if (sysc(i)%seed%name(1:idx) /= str) cycle main
             end do

             ! remove the string
             removed = removed // sysc(i1)%seed%name(1:idx)
             do i = 1, nsys
                if (sysc(i)%status == sys_empty.or.sysc(i)%renamed) cycle
                sysc(i)%seed%name = sysc(i)%seed%name(idx+1:)
             end do
          end do main
       end if
    end if

  end subroutine system_shorten_names

  !> Add systems from the given seeds. If collapse is present and
  !> true, reorder to make the last system be first and collapse the
  !> other seeds in the tree view. If iafield, load a field from that
  !> seed. If iavib, load vibrational data from that seed. If forceidx
  !> is present, force the system to be in the provided index.
  module subroutine add_systems_from_seeds(nseed,seed,collapse,iafield,iavib,forceidx)
    use utils, only: get_current_working_dir
    use grid1mod, only: grid1_register_ae
    use gui_main, only: reuse_mid_empty_systems
    use windows, only: regenerate_window_pointers
    use interfaces_threads, only: allocate_mtx, mtx_init, mtx_plain
    use crystalseedmod, only: crystalseed
    use tools_io, only: uout
    use types, only: realloc
    use param, only: dirsep, atmcov0
    integer, intent(in) :: nseed
    type(crystalseed), allocatable, intent(in) :: seed(:)
    logical, intent(in), optional :: collapse
    integer, intent(in), optional :: iafield, iavib
    integer, intent(in), optional :: forceidx

    integer :: i, j, nid, idum
    integer :: iafield_, iavib_, forceidx_
    integer :: iseed, iseed_, idx
    character(len=:), allocatable :: errmsg, str
    type(system), allocatable :: syaux(:)
    type(sysconf), allocatable :: syscaux(:)
    logical :: collapse_, isrun
    integer, allocatable :: id(:)

    if (nseed == 0) then
       write (uout,'("WARNING : Could not read structures.")')
       return
    end if

    ! initialize
    errmsg = ""
    collapse_ = .false.
    if (present(collapse)) collapse_ = collapse
    iafield_ = 0
    if (present(iafield)) iafield_ = iafield
    iavib_ = 0
    if (present(iafield)) iavib_ = iavib
    forceidx_ = -1
    if (present(forceidx)) forceidx_ = forceidx

    ! find contiguous IDs for the new systems
    if (nseed == 1 .and. forceidx_ > 0) then
       allocate(id(nseed))
       id(1) = forceidx_
    else
       allocate(id(nseed))
       if (reuse_mid_empty_systems) then
          ! may reuse old IDs
          nid = 0
          do i = 1, nsys
             if (sysc(i)%status == sys_empty) then
                nid = nid + 1
                id(nid) = i
                if (nid == nseed) exit
             else
                nid = 0
             end if
          end do
          do i = nid+1, nseed
             id(i) = nsys + i - nid
          end do
       else
          ! append IDs at the end, discard empty systems
          nid = 0
          do i = nsys, 1, -1
             if (sysc(i)%status /= sys_empty) then
                nid = i
                exit
             end if
          end do
          do i = 1, nseed
             id(i) = nid + i
          end do
       end if
    end if

    ! increment and reallocate if necessary
    nsys = max(nsys,id(nseed))
    if (nsys > size(sys,1)) then
       ! stop the initialization if it's running
       isrun = are_threads_running()
       if (isrun) &
          call kill_initialization_thread()

       allocate(syaux(2*nsys))
       syaux(1:size(sys,1)) = sys
       call move_alloc(syaux,sys)

       allocate(syscaux(2*nsys))
       syscaux(1:size(sysc,1)) = sysc
       call move_alloc(syscaux,sysc)

       ! refresh system and window pointers
       call regenerate_system_pointers()
       call regenerate_window_pointers()

       ! start or restart initializations
       if (isrun) &
          call launch_initialization_thread()
    end if

    do iseed = 1, nseed
       ! re-order if collapse to have the final structure be first, flag them
       if (collapse_) then
          iseed_ = mod(iseed,nseed)+1
       else
          iseed_ = iseed
       end if

       ! make a new system for this seed
       idx = id(iseed_)
       sys(idx)%isinit = .false.
       sysc(idx)%id = idx
       sysc(idx)%seed = seed(iseed)
       sysc(idx)%idseed = iseed
       sysc(idx)%has_field = .false.
       sysc(idx)%has_vib = .false.
       sysc(idx)%renamed = .false.
       sysc(idx)%showfields = .false.
       if (allocated(sysc(idx)%highlight_rgba)) deallocate(sysc(idx)%highlight_rgba)
       if (allocated(sysc(idx)%highlight_rgba_transient)) deallocate(sysc(idx)%highlight_rgba_transient)

       ! write down the full name
       str = trim(adjustl(sysc(idx)%seed%name))
       if (str(1:1) == dirsep) then
          sysc(idx)%fullname = str
       else
          sysc(idx)%fullname = get_current_working_dir() // dirsep // str
       end if

       ! initialization status
       sysc(idx)%status = sys_loaded_not_init
       if (collapse_.and.iseed_ == 1) then
          ! master
          sysc(idx)%collapse = -1
          sysc(idx)%hidden = .false.
       elseif (collapse_) then
          ! dependent
          sysc(idx)%collapse = id(1)
          sysc(idx)%hidden = .true.
       else
          ! independent
          sysc(idx)%collapse = 0
          sysc(idx)%hidden = .false.
       end if

       ! initialize the mutex
#ifdef _THREADS
       if (.not.c_associated(sysc(idx)%thread_lock)) then
          sysc(idx)%thread_lock = allocate_mtx()
          idum = mtx_init(sysc(idx)%thread_lock,mtx_plain)
       end if
#endif

       ! register all all-electron densities (global - should not
       ! be done by threads because agrid and cgrid are common)
       do j = 1, seed(iseed)%nspc
          call grid1_register_ae(seed(iseed)%spc(j)%z)
       end do

       ! set the iafield
       if (iseed == iafield_) sysc(idx)%has_field = .true.

       ! set the iafield
       if (iseed == iavib_) sysc(idx)%has_vib = .true.

       ! redo everything
       call sysc(idx)%post_event(lastchange_geometry)

       ! set the covalent radii
       sysc(idx)%atmcov = atmcov0
    end do
    deallocate(id)

  end subroutine add_systems_from_seeds

  !> Add systems by reading them from a file, passed by name. mol = 1
  !> read as crystal, 0 read as molecule, -1 autodetect. isformat,
  !> force file format if /= 0. If forceidx is present, force the
  !> system to be in the provided index.
  module subroutine add_systems_from_name(name,mol,isformat,readlastonly,rborder,molcubic,&
     forceidx)
    use crystalseedmod, only: crystalseed, read_seeds_from_file
    use tools_io, only: uout
    character(len=*), intent(in) :: name
    integer, intent(in) :: mol
    integer, intent(in) :: isformat
    logical, intent(in) :: readlastonly
    real*8, intent(in) :: rborder
    logical, intent(in) :: molcubic
    integer, intent(in), optional :: forceidx

    integer :: i, nseed
    type(crystalseed), allocatable :: seed(:)
    integer :: iafield, iavib
    logical :: collapse
    character(len=:), allocatable :: errmsg

    ! read all seeds from the file
    errmsg = ""
    nseed = 0
    allocate(seed(1))
    call read_seeds_from_file(name,mol,isformat,readlastonly,nseed,seed,collapse,errmsg,&
       iafield,iavib)
    if (len_trim(errmsg) > 0 .or. nseed == 0) goto 999

    ! set the border and molcubic
    do i = 1, nseed
       seed(i)%border = rborder
       seed(i)%cubic = molcubic
    end do

    ! add the systems
    call add_systems_from_seeds(nseed,seed,collapse,iafield,iavib,forceidx=forceidx)

    return
999 continue
    write (uout,'("WARNING : Could not read structures from: ",A)') trim(name)
    write (uout,'("WARNING : ",A/)') trim(errmsg)

  end subroutine add_systems_from_name

  ! Remove system with index idx and leave behind a sys_empty spot. If
  ! master and collapsed, kill all dependents. If master and extended,
  ! make all dependents master.
  recursive module subroutine remove_system(idx,kill_dependents_if_extended)
    use interfaces_threads, only: deallocate_mtx
    integer, intent(in) :: idx
    logical, intent(in), optional :: kill_dependents_if_extended

    integer :: i
    logical :: kdie

    kdie = .false.
    if (present(kill_dependents_if_extended)) kdie = kill_dependents_if_extended

    if (.not.ok_system(idx,sys_loaded_not_init)) return
    call sys(idx)%end()
    call sysc(idx)%seed%end()
#ifdef _THREADS
    if (c_associated(sysc(idx)%thread_lock)) then
       call deallocate_mtx(sysc(idx)%thread_lock)
       sysc(idx)%thread_lock = c_null_ptr
    end if
#endif
    sysc(idx)%status = sys_empty
    sysc(idx)%hidden = .false.
    sysc(idx)%showfields = .false.
    call sysc(idx)%post_event(lastchange_geometry)

    if (sysc(idx)%collapse == -1 .or. kdie) then
       ! kill all dependents if collapsed
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty .and. sysc(i)%collapse == idx) &
             call remove_system(i)
       end do
    else
       ! make all dependents their own master if extended
       do i = 1, nsys
          if (sysc(i)%status /= sys_empty .and. sysc(i)%collapse == idx) &
             sysc(i)%collapse = 0
       end do
    end if

  end subroutine remove_system

  !> Remove systems given by indices idx.
  module subroutine remove_systems(idx)
    integer, intent(in) :: idx(:)

    integer :: k

    do k = 1, size(idx,1)
       call remove_system(idx(k))
    end do

  end subroutine remove_systems

  !> Re-read the system from file name of system idx
  module subroutine reread_system_from_file(idx)
    use crystalseedmod, only: crystalseed, read_seeds_from_file, realloc_crystalseed
    use param, only: isformat_r_from_library
    integer, intent(in) :: idx

    character(len=:), allocatable :: file, errmsg
    logical :: exist
    integer :: isformat, mol, iafield, iavib, nseed, icol, idseed
    type(crystalseed), allocatable :: seed(:)
    logical :: collapse, ihid, renamed
    character(len=:), allocatable :: fullname

    ! the seed is available
    if (.not.ok_system(idx,sys_loaded_not_init)) return

    ! save the seed info
    if (sysc(idx)%seed%ismolecule) then
       mol = 1
    else
       mol = 0
    end if
    isformat = sysc(idx)%seed%isformat
    icol = sysc(idx)%collapse
    ihid = sysc(idx)%hidden
    renamed = sysc(idx)%renamed
    fullname = sysc(idx)%fullname
    idseed = sysc(idx)%idseed

    ! make sure the file exists
    file = sysc(idx)%seed%file
    if (isformat /= isformat_r_from_library) then
       inquire(file=file,exist=exist)
       if (.not.exist) return
    end if

    ! terminate the system
    call sys(idx)%end()
    call sysc(idx)%seed%end()
    sysc(idx)%status = sys_empty

    ! re-read the seeds from file
    errmsg = ""
    nseed = 0
    allocate(seed(1))
    if (isformat == isformat_r_from_library) then
       call seed(1)%read_library(file,exist)
       if (.not.exist) return
    else
       call read_seeds_from_file(file,mol,isformat,.false.,nseed,seed,collapse,errmsg,&
          iafield,iavib)
       if (len_trim(errmsg) > 0 .or. nseed == 0) return

       ! move the relevant seed to the first position
       if (nseed > 1) then
          seed(1) = seed(sysc(idx)%idseed)
          call realloc_crystalseed(seed,1)
       end if
    end if

    ! add the system again
    call add_systems_from_seeds(1,seed,.false.,iafield,iavib,idx)
    sysc(idx)%collapse = icol
    sysc(idx)%hidden = ihid
    sysc(idx)%renamed = renamed
    sysc(idx)%fullname = fullname
    sysc(idx)%idseed = idseed

    ! initialize
    call launch_initialization_thread()
    call system_shorten_names()

  end subroutine reread_system_from_file

  !> Duplicate the given system.
  module subroutine duplicate_system(idx)
    use crystalseedmod, only: crystalseed
    integer, intent(in) :: idx

    type(crystalseed), allocatable :: seed(:)

    if (.not.ok_system(idx,sys_loaded_not_init)) return
    allocate(seed(1))
    seed(1) = sysc(idx)%seed
    call add_systems_from_seeds(1,seed)
    call launch_initialization_thread()

  end subroutine duplicate_system

  ! Write the system idx to the same file it was read from. Only for
  ! those file formats that are not outputs.
  module subroutine write_system(idx)
    use tools_io, only: uout
    use param, only: isformat_write_from_read, isformat_w_unknown
    integer, intent(in) :: idx

    integer :: iwformat
    character(len=:), allocatable :: errmsg

    ! only for initialized systems
    if (.not.ok_system(idx,sys_init)) return

    ! convert to the write format
    iwformat = isformat_write_from_read(sysc(idx)%seed%isformat)
    if (iwformat == isformat_w_unknown) return

    ! write the file
    call sys(idx)%c%write_any_file(sysc(idx)%fullname,errmsg,iwformat=iwformat)
    if (len_trim(errmsg) > 0) then
       write (uout,'("WARNING : Could not write structures: ",A)') trim(sysc(idx)%fullname)
       write (uout,'("WARNING : ",A/)') trim(errmsg)
    end if

  end subroutine write_system

  !> This routine regenerates all pointers to the sytsems in the
  !> sys(:) and sysc(:) structures and components. It is used when an
  !> array size is exceeded and move_alloc needs to be used to
  !> allocate more memory.
  module subroutine regenerate_system_pointers()
    use fieldmod, only: type_wfn, type_pi, type_dftb, type_elk, type_grid

    integer :: i, j

    ! refresh pointers in the systems after the move_alloc
    do i = 1, nsys
       do j = 1, sys(i)%nf
          sys(i)%f(j)%sptr = c_loc(sys(i))
          sys(i)%f(j)%c => sys(i)%c

          ! reallocate environments; clean this up when that part changes
          if (sys(i)%f(j)%type == type_wfn) then
             sys(i)%f(j)%wfn%cptr = c_loc(sys(i)%c)
          elseif (sys(i)%f(j)%type == type_pi) then
             sys(i)%f(j)%pi%cptr = c_loc(sys(i)%c)
          elseif (sys(i)%f(j)%type == type_dftb) then
             sys(i)%f(j)%dftb%cptr = c_loc(sys(i)%c)
          elseif (sys(i)%f(j)%type == type_elk) then
             sys(i)%f(j)%elk%cptr = c_loc(sys(i)%c)
          elseif (sys(i)%f(j)%type == type_grid) then
             sys(i)%f(j)%grid%cptr = c_loc(sys(i)%c)
          end if
       end do
    end do

  end subroutine regenerate_system_pointers

  !> Check that the system ID is sane and has at least the requested
  !> initialization level.
  module function ok_system(isys,level)
    integer, intent(in) :: isys, level
    logical :: ok_system

    ok_system = (isys >= 1 .and. isys <= nsys)
    if (ok_system) ok_system = (sysc(isys)%status >= level)

  end function ok_system

  !> Set the time for last change at level level. If keepfields is
  !> present and true, do not reset the associated fields.
  module subroutine post_event(sysc,level,keepfields)
    use interfaces_glfw, only: glfwGetTime
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: level
    logical, intent(in), optional :: keepfields

    real*8 :: time
    logical :: keepfields_

    keepfields_ = .false.
    if (present(keepfields)) keepfields_ = keepfields

    time = glfwGetTime()
    if (level >= lastchange_render) sysc%timelastchange_render = time
    if (level >= lastchange_buildlists) sysc%timelastchange_buildlists = time
    if (level >= lastchange_rebond) sysc%timelastchange_rebond = time
    if (level >= lastchange_geometry) then
       if (.not.keepfields_) &
          call sys(sysc%id)%reset_fields()
       call sysc%highlight_clear(.true.)
       call sysc%highlight_clear(.false.)
       sysc%timelastchange_geometry = time
    end if

  end subroutine post_event

  !> Highlight atoms in the system. If transient, add the highlighted
  !> atom to the transient list and clear the list before adding. Add
  !> atoms with indices idx of the given type (atlisttype_*) and color
  !> rgba.
  module subroutine highlight_atoms(sysc,transient,idx,type,rgba)
    class(sysconf), intent(inout) :: sysc
    logical, intent(in) :: transient
    integer, intent(in) :: idx(:)
    integer, intent(in) :: type
    real(c_float), intent(in) :: rgba(:,:)

    integer :: nat, iat
    real(c_float), allocatable :: highlight_aux(:,:)
    logical :: changed
    integer, allocatable :: imask(:)

    ! input checks
    if (sysc%status < sys_init) return
    nat = sys(sysc%id)%c%ncel
    if (nat <= 0) return

    ! check the size of the highlight arrays, allocate if not already done
    if (.not.transient) then
       if (allocated(sysc%highlight_rgba)) then
          if (size(sysc%highlight_rgba,2) /= nat) deallocate(sysc%highlight_rgba)
       end if
       if (.not.allocated(sysc%highlight_rgba)) then
          allocate(sysc%highlight_rgba(4,nat))
          sysc%highlight_rgba = -1._c_float
          call sysc%post_event(lastchange_render)
       end if
    else
       if (allocated(sysc%highlight_rgba_transient)) then
          if (size(sysc%highlight_rgba_transient,2) /= nat) deallocate(sysc%highlight_rgba_transient)
       end if
       if (.not.allocated(sysc%highlight_rgba_transient)) then
          allocate(sysc%highlight_rgba_transient(4,nat))
          sysc%highlight_rgba_transient = -1._c_float
          call sysc%post_event(lastchange_render)
       end if
    end if

    ! save a copy of the current highlight
    allocate(highlight_aux(4,nat))
    if (transient) then
       highlight_aux = sysc%highlight_rgba_transient
    else
       highlight_aux = sysc%highlight_rgba
    end if

    ! if transient, clear the list
    if (transient) sysc%highlight_rgba_transient = -1

    ! highlight the atoms
    call sysc%attype_celatom_mask(type,idx,imask=imask)
    do iat = 1, nat
       if (imask(iat) == 0) cycle

       if (.not.transient) then
          sysc%highlight_rgba(:,iat) = rgba(:,imask(iat))
       else
          sysc%highlight_rgba_transient(:,iat) = rgba(:,imask(iat))
       end if
    end do

    ! set highlighted, the transient flag, and time for render if highlight has changed
    if (transient) then
       changed = any(sysc%highlight_rgba_transient /= highlight_aux)
       sysc%highlight_transient_set = .true.
    else
       changed = any(sysc%highlight_rgba /= highlight_aux)
    end if
    if (changed) then
       call sysc%post_event(lastchange_render)
    end if

  end subroutine highlight_atoms

  !> Clear highlights in the system. If transient, clear the
  !> highlighted atoms in the transient list. If idx and type are
  !> present, clear atoms with indices idx of the given type
  !> (0=species,1=nneq,2=ncel,3=nmol).
  module subroutine highlight_clear(sysc,transient,idx,type)
    class(sysconf), intent(inout) :: sysc
    logical, intent(in) :: transient
    integer, intent(in), optional :: idx(:)
    integer, intent(in), optional :: type

    integer :: nat, iat, id, i
    logical :: changed
    real(c_float), allocatable :: highlight_aux(:,:)

    ! input checks
    if (sysc%status < sys_init) return
    nat = sys(sysc%id)%c%ncel
    if (nat <= 0) return

    ! check the size of the highlight arrays, allocate if not already done
    if (.not.transient) then
       if (allocated(sysc%highlight_rgba)) then
          if (size(sysc%highlight_rgba,2) /= nat) deallocate(sysc%highlight_rgba)
       end if
       if (.not.allocated(sysc%highlight_rgba)) then
          allocate(sysc%highlight_rgba(4,nat))
          sysc%highlight_rgba = -1._c_float
          call sysc%post_event(lastchange_render)
       end if
    else
       if (allocated(sysc%highlight_rgba_transient)) then
          if (size(sysc%highlight_rgba_transient,2) /= nat) deallocate(sysc%highlight_rgba_transient)
       end if
       if (.not.allocated(sysc%highlight_rgba_transient)) then
          allocate(sysc%highlight_rgba_transient(4,nat))
          sysc%highlight_rgba_transient = -1._c_float
          call sysc%post_event(lastchange_render)
       end if
    end if

    ! clear highlights
    if (transient) then
       changed = any(sysc%highlight_rgba_transient >= 0._c_float)
       sysc%highlight_rgba_transient = -1._c_float
    else
       if (present(idx) .and. present(type)) then
          ! save a copy of the current highlight
          allocate(highlight_aux(4,nat))
          highlight_aux = sysc%highlight_rgba

          ! highlight the atoms
          do iat = 1, nat
             if (type == 0) then ! species
                id = sys(sysc%id)%c%atcel(iat)%is
             elseif (type == 1) then ! nneq
                id = sys(sysc%id)%c%atcel(iat)%idx
             elseif (type == 2) then ! ncel
                id = iat
             else ! nmol
                id = sys(sysc%id)%c%idatcelmol(1,iat)
             end if

             do i = 1, size(idx,1)
                if (idx(i) == id) then
                   if (.not.transient) then
                      sysc%highlight_rgba(:,iat) = -1._c_float
                   else
                      sysc%highlight_rgba_transient(:,iat) = -1._c_float
                   end if
                end if
             end do
          end do
          changed = any(sysc%highlight_rgba /= highlight_aux)
       else
          changed = any(sysc%highlight_rgba >= 0._c_float)
          sysc%highlight_rgba = -1._c_float
       end if
    end if

    ! set highlighted and time for render
    if (changed) then
       call sysc%post_event(lastchange_render)
    end if

  end subroutine highlight_clear

  !> Remove the highlighted atoms from the system.
  module subroutine remove_highlighted_atoms(sysc)
    class(sysconf), intent(inout) :: sysc

    integer :: i, nat, id
    integer, allocatable :: iat(:)

    ! consistency checks
    id = sysc%id
    if (.not.ok_system(id,sys_init)) return

    ! build the list of atoms to remove
    allocate(iat(sys(id)%c%ncel))
    nat = 0
    do i = 1, sys(id)%c%ncel
       if (any(sysc%highlight_rgba(:,i) >= 0._c_float)) then
          nat = nat + 1
          iat(nat) = i
       end if
    end do

    ! return if nothing to do
    if (nat == 0) return

    ! remove the atoms, reset fields
    call sys(id)%c%delete_atoms(nat,iat(1:nat))

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine remove_highlighted_atoms

  ! Show a simple combo using the allowed atom types. Type is the
  ! initial and final atom type in the combo.
  module function attype_combo_simple(sysc,label,type,allowed,units)
    use utils, only: iw_combo_simple
    class(sysconf), intent(inout) :: sysc
    character(len=*), intent(in) :: label
    integer, intent(inout) :: type
    integer, intent(in) :: allowed(:)
    logical, intent(in), optional :: units
    logical :: attype_combo_simple

    character(len=:), allocatable :: strcombo
    logical :: allowedin(atlisttype_NUM)
    integer :: i, id, icombo, icount, itype
    integer :: decode(0:atlisttype_NUM-1)
    logical(c_bool) :: ch
    logical :: units_

    ! consistency checks and initialize
    attype_combo_simple = .false.
    id = sysc%id
    if (.not.ok_system(id,sys_init)) return
    units_ = .true.
    if (present(units)) units_ = units

    ! build the allowed mask
    allowedin = .false.
    do i = 1, size(allowed,1)
       allowedin(allowed(i)) = .true.
    end do

    ! deactivate the unallowed types
    if (sys(id)%c%ismolecule) then
       allowedin(atlisttype_nneq) = .false.
       allowedin(atlisttype_ncel_frac) = .false.
    end if

    ! if the input type is not allowed, reset it
    itype = -1
    if (.not.allowedin(type)) then
       if (sys(id)%c%ismolecule) then
          if (allowedin(atlisttype_ncel_ang)) then
             itype = atlisttype_ncel_ang
          elseif (allowedin(atlisttype_ncel_bohr)) then
             itype = atlisttype_nneq
          end if
       else
          if (allowedin(atlisttype_ncel_frac)) then
             itype = atlisttype_ncel_frac
          elseif (allowedin(atlisttype_ncel_ang)) then
             itype = atlisttype_ncel_ang
          elseif (allowedin(atlisttype_ncel_ang)) then
             itype = atlisttype_nneq
          end if
       end if
    else
       itype = type
    end if
    if (itype < 0) then
       do i = 1, size(allowedin,1)
          if (allowedin(i)) then
             itype = i
             exit
          end if
       end do
    end if
    if (itype < 0) return

    ! build the combo string
    strcombo = ""
    if (allowedin(atlisttype_species)) strcombo = strcombo // "Species" // c_null_char
    if (allowedin(atlisttype_nneq)) strcombo = strcombo // "Symmetry unique" // c_null_char
    if (sys(id)%c%ismolecule) then
       if (allowedin(atlisttype_ncel_bohr)) then
          if (units_) then
             strcombo = strcombo // "Atoms (bohr)" // c_null_char
          else
             strcombo = strcombo // "Atoms" // c_null_char
          end if
       end if
       if (allowedin(atlisttype_ncel_ang)) then
          if (units_) then
             strcombo = strcombo // "Atoms (Å)" // c_null_char
          else
             strcombo = strcombo // "Atoms" // c_null_char
          end if
       end if
    else
       if (allowedin(atlisttype_ncel_frac)) then
          if (units_) then
             strcombo = strcombo // "Cell (fractional)" // c_null_char
          else
             strcombo = strcombo // "Cell" // c_null_char
          end if
       end if
       if (allowedin(atlisttype_ncel_bohr)) then
          if (units_) then
             strcombo = strcombo // "Cell (bohr)" // c_null_char
          else
             strcombo = strcombo // "Cell" // c_null_char
          end if
       end if
       if (allowedin(atlisttype_ncel_ang)) then
          if (units_) then
             strcombo = strcombo // "Cell (Å)" // c_null_char
          else
             strcombo = strcombo // "Cell" // c_null_char
          end if
       end if
    end if
    if (allowedin(atlisttype_nmol)) strcombo = strcombo // "Molecules" // c_null_char
    strcombo = strcombo // c_null_char // c_null_char

    ! encode the combo
    icombo = -1
    icount = -1
    do i = 1, size(allowedin,1)
       if (allowedin(i)) then
          icount = icount + 1
          decode(icount) = i
          if (itype == i) icombo = icount
       end if
    end do

    ! the actual combo
    call iw_combo_simple(label,strcombo,icombo,changed=ch)
    attype_combo_simple = ch

    ! decode the combo
    type = decode(icombo)

  end function attype_combo_simple

  ! For the given atom type, return the number of entities (atoms,
  ! species, molecules, etc.)
  module function attype_number(sysc,type)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer :: attype_number

    integer :: isys

    ! initialize
    attype_number = 0
    isys = sysc%id

    ! get the atom number
    if (type == atlisttype_species) then
       attype_number = sys(isys)%c%nspc
    elseif (type == atlisttype_nneq) then
       attype_number = sys(isys)%c%nneq
    elseif (type == atlisttype_nmol) then
       attype_number = sys(isys)%c%nmol
    else
       attype_number = sys(isys)%c%ncel
    end if

  end function attype_number

  ! For the given atom type, return the ID of the corresponding
  ! atomic species. If this cannot be satisfied (e.g. type is nmol),
  ! return zero.
  module function attype_species(sysc,type,id)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    integer :: attype_species

    integer :: isys

    ! initialize
    attype_species = 0
    isys = sysc%id

    ! get the species ID
    if (type == atlisttype_species) then
       attype_species = id
    elseif (type == atlisttype_nneq) then
       attype_species = sys(isys)%c%at(id)%is
    elseif (type /= atlisttype_nmol) then
       attype_species = sys(isys)%c%atcel(id)%is
    end if

  end function attype_species

  ! For the given atom type, return the corresponding atom name.
  module function attype_name(sysc,type,id)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    character(len=:), allocatable :: attype_name

    integer :: isys

    ! initialize
    attype_name = ""
    isys = sysc%id

    ! get the atom name
    if (type == atlisttype_species) then
       attype_name = trim(sys(isys)%c%spc(id)%name)
    elseif (type == atlisttype_nneq) then
       attype_name = trim(sys(isys)%c%at(id)%name)
    elseif (type /= atlisttype_nmol) then
       attype_name = trim(sys(isys)%c%at(sys(isys)%c%atcel(id)%idx)%name)
    end if

  end function attype_name

  ! For the given atom type, set the corresponding atom name.
  module subroutine set_attype_name(sysc,type,id,str)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    character(len=*), intent(in) :: str

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    if (type == atlisttype_species) then
       sys(isys)%c%spc(id)%name = str
    elseif (type == atlisttype_nneq) then
       sys(isys)%c%at(id)%name = str
    elseif (type /= atlisttype_nmol) then
       sys(isys)%c%at(sys(isys)%c%atcel(id)%idx)%name = str
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_attype_name

  ! For the given atom type, return the corresponding atomic
  ! coordinates.
  module function attype_coordinates(sysc,type,id)
    use param, only: bohrtoa
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    real*8 :: attype_coordinates(3)

    integer :: isys

    ! initialize
    attype_coordinates = 0d0
    isys = sysc%id

    ! get the atom name
    if (type == atlisttype_nneq) then
       attype_coordinates = sys(isys)%c%at(id)%x
    elseif (type == atlisttype_ncel_frac) then
       attype_coordinates = sys(isys)%c%atcel(id)%x
    elseif (type == atlisttype_ncel_bohr) then
       if (sys(isys)%c%ismolecule) then
          attype_coordinates = sys(isys)%c%atcel(id)%r + sys(isys)%c%molx0
       else
          attype_coordinates = sys(isys)%c%atcel(id)%r
       end if
    elseif (type == atlisttype_ncel_ang) then
       if (sys(isys)%c%ismolecule) then
          attype_coordinates = (sys(isys)%c%atcel(id)%r + sys(isys)%c%molx0) * bohrtoa
       else
          attype_coordinates = sys(isys)%c%atcel(id)%r * bohrtoa
       end if
    end if

  end function attype_coordinates

  ! For the given atom type, return the decimals with which to
  ! represent the atomic coordinates.
  module function attype_coordinates_decimals(sysc,type)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer :: attype_coordinates_decimals

    if (type == atlisttype_nneq .or. type == atlisttype_ncel_frac) then
       attype_coordinates_decimals = 6
    else
       attype_coordinates_decimals = 4
    end if

  end function attype_coordinates_decimals

  ! If mask is present, return a logical mask of length equal to the
  ! number of atoms in the cell, with .true. if the atom matches one
  ! of the given ids for the corresponding atom type. If imask is
  ! present, return a mask of length equal to the number of atoms in the
  ! cell. For a given element i in the mask, imask(i) is zero if this is
  ! not one of the atoms in ids for the corresponding type. If it is,
  ! imask(i) has type ids(imask(i)).
  module subroutine attype_celatom_mask(sysc,type,ids,mask,imask)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: ids(:)
    logical, allocatable, intent(inout), optional :: mask(:)
    integer, allocatable, intent(inout), optional :: imask(:)

    integer :: nat, i, j, id
    logical :: ok

    ! initialize
    if (present(mask)) then
       if (allocated(mask)) deallocate(mask)
    end if
    if (present(imask)) then
       if (allocated(imask)) deallocate(imask)
    end if

    ! consistency checks
    id = sysc%id
    if (.not.ok_system(id,sys_init)) return

    ! initialize mask
    nat = sys(id)%c%ncel
    if (nat <= 0) return
    if (present(mask)) then
       allocate(mask(nat))
       mask = .false.
    end if
    if (present(imask)) then
       allocate(imask(nat))
       imask = 0
    end if

    do i = 1, size(ids,1)
       do j = 1, nat
          if (type == atlisttype_species) then
             ok = sys(id)%c%atcel(j)%is == ids(i)
          elseif (type == atlisttype_nneq) then
             ok = sys(id)%c%atcel(j)%idx == ids(i)
          elseif (type == atlisttype_nmol) then
             ok = sys(id)%c%idatcelmol(1,j) == ids(i)
          else
             ok = (j == ids(i))
          end if
          if (ok) then
             if (present(mask)) mask(j) = .true.
             if (present(imask)) imask(j) = i
          end if
       end do
    end do

  end subroutine attype_celatom_mask

  ! Given the celatom number id, return the corresponding id
  ! corresponding to the given type.
  module function attype_celatom_to_id(sysc,type,id)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    integer :: attype_celatom_to_id

    integer :: isys

    ! initialize
    attype_celatom_to_id = 0
    isys = sysc%id
    if (id < 1 .or. id > sys(isys)%c%ncel) return

    if (type == atlisttype_species) then
       attype_celatom_to_id = sys(isys)%c%atcel(id)%is
    elseif (type == atlisttype_nneq) then
       attype_celatom_to_id = sys(isys)%c%atcel(id)%idx
    elseif (type == atlisttype_nmol) then
       attype_celatom_to_id = sys(isys)%c%idatcelmol(1,id)
    else
       attype_celatom_to_id = id
    end if

  end function attype_celatom_to_id

  ! Given the identifier id belonging to type typein, return the only
  ! ID corresponding to type type that matches it. Example: if given
  ! non-equivalent atom ID 3, return the species ID of the species
  ! corresponding to that atom. If the result is ambiguous, return
  ! 0. Example of ambiguous result: return the non-equivalent atom ID
  ! corresponding to species 3.
  module function attype_type_id_to_id(sysc,typein,id,typeout)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: typein
    integer, intent(in) :: id
    integer, intent(in) :: typeout
    integer :: attype_type_id_to_id

    ! initialize and return the same ID if the types are equal
    attype_type_id_to_id = 0
    if (typein == typeout) then
       attype_type_id_to_id = id
       return
    end if
    if (typeout == atlisttype_species) then
       attype_type_id_to_id = sysc%attype_species(typein,id)
    elseif (typein == atlisttype_ncel_frac.or.typein == atlisttype_ncel_bohr.or.&
       typein == atlisttype_ncel_ang) then
       if (typeout == atlisttype_nneq) then
          attype_type_id_to_id = sys(sysc%id)%c%atcel(id)%idx
       elseif (typeout == atlisttype_nmol) then
          attype_type_id_to_id = sys(sysc%id)%c%idatcelmol(1,id)
       else
          attype_type_id_to_id = id
       end if
    end if

  end function attype_type_id_to_id

  ! For the atom identifier id corresponding to the given atom type,
  ! set the atomic position(s) in the system.
  module subroutine set_atom_position(sysc,type,id,x,forcewyc)
    use global, only: iunit_bohr, iunit_fractional
    use param, only: bohrtoa
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    real*8, intent(in) :: x(3)
    logical, intent(in) :: forcewyc

    integer :: isys, leqv, i, ichange
    real*8 :: x_(3)
    character*3 :: pg
    real*8 :: xd(3), xdo(3), lrotm(3,3,48), ravg(3,3)

    real*8, parameter :: tighteps = 1d-7

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! initialize
    x_ = x

    ! move the atom
    if (type == atlisttype_nneq) then
       if (forcewyc) then
          ! calculate the symmetrizing matrix
          pg = sys(isys)%c%sitesymm(sys(isys)%c%at(id)%x,tighteps,leqv,lrotm)
          ravg = 0d0
          do i = 1, leqv
             ravg = ravg + lrotm(:,:,i)
          end do
          ravg = ravg / leqv

          ! calculate and project the displacement vector (fractional) onto the symmetry element
          xd = x - sys(isys)%c%at(id)%x
          xd = xd - nint(xd)

          ! identify which coordinate has changed
          ichange = -1
          xdo = xd
          if (count(abs(xdo) < tighteps) == 2) then
             do i = 1, 3
                if (abs(xdo(i)) > tighteps) then
                   ichange = i
                   exit
                end if
             end do
          end if

          ! calculate the transformed difference
          xd = matmul(ravg,xd)

          ! if no displacements, do nothing
          if (norm2(xd) < tighteps) return

          ! scale to have the correct number on the relevant coordinate
          if (ichange > 0 .and. abs(xd(ichange)) > tighteps) &
             xd = xd / xd(ichange) * xdo(ichange)

          ! apply the displacement
          x_ = sys(isys)%c%at(id)%x + xd
       end if

       ! displace
       call sys(isys)%c%move_atom(id,x_,iunit_fractional,.true.,.false.)
    elseif (type == atlisttype_ncel_frac) then
       call sys(isys)%c%move_atom(id,x_,iunit_fractional,.false.,.false.)
    elseif (type == atlisttype_ncel_bohr) then
       if (sys(isys)%c%ismolecule) &
          x_ = x - sys(isys)%c%molx0
       call sys(isys)%c%move_atom(id,x_,iunit_bohr,.false.,.false.)
    elseif (type == atlisttype_ncel_ang) then
       if (sys(isys)%c%ismolecule) then
          x_ = x/bohrtoa - sys(isys)%c%molx0
       else
          x_ = x/bohrtoa
       end if
       call sys(isys)%c%move_atom(id,x_,iunit_bohr,.false.,.false.)
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_atom_position

  ! For the atom identifier id corresponding to the given atom type,
  ! set the atomic position(s) in the system.
  module subroutine reread_geometry_from_file(sysc)
    use crystalseedmod, only: crystalseed, read_seeds_from_file, realloc_crystalseed
    use param, only: isformat_r_from_library
    class(sysconf), intent(inout) :: sysc

    integer :: isys
    character(len=:), allocatable :: file, errmsg
    logical :: exist
    integer :: isformat, mol, nseed
    type(crystalseed), allocatable :: seed(:)
    logical :: collapse

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! make sure the file exists
    if (sysc%seed%ismolecule) then
       mol = 1
    else
       mol = 0
    end if
    isformat = sysc%seed%isformat
    file = sysc%seed%file
    if (isformat /= isformat_r_from_library) then
       inquire(file=file,exist=exist)
       if (.not.exist) return
    end if

    ! re-read the seeds from file
    errmsg = ""
    nseed = 0
    allocate(seed(1))
    if (isformat == isformat_r_from_library) then
       call seed(1)%read_library(file,exist)
       if (.not.exist) return
    else
       call read_seeds_from_file(file,mol,isformat,.false.,nseed,seed,collapse,errmsg)
       if (len_trim(errmsg) > 0 .or. nseed == 0) return

       ! move the relevant seed to the first position
       if (nseed > 1) then
          seed(1) = seed(sysc%idseed)
          call realloc_crystalseed(seed,1)
       end if
    end if

    ! reset the geometry from the seed
    call sys(isys)%c%struct_new(seed(1),crashfail=.true.)

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine reread_geometry_from_file

  !xx! private procedures

  ! Thread worker: run over all systems and initialize the ones that are not locked
  function initialization_thread_worker(arg)
    use gui_main, only: force_quit_threads
    use interfaces_threads, only: thrd_success, mtx_unlock, mtx_trylock
    use tools_io, only: string, uout
    use param, only: isformat_r_crystal, isformat_r_fchk, isformat_r_gaussian, ivformat_crystal_out,&
       ivformat_gaussian_log, ivformat_gaussian_fchk, ivformat_unknown
    type(c_ptr), value :: arg
    integer(c_int) :: initialization_thread_worker
    type(thread_info), pointer :: ti

    integer :: i, i0, i1
    integer(c_int) :: idum
    integer :: iff, ivformat
    character(len=:), allocatable :: errmsg

    ! recover the thread info pointer
    call c_f_pointer(arg,ti)
    ti%active = .true.

    ! run over systems
    i0 = 1
    i1 = nsys
    do i = i0, i1
       ! check if they want us to quit
       if (force_quit_threads) exit
       ! try to grab the lock for this system
#ifdef _THREADS
       if (c_associated(sysc(i)%thread_lock)) then
          idum = mtx_trylock(sysc(i)%thread_lock)
          if (idum == thrd_success) then
#endif
             ! see if we can load it - only uninitialized and not hidden
             if (sysc(i)%status == sys_loaded_not_init.and..not.sysc(i)%hidden) then
                sysc(i)%status = sys_initializing
                ! load the seed
                call sys(i)%new_from_seed(sysc(i)%seed,ti=ti)

                ! load any fields
                if (sysc(i)%has_field) then
                   call sys(i)%load_field_string(sysc(i)%seed%file,.false.,iff,errmsg,ti=ti)
                   if (len_trim(errmsg) > 0) then
                      write (uout,'("!! Warning !! Could not read field for system: ",A)') string(i)
                      write (uout,'("!! Warning !! Error message: ",A)') trim(errmsg)
                   end if
                end if

                ! load vibration data
                if (sysc(i)%has_vib) then
                   if (sysc(i)%seed%isformat == isformat_r_crystal) then
                      ivformat = ivformat_crystal_out
                   elseif (sysc(i)%seed%isformat == isformat_r_gaussian) then
                      ivformat = ivformat_gaussian_log
                   elseif (sysc(i)%seed%isformat == isformat_r_fchk) then
                      ivformat = ivformat_gaussian_fchk
                   else
                      ivformat = ivformat_unknown
                   end if
                   call sys(i)%c%vib%read_file(sys(i)%c,sysc(i)%seed%file,"",ivformat,errmsg,ti=ti)
                   if (len_trim(errmsg) > 0) then
                      write (uout,'("!! Warning !! Could not read vibration data for system: ",A)') string(i)
                      write (uout,'("!! Warning !! Error message: ",A)') trim(errmsg)
                   end if
                end if

                ! all data is in-place, so we flag the system as ready
                sysc(i)%status = sys_ready

                ! initialize the scene and build the initial draw list
                call sysc(i)%sc%init(i)

                ! this system has been initialized
                call sysc(i)%post_event(lastchange_geometry,keepfields=.true.)
                sysc(i)%status = sys_init
             end if

#ifdef _THREADS
             ! unlock
             idum = mtx_unlock(sysc(i)%thread_lock)
          end if
       end if
#endif
    end do
    initialization_thread_worker = 0
    ti%active = .false.

  end function initialization_thread_worker

end submodule proc
