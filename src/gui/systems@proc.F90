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

  ! Re-write the seed names from the full-path names. Each system is
  ! shown with only its file name (last path component) unless that would
  ! collide with another system's name; in that case, the minimum number
  ! of parent directories needed to tell the clashing systems apart is
  ! prepended. Renamed systems keep their user-given name and are ignored.
  module subroutine system_shorten_names()
    use tools_io, only: string
    use param, only: dirsep

    integer :: i, a, b, np, d, ndup, si
    integer, allocatable :: plist(:)
    logical :: unique
    character(len=:), allocatable :: fni, fnj

    ! collect the participating systems (non-empty and not renamed)
    allocate(plist(nsys))
    np = 0
    do i = 1, nsys
       if (sysc(i)%status /= sys_empty .and..not.sysc(i)%renamed) then
          np = np + 1
          plist(np) = i
       end if
    end do

    ! for each system, use the shortest run of trailing path components
    ! that is unique among the other (non-duplicate) participating systems
    do a = 1, np
       i = plist(a)
       fni = trim(sysc(i)%fullname)
       d = 0
       do
          d = d + 1
          si = laststart(fni,d) ! start of this system's last d components
          unique = .true.
          do b = 1, np
             if (b == a) cycle
             fnj = trim(sysc(plist(b))%fullname)
             if (fnj == fni) cycle ! an exact duplicate, handled below
             ! not unique if the last d path components coincide
             if (fni(si:) == fnj(laststart(fnj,d):)) then
                unique = .false.
                exit
             end if
          end do
          if (unique) then
             sysc(i)%seed%name = fni(si:)
             exit
          elseif (si == 1) then
             ! the whole path is already used (should only happen for exact
             ! duplicates); stop here and let the duplicate pass disambiguate
             sysc(i)%seed%name = fni
             exit
          end if
       end do
    end do

    ! exact duplicates (same full-path name, e.g. from Duplicate or reload)
    ! cannot be told apart by the path: append a counter
    do a = 1, np
       i = plist(a)
       ndup = 0
       do b = 1, a-1
          if (trim(sysc(plist(b))%fullname) == trim(sysc(i)%fullname)) ndup = ndup + 1
       end do
       if (ndup > 0) &
          sysc(i)%seed%name = trim(sysc(i)%seed%name) // " (" // string(ndup+1) // ")"
    end do

  contains
    !> Starting character index of the last d path components of str
    !> (the part after the d-th dirsep counted from the end; 1 if str has
    !> fewer than d components).
    integer function laststart(str,d)
      character(len=*), intent(in) :: str
      integer, intent(in) :: d
      integer :: p, c
      laststart = 1
      c = 0
      do p = len_trim(str), 1, -1
         if (str(p:p) == dirsep) then
            c = c + 1
            if (c == d) then
               laststart = p + 1
               return
            end if
         end if
      end do
    end function laststart
  end subroutine system_shorten_names

  !> Add systems from the given seeds. If collapse is present and
  !> true, reorder to make the last system be first and collapse the
  !> other seeds in the tree view. If iafield, load a field from that
  !> seed. If iavib, load vibrational data from that seed. If forceidx
  !> is present, force the system to be in the provided index.
  module subroutine add_systems_from_seeds(nseed,seed,collapse,iafield,iavib,forceidx,idlist)
    use utils, only: get_current_working_dir
    use grid1mod, only: grid1_register_ae
    use gui_main, only: reuse_mid_empty_systems
    use windows, only: regenerate_window_pointers, win, iwin_tree
    use interfaces_threads, only: allocate_mtx, mtx_init, mtx_plain
    use crystalseedmod, only: crystalseed
    use global, only: symprec
    use tools_io, only: uout
    use types, only: realloc
    use param, only: dirsep, atmcov0
    integer, intent(in) :: nseed
    type(crystalseed), allocatable, intent(in) :: seed(:)
    logical, intent(in), optional :: collapse
    integer, intent(in), optional :: iafield, iavib
    integer, intent(in), optional :: forceidx
    integer, allocatable, intent(out), optional :: idlist(:)

    integer :: i, j, nid, idum
    integer :: iafield_, iavib_, forceidx_
    integer :: iseed, iseed_, idx
    character(len=:), allocatable :: errmsg, str
    type(system), allocatable :: syaux(:)
    type(sysconf), allocatable :: syscaux(:)
    logical :: collapse_, isrun, isabspath
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
       isabspath = (str(1:1) == dirsep .or. str(1:1) == "/")
#ifdef _WIN32
       if (.not.isabspath .and. len(str) >= 2) then
          isabspath = (str(2:2) == ":")
       end if
#endif
       if (isabspath) then
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

       ! set the covalent radii
       sysc(idx)%atmcov = atmcov0

       ! set the symmetry epsilon
       sysc(idx)%symeps = symprec

       ! redo everything
       call sysc(idx)%post_event(lastchange_geometry)
    end do

    ! select the first new system in the tree
    if (allocated(win)) then
       if (iwin_tree > 0 .and. iwin_tree <= size(win)) &
          win(iwin_tree)%forceselect = id(1)
    end if

    ! return the IDs of the new systems if requested
    if (present(idlist)) then
       idlist = id
    end if
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
    sysc(idx)%md_run = .false.
    call sysc(idx)%md%free()
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
    use param, only: isformat_r_from_library, isformat_r_derived
    integer, intent(in) :: idx

    character(len=:), allocatable :: file, errmsg
    logical :: exist
    integer :: isformat, mol, iafield, iavib, nseed, icol, idseed
    type(crystalseed), allocatable :: seed(:)
    logical :: collapse, ihid, renamed, derived
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

    ! derived (in-memory) systems are not backed by a file: keep a copy of
    ! the initial seed and rebuild the system from it
    derived = (isformat == isformat_r_derived)
    iafield = 0
    iavib = 0
    allocate(seed(1))
    if (derived) then
       seed(1) = sysc(idx)%seed
    else
       ! make sure the file exists
       file = sysc(idx)%seed%file
       if (isformat /= isformat_r_from_library) then
          inquire(file=file,exist=exist)
          if (.not.exist) return
       end if
    end if

    ! terminate the system
    call sys(idx)%end()
    call sysc(idx)%seed%end()
    sysc(idx)%status = sys_empty

    ! re-read the seeds from file (derived systems already have the seed)
    if (.not.derived) then
       errmsg = ""
       nseed = 0
       if (isformat == isformat_r_from_library) then
          call seed(1)%read_library(file,exist)
          if (.not.exist) return
       else
          call read_seeds_from_file(file,mol,isformat,.false.,nseed,seed,collapse,errmsg,&
             iafield,iavib)
          if (len_trim(errmsg) > 0 .or. nseed == 0) return

          ! move the relevant seed to the first position
          if (nseed > 1) then
             seed(1) = seed(idseed)
             call realloc_crystalseed(seed,1)
          end if
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

  !> Return .true. if the (optional) errmsg is present and non-empty.
  function has_errmsg(errmsg)
    character(len=:), allocatable, intent(in), optional :: errmsg
    logical :: has_errmsg

    has_errmsg = .false.
    if (present(errmsg)) has_errmsg = (len_trim(errmsg) > 0)

  end function has_errmsg

  !> Set the time for last change at level level. If keepfields is
  !> present and true, do not reset the associated fields.
  module subroutine post_event(sysc,level,keepfields,nocapture)
    use interfaces_glfw, only: glfwGetTime
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: level
    logical, intent(in), optional :: keepfields
    logical, intent(in), optional :: nocapture

    real*8 :: time
    logical :: keepfields_, nocapture_

    keepfields_ = .false.
    if (present(keepfields)) keepfields_ = keepfields
    nocapture_ = .false.
    if (present(nocapture)) nocapture_ = nocapture

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
       ! record in the undo/redo history, unless the caller is restoring a
       ! previous state (which also posts a geometry event)
       if (.not.nocapture_ .and. ok_system(sysc%id,sys_init)) &
          call sysc%undo_capture(time)
    end if

  end subroutine post_event

  !> Recompute the bonds/connectivity for this system using its configured
  !> bonding parameters (per-species covalent radii, bond factor, metal-bond
  !> tolerance), then signal the scene/representations to refresh. This is the
  !> single entry point for the "Recalculate bonds" action (Tools menu, Ctrl+B)
  !> and the View/Edit Geometry "Apply" button.
  module subroutine rebond(sysc)
    class(sysconf), intent(inout) :: sysc

    call sys(sysc%id)%c%rebond(sysc%atmcov,sysc%bondfactor,bonddelta=sysc%bonddelta,&
       allowed=sysc%bondallowed)
    call sysc%post_event(lastchange_rebond)

  end subroutine rebond

  !> Reset the undo/redo history for this system: discard all saved states
  !> and re-seed the history with the current geometry (if the system is
  !> initialized). Called when a system is first initialized or reloaded.
  module subroutine undo_reset(sysc)
    class(sysconf), intent(inout) :: sysc

    sysc%undo_n = 0
    sysc%undo_icur = 0
    sysc%undo_ibase = 1
    ! seed the history with the current geometry (undo_capture is a no-op if
    ! the system is not yet initialized); the time is irrelevant here
    call sysc%undo_capture(0d0)
    ! the next edit must start a new entry, not coalesce with the seed
    sysc%undo_lasttime = -1d30

  end subroutine undo_reset

  !> Capture the current geometry as a new state in the undo/redo history.
  !> time is the current GLFW time, used to coalesce rapid consecutive
  !> captures (a continuous drag produces one capture per frame) into a
  !> single history entry. Discards any redo states ahead of the current.
  module subroutine undo_capture(sysc,time)
    class(sysconf), intent(inout) :: sysc
    real*8, intent(in) :: time

    integer :: isys
    logical :: copysym, coalesce

    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (.not.allocated(sysc%undo_seed)) allocate(sysc%undo_seed(undo_maxdepth))

    ! drop any redo states ahead of the current one
    sysc%undo_n = sysc%undo_icur

    ! coalesce consecutive captures (frames of an interactive drag) by
    ! overwriting the current top state instead of appending a new one
    coalesce = (sysc%undo_icur >= 1) .and. ((time - sysc%undo_lasttime) < undo_coalesce_time)
    if (.not.coalesce) then
       ! drop the oldest state if the history is full by advancing the base of
       ! the ring buffer; the new state then reuses the freed (oldest) slot
       if (sysc%undo_n >= undo_maxdepth) then
          sysc%undo_ibase = modulo(sysc%undo_ibase,undo_maxdepth) + 1
          sysc%undo_n = undo_maxdepth - 1
       end if
       sysc%undo_n = sysc%undo_n + 1
       sysc%undo_icur = sysc%undo_n
    end if

    ! save the current geometry into the (possibly new) top slot, inheriting
    ! both the symmetry (when available) and the bonding, so that restoring
    ! recovers the exact same system
    copysym = (.not.sys(isys)%c%ismolecule .and. sys(isys)%c%spgavail)
    call sys(isys)%c%makeseed(sysc%undo_seed(undo_slot(sysc,sysc%undo_icur)),copysym=copysym,&
       copybonding=.true.)
    sysc%undo_lasttime = time

  end subroutine undo_capture

  !> Restore the previous geometry state from the undo history.
  module subroutine undo(sysc)
    class(sysconf), intent(inout) :: sysc

    if (.not.sysc%can_undo()) return
    sysc%undo_icur = sysc%undo_icur - 1
    call undo_restore(sysc)

  end subroutine undo

  !> Restore the next geometry state from the redo history.
  module subroutine redo(sysc)
    class(sysconf), intent(inout) :: sysc

    if (.not.sysc%can_redo()) return
    sysc%undo_icur = sysc%undo_icur + 1
    call undo_restore(sysc)

  end subroutine redo

  !> .true. if there is a previous state to undo to.
  module function can_undo(sysc)
    class(sysconf), intent(in) :: sysc
    logical :: can_undo

    can_undo = (sysc%undo_icur > 1)

  end function can_undo

  !> .true. if there is a next state to redo to.
  module function can_redo(sysc)
    class(sysconf), intent(in) :: sysc
    logical :: can_redo

    can_redo = (sysc%undo_icur < sysc%undo_n)

  end function can_redo

  !> Rebuild the system's crystal structure from the current state in the
  !> undo history and refresh the scene. Helper for undo/redo; the restore
  !> posts with nocapture so it is not itself recorded as a new history entry.
  subroutine undo_restore(sysc)
    class(sysconf), intent(inout) :: sysc

    integer :: isys

    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sysc%undo_icur < 1 .or. sysc%undo_icur > sysc%undo_n) return

    call sys(isys)%c%struct_new(sysc%undo_seed(undo_slot(sysc,sysc%undo_icur)),crashfail=.true.)
    sysc%sc%nextbuildlists_fixcam = .true.
    call sysc%post_event(lastchange_geometry,nocapture=.true.)

    ! force the next edit to start a new history entry instead of coalescing
    ! with (overwriting) the state we just restored
    sysc%undo_lasttime = -1d30

  end subroutine undo_restore

  !> Map a logical history index j (1:undo_n) to its physical slot in the
  !> undo_seed ring buffer, accounting for the moving base pointer.
  pure function undo_slot(sysc,j)
    class(sysconf), intent(in) :: sysc
    integer, intent(in) :: j
    integer :: undo_slot

    undo_slot = modulo(sysc%undo_ibase - 1 + (j - 1),undo_maxdepth) + 1

  end function undo_slot

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
  !> present, clear atoms with indices idx of the given atlisttype.
  module subroutine highlight_clear(sysc,transient,idx,type)
    class(sysconf), intent(inout) :: sysc
    logical, intent(in) :: transient
    integer, intent(in), optional :: idx(:)
    integer, intent(in), optional :: type

    integer :: nat, iat
    logical :: changed
    real(c_float), allocatable :: highlight_aux(:,:)
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

    ! clear highlights
    if (transient) then
       changed = any(sysc%highlight_rgba_transient >= 0._c_float)
       sysc%highlight_rgba_transient = -1._c_float
    else
       if (present(idx) .and. present(type)) then
          ! save a copy of the current highlight
          allocate(highlight_aux(4,nat))
          highlight_aux = sysc%highlight_rgba

          ! clear the highlight for the atoms with the given ids/type
          call sysc%attype_celatom_mask(type,idx,imask=imask)
          do iat = 1, nat
             if (imask(iat) /= 0) sysc%highlight_rgba(:,iat) = -1._c_float
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

  !> Select (persistent highlight) all cell atoms in the system.
  module subroutine highlight_all(sysc)
    use gui_main, only: ColorHighlightSelectScene
    class(sysconf), intent(inout) :: sysc

    integer :: i, nat
    integer, allocatable :: iat(:)

    if (sysc%status < sys_init) return
    nat = sys(sysc%id)%c%ncel
    if (nat <= 0) return

    allocate(iat(nat))
    do i = 1, nat
       iat(i) = i
    end do
    call sysc%highlight_atoms(.false.,iat,atlisttype_ncel_frac,&
       spread(ColorHighlightSelectScene,2,nat))

  end subroutine highlight_all

  !> Invert the persistent selection of all cell atoms in the system. Atoms
  !> currently selected become unselected and vice versa.
  module subroutine highlight_invert(sysc)
    use gui_main, only: ColorHighlightSelectScene
    class(sysconf), intent(inout) :: sysc

    integer :: i, nat, nsel, nuns
    integer, allocatable :: sel(:), uns(:)

    if (sysc%status < sys_init) return
    nat = sys(sysc%id)%c%ncel
    if (nat <= 0) return

    ! split the cell atoms into selected and unselected
    call sysc%highlighted_atom_list(nsel,sel)
    allocate(uns(nat))
    nuns = 0
    do i = 1, nat
       if (.not.any(sel(1:nsel) == i)) then
          nuns = nuns + 1
          uns(nuns) = i
       end if
    end do

    if (nsel > 0) &
       call sysc%highlight_clear(.false.,sel(1:nsel),atlisttype_ncel_frac)
    if (nuns > 0) &
       call sysc%highlight_atoms(.false.,uns(1:nuns),atlisttype_ncel_frac,&
          spread(ColorHighlightSelectScene,2,nuns))

  end subroutine highlight_invert

  !> Return the list of highlighted (selected) cell atoms. On output, nat is the
  !> number of selected atoms and iat(1:nat) their cell-atom indices. iat is
  !> allocated to ncel; only the first nat entries are meaningful. Returns
  !> nat = 0 if the system is not initialized or has no selection.
  module subroutine highlighted_atom_list(sysc,nat,iat)
    class(sysconf), intent(in) :: sysc
    integer, intent(out) :: nat
    integer, allocatable, intent(inout) :: iat(:)

    integer :: i, ncel

    nat = 0
    if (sysc%status < sys_init) return
    ncel = sys(sysc%id)%c%ncel
    if (ncel <= 0) return

    if (allocated(iat)) then
       if (size(iat,1) < ncel) deallocate(iat)
    end if
    if (.not.allocated(iat)) allocate(iat(ncel))

    if (.not.allocated(sysc%highlight_rgba)) return
    do i = 1, min(ncel,size(sysc%highlight_rgba,2))
       if (any(sysc%highlight_rgba(:,i) >= 0._c_float)) then
          nat = nat + 1
          iat(nat) = i
       end if
    end do

  end subroutine highlighted_atom_list

  !> Create a new system containing only the highlighted (selected)
  !> atoms of this system. By default the new system preserves the
  !> parent type (crystal or molecule). If forcemolecule is present
  !> and true, the selection is always turned into a non-periodic
  !> molecule.
  module subroutine new_system_from_highlighted(sysc,forcemolecule)
    use crystalseedmod, only: crystalseed
    use types, only: realloc
    use global, only: rborder_def
    use param, only: isformat_r_derived
    class(sysconf), intent(inout) :: sysc
    logical, intent(in), optional :: forcemolecule

    integer :: i, nat, id
    integer, allocatable :: iat(:)
    type(crystalseed), allocatable :: seed(:)
    logical :: molecule

    ! consistency checks
    id = sysc%id
    if (.not.ok_system(id,sys_init)) return

    ! list the selected cell atoms; no-op if nothing is selected
    call sysc%highlighted_atom_list(nat,iat)
    if (nat == 0) return

    ! whether to force the new system to be a molecule
    molecule = .false.
    if (present(forcemolecule)) molecule = forcemolecule

    ! build a seed from the crystal, then keep only the selected atoms
    allocate(seed(1))
    call sys(id)%c%makeseed(seed(1),copysym=.false.)
    call realloc(seed(1)%x,3,nat)
    call realloc(seed(1)%is,nat)
    call realloc(seed(1)%atname,nat)
    do i = 1, nat
       if (molecule) then
          ! absolute Cartesian (bohr) coordinates for a molecule seed
          seed(1)%x(:,i) = sys(id)%c%atcel(iat(i))%r + sys(id)%c%molx0
       else
          seed(1)%x(:,i) = sys(id)%c%atcel(iat(i))%x
       end if
       seed(1)%is(i) = sys(id)%c%atcel(iat(i))%is
       seed(1)%atname(i) = sys(id)%c%at(sys(id)%c%atcel(iat(i))%idx)%name
    end do
    seed(1)%nat = nat

    ! turn the seed into a non-periodic molecule if requested
    if (molecule) then
       seed(1)%ismolecule = .true.
       seed(1)%useabr = 0
       seed(1)%havesym = 0
       seed(1)%neqv = 0
       seed(1)%ncv = 0
       seed(1)%havex0 = .false.
       seed(1)%molx0 = 0d0
       seed(1)%border = rborder_def
       seed(1)%cubic = .false.
    end if
    seed(1)%name = trim(sysc%seed%name) // " (selection)"

    ! derived system
    seed(1)%file = ""
    seed(1)%isformat = isformat_r_derived

    ! create the new system (add_systems_from_seeds selects it in the tree)
    call add_systems_from_seeds(1,seed)
    call launch_initialization_thread()

  end subroutine new_system_from_highlighted

  !> Remove, merge or duplicate the highlighted atoms in the system.
  module subroutine edit_highlighted_atoms(sysc,remove,merge,duplicate,errmsg)
    class(sysconf), intent(inout) :: sysc
    logical, intent(in), optional :: remove, merge, duplicate
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: nat, id
    integer, allocatable :: iat(:)

    errmsg = ""

    ! consistency checks
    id = sysc%id
    if (.not.ok_system(id,sys_init)) return

    ! nothing to do if no atom is highlighted (the highlight array may not even
    ! be allocated, e.g. when the Delete key is pressed with only a
    ! transient/symmetry highlight)
    call sysc%highlighted_atom_list(nat,iat)
    if (nat == 0) return

    ! remove/merge/duplicate the atoms
    call sys(id)%c%edit_atom_list(nat,iat(1:nat),remove,merge,duplicate,errmsg=errmsg)
    if (has_errmsg(errmsg)) return

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine edit_highlighted_atoms

  !> Remove, merge or duplicate the highlighted species in the system.
  module subroutine edit_highlighted_species(sysc,selected,remove,merge,duplicate,errmsg)
    use crystalseedmod, only: crystalseed
    use types, only: species
    class(sysconf), intent(inout) :: sysc
    logical, intent(in) :: selected(:)
    logical, intent(in), optional :: remove, merge, duplicate
    character(len=:), allocatable, intent(inout) :: errmsg

    type(crystalseed) :: seed
    integer :: i, nat, id, ipres, nspc
    integer, allocatable :: iat(:), map(:)
    logical :: remove_, merge_, duplicate_
    type(species), allocatable :: spc(:)

    errmsg = ""

    ! check input options
    remove_ = .false.
    merge_ = .false.
    duplicate_ = .false.
    ipres = 0
    if (present(remove)) then
       remove_ = remove
       if (remove_) ipres = ipres + 1
    end if
    if (present(merge)) then
       merge_ = merge
       if (merge_) ipres = ipres + 1
    end if
    if (present(duplicate)) then
       duplicate_ = duplicate
       if (duplicate_) ipres = ipres + 1
    end if
    if (ipres == 0) return
    if (ipres > 1) then
       errmsg = 'more than one of merge/remove/duplicate'
       return
    end if

    ! consistency checks
    id = sysc%id
    if (.not.ok_system(id,sys_init)) return
    if (size(selected,1) /= sys(id)%c%nspc) return

    ! remove/merge the highlighted atoms
    if (remove_ .or. merge_) then
       call sysc%highlighted_atom_list(nat,iat)

       ! remove/merge/duplicate the atoms
       if (nat > 0) then
          call sys(id)%c%edit_atom_list(nat,iat(1:nat),remove,merge,duplicate,errmsg=errmsg)
          if (has_errmsg(errmsg)) return
       end if
    end if

    if (remove_ .or. duplicate_) then
       ! make seed from this crystal
       call sys(id)%c%makeseed(seed,copysym=.false.)

       if (remove_) then
          ! remove the species
          allocate(spc(count(.not.selected)),map(sys(id)%c%nspc))
          map = 0
          nspc = 0
          do i = 1, sys(id)%c%nspc
             if (.not.selected(i)) then
                nspc = nspc + 1
                spc(nspc) = sys(id)%c%spc(i)
                map(i) = nspc
             end if
          end do
       elseif (duplicate_) then
          ! duplicate the species
          allocate(spc(sys(id)%c%nspc+count(selected)),map(sys(id)%c%nspc))
          map = 0
          nspc = 0
          do i = 1, sys(id)%c%nspc
             nspc = nspc + 1
             spc(nspc) = sys(id)%c%spc(i)
             map(i) = nspc
             if (selected(i)) then
                nspc = nspc + 1
                spc(nspc) = sys(id)%c%spc(i)
             end if
          end do
       end if

       ! remap the species in the seed
       seed%nspc = nspc
       seed%spc = spc
       do i = 1, seed%nat
          seed%is(i) = map(seed%is(i))
       end do

       ! build the new crystal
       call sys(id)%c%struct_new(seed,crashfail=.false.)
       if (.not.sys(id)%c%isinit) then
          errmsg = "Could not rebuild the structure after editing the species"
          return
       end if
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine edit_highlighted_species

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
    if (.not.sys(id)%c%ismolecule .and. allowedin(atlisttype_ncel_frac)) &
       strcombo = strcombo // "Fractional" // c_null_char
    if (allowedin(atlisttype_ncel_bohr)) then
       if (units_) then
          strcombo = strcombo // "Cartesian (bohr)" // c_null_char
       else
          strcombo = strcombo // "Cartesian" // c_null_char
       end if
    end if
    if (allowedin(atlisttype_ncel_ang)) then
       if (units_) then
          strcombo = strcombo // "Cartesian (Å)" // c_null_char
       else
          strcombo = strcombo // "Cartesian" // c_null_char
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
    call iw_combo_simple(label,strcombo,icombo,changed=attype_combo_simple)

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

  ! For the given atom type, set the corresponding atomic species.
  module subroutine set_attype_species(sysc,type,id,is,copybonding)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    integer, intent(in) :: is
    logical, intent(in), optional :: copybonding

    integer :: isys, i, nat
    integer, allocatable :: iat(:)
    logical :: ok, copybonding_

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

    ! initialize
    allocate(iat(sys(isys)%c%ncel))
    nat = 0

    ! pick the atoms
    do i = 1, sys(isys)%c%ncel
       if (type == atlisttype_nneq) then
          ok = (sys(isys)%c%atcel(i)%idx == id)
       elseif (type == atlisttype_ncel_frac .or. type == atlisttype_ncel_bohr .or. type == atlisttype_ncel_ang) then
          ok = (i == id)
       end if
       if (ok) then
          nat = nat + 1
          iat(nat) = i
       end if
    end do

    ! execute
    call sys(isys)%c%change_atom_species(nat,iat,is,copybonding=copybonding_)

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_attype_species

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

  ! For the given atom type and id, set the site occupancy to occ. The
  ! occupancy lives on the non-equivalent atom, so a change in the cell list
  ! applies to all symmetry-equivalent atoms.
  module subroutine set_attype_occupancy(sysc,type,id,occ)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    real*8, intent(in) :: occ

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! occupancy is clamped to (0,1] (0 would be a vacancy, rejected by seed_check)
    if (type == atlisttype_nneq) then
       sys(isys)%c%at(id)%occ = min(max(occ,1d-10),1d0)
    elseif (type /= atlisttype_species .and. type /= atlisttype_nmol) then
       sys(isys)%c%at(sys(isys)%c%atcel(id)%idx)%occ = min(max(occ,1d-10),1d0)
    else
       return
    end if

    ! update the partial-occupancy flag
    call sys(isys)%c%set_haveocc()

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_attype_occupancy

  ! For the given atom type and id, return the site occupancy. The occupancy
  ! lives on the non-equivalent atom, reached directly (nneq) or via the
  ! complete-list index (cell lists).
  module function attype_occupancy(sysc,type,id)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    real*8 :: attype_occupancy

    integer :: isys

    attype_occupancy = 1d0
    isys = sysc%id
    if (type == atlisttype_nneq) then
       attype_occupancy = sys(isys)%c%at(id)%occ
    elseif (type /= atlisttype_species .and. type /= atlisttype_nmol) then
       attype_occupancy = sys(isys)%c%at(sys(isys)%c%atcel(id)%idx)%occ
    end if

  end function attype_occupancy

  ! For the given atom type and id, return the list of occupants of a mixed
  ! site ("Ta 0.500 / Mg 0.250 / Na 0.250"), or an empty string if not mixed.
  module function attype_mixed(sysc,type,id)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    character(len=:), allocatable :: attype_mixed

    integer :: isys, ineq

    attype_mixed = ""
    isys = sysc%id
    if (type == atlisttype_nneq) then
       ineq = id
    elseif (type /= atlisttype_species .and. type /= atlisttype_nmol) then
       ineq = sys(isys)%c%atcel(id)%idx
    else
       return
    end if
    attype_mixed = sys(isys)%c%mix_string(ineq,decimal=3)

  end function attype_mixed

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

  ! For the given atom type, return the units of the atomic coordinates
  ! ("Å", "bohr", or "fractional"), for use in column headers and labels.
  module function attype_coordinates_units(sysc,type)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    character(len=:), allocatable :: attype_coordinates_units

    if (type == atlisttype_ncel_ang) then
       attype_coordinates_units = "Å"
    elseif (type == atlisttype_ncel_bohr) then
       attype_coordinates_units = "bohr"
    else
       attype_coordinates_units = "fractional"
    end if

  end function attype_coordinates_units

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

  !> Add atom of the given atom type with species is and coordinates
  !> x. If is <= 0, add the species with Z = abs(is) to the system.
  module subroutine attype_add_atom(sysc,type,is,x)
    use global, only: iunit_bohr, iunit_fractional
    use param, only: bohrtoa
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: is
    real*8, intent(in) :: x(3)

    integer :: isys
    real*8 :: x_(3)

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! initialize
    x_ = x

    ! add the atom
    if (type == atlisttype_nneq) then
       call sys(isys)%c%add_atom(is,x_,iunit_fractional,.true.)
    elseif (type == atlisttype_ncel_frac) then
       call sys(isys)%c%add_atom(is,x_,iunit_fractional,.false.)
    elseif (type == atlisttype_ncel_bohr) then
       if (sys(isys)%c%ismolecule) x_ = x - sys(isys)%c%molx0
       call sys(isys)%c%add_atom(is,x_,iunit_bohr,.false.)
    elseif (type == atlisttype_ncel_ang) then
       if (sys(isys)%c%ismolecule) then
          x_ = x/bohrtoa - sys(isys)%c%molx0
       else
          x_ = x/bohrtoa
       end if
       call sys(isys)%c%add_atom(is,x_,iunit_bohr,.false.)
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine attype_add_atom

  !> Reorder the atoms in the system for the given atom types using
  !> the permutation iord.
  module subroutine attype_reorder(sysc,type,iord)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: iord(:)

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! reorder
    if (type == atlisttype_nneq) then
       call sys(isys)%c%reorder_atoms(iord,.true.)
    elseif (type == atlisttype_species) then
       call sys(isys)%c%reorder_species(iord)
    elseif (type == atlisttype_nmol) then
       call sys(isys)%c%reorder_molecules(iord)
    else
       call sys(isys)%c%reorder_atoms(iord,.false.)
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine attype_reorder

  !> Swap two of the atoms in the system with IDs i1 and i2 for the
  !> given atom type.
  module subroutine attype_swap_atoms(sysc,type,i1,i2)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: i1, i2

    integer :: isys, ntype, i
    integer, allocatable :: iord(:)

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! get the number of entities
    ntype = sysc%attype_number(type)
    if (i1 < 1 .or. i1 > ntype .or. i2 < 1 .or. i2 > ntype) return

    ! build the permutation
    allocate(iord(ntype))
    do i = 1, ntype
       iord(i) = i
    end do
    iord(i2) = i1
    iord(i1) = i2

    ! reorder
    call sysc%attype_reorder(type,iord)
    deallocate(iord)

  end subroutine attype_swap_atoms

  !> Swap the two molecular fragments in the system with IDs i1 and i2.
  module subroutine swap_molecules(sysc,i1,i2)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: i1, i2

    integer :: isys, nmol, i
    integer, allocatable :: iperm(:)

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! check the molecule IDs are in range
    nmol = sys(isys)%c%nmol
    if (i1 < 1 .or. i1 > nmol .or. i2 < 1 .or. i2 > nmol) return

    ! build the permutation that swaps the two molecules
    allocate(iperm(nmol))
    do i = 1, nmol
       iperm(i) = i
    end do
    iperm(i2) = i1
    iperm(i1) = i2

    ! reorder the molecules and signal the geometry change
    call sys(isys)%c%reorder_molecules(iperm)
    deallocate(iperm)
    call sysc%post_event(lastchange_geometry)

  end subroutine swap_molecules

  ! For the atom identifier id corresponding to the given atom type,
  ! set the atomic position(s) in the system.
  module subroutine set_atom_position(sysc,type,id,x,forcewyc,copybonding)
    use global, only: iunit_bohr, iunit_fractional
    use param, only: bohrtoa
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    real*8, intent(in) :: x(3)
    logical, intent(in) :: forcewyc
    logical, intent(in), optional :: copybonding

    integer :: isys, leqv, i, ichange
    real*8 :: x_(3)
    character*3 :: pg
    real*8 :: xd(3), xdo(3), lrotm(3,3,48), ravg(3,3)
    logical :: copybonding_

    real*8, parameter :: tighteps = 1d-7

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! initialize
    x_ = x
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding

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
       call sys(isys)%c%move_atom(id,x_,iunit_fractional,.true.,.false.,copybonding=copybonding_)
    elseif (type == atlisttype_ncel_frac) then
       call sys(isys)%c%move_atom(id,x_,iunit_fractional,.false.,.false.,copybonding=copybonding_)
    elseif (type == atlisttype_ncel_bohr) then
       if (sys(isys)%c%ismolecule) &
          x_ = x - sys(isys)%c%molx0
       call sys(isys)%c%move_atom(id,x_,iunit_bohr,.false.,.false.,copybonding=copybonding_)
    elseif (type == atlisttype_ncel_ang) then
       if (sys(isys)%c%ismolecule) then
          x_ = x/bohrtoa - sys(isys)%c%molx0
       else
          x_ = x/bohrtoa
       end if
       call sys(isys)%c%move_atom(id,x_,iunit_bohr,.false.,.false.,copybonding=copybonding_)
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_atom_position

  ! For the molecule identifier id corresponding to the given molecule
  ! coordinate type, rigidly translate the fragment so that its center
  ! of mass is at position x.
  module subroutine set_molecule_position(sysc,type,id,x,copybonding)
    use global, only: iunit_bohr, iunit_fractional
    use param, only: bohrtoa
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    real*8, intent(in) :: x(3)
    logical, intent(in), optional :: copybonding

    integer :: isys
    real*8 :: x_(3)
    logical :: copybonding_

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! move the molecule (units interpreted in the internal frame, as for move_atom)
    x_ = x
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    if (type == atlisttype_ncel_frac) then
       call sys(isys)%c%move_molecule(id,x_,iunit_fractional,.false.,copybonding=copybonding_)
    elseif (type == atlisttype_ncel_bohr) then
       if (sys(isys)%c%ismolecule) &
          x_ = x - sys(isys)%c%molx0
       call sys(isys)%c%move_molecule(id,x_,iunit_bohr,.false.,copybonding=copybonding_)
    elseif (type == atlisttype_ncel_ang) then
       if (sys(isys)%c%ismolecule) then
          x_ = x/bohrtoa - sys(isys)%c%molx0
       else
          x_ = x/bohrtoa
       end if
       call sys(isys)%c%move_molecule(id,x_,iunit_bohr,.false.,copybonding=copybonding_)
    end if

    ! the geometry changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_molecule_position

  ! Rigidly rotate molecule id about its center of mass so that the
  ! Euler angles (ZYZ, radians) of its standard orientation become euler.
  module subroutine set_molecule_rotation(sysc,id,euler,copybonding)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: id
    real*8, intent(in) :: euler(3)
    logical, intent(in), optional :: copybonding

    integer :: isys
    logical :: copybonding_

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! rotate the molecule
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call sys(isys)%c%rotate_molecule(id,euler=euler,copybonding=copybonding_)

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_molecule_rotation

  ! For the atom identifier id corresponding to the given atom type,
  ! set the atomic number and the name of the corresponding species.
  module subroutine set_atomic_number(sysc,type,id,iz,setatomnames)
    use tools_io, only: nameguess
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: type
    integer, intent(in) :: id
    integer, intent(in) :: iz
    logical, intent(in), optional :: setatomnames

    integer :: isys, ispc, i
    logical :: setatomnames_

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! optional arguments
    setatomnames_ = .false.
    if (present(setatomnames)) setatomnames_ = setatomnames

    ! set the atomic number
    ispc = sysc%attype_species(type,id)
    sys(isys)%c%spc(ispc)%z = iz
    sys(isys)%c%spc(ispc)%name = nameguess(iz,.true.)

    ! reset the atom names
    if (setatomnames_) then
       do i = 1, sys(isys)%c%nneq
          if (sys(isys)%c%at(i)%is == ispc) then
             sys(isys)%c%at(i)%name = sys(isys)%c%spc(ispc)%name
          end if
       end do
    end if

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine set_atomic_number

  !> Add atomic species with atomic number iz to the system.
  module subroutine add_species(sysc,iz)
    use tools_io, only: nameguess
    use types, only: realloc
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: iz

    integer :: ispc, isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! add the species
    sys(isys)%c%nspc = sys(isys)%c%nspc + 1
    call realloc(sys(isys)%c%spc,sys(isys)%c%nspc)
    ispc = sys(isys)%c%nspc
    sys(isys)%c%spc(ispc)%z = iz
    sys(isys)%c%spc(ispc)%name = nameguess(iz,.true.)

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine add_species

  ! Remove the bond between cell atoms iat1 and iat2 (iat2 at lattice vector
  ! lvec relative to iat1) by editing the connectivity in place. Posts a rebond
  ! event (not a geometry change), so fields and selections remain valid.
  module subroutine remove_bond(sysc,iat1,iat2,lvec)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! remove the bond from the connectivity and signal a rebond
    call sys(isys)%c%remove_bond(iat1,iat2,lvec)
    call sysc%post_event(lastchange_rebond)

  end subroutine remove_bond

  ! Set the bond order (0=dashed, 1=single, 2=double, 3=triple) of the bond
  ! between cell atoms iat1 and iat2 (iat2 at lattice vector lvec relative to
  ! iat1). Posts a rebond event (not a geometry change).
  module subroutine set_bond_order(sysc,iat1,iat2,lvec,order)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer, intent(in) :: order

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! set the bond order in the connectivity and signal a rebond
    call sys(isys)%c%set_bond_order(iat1,iat2,lvec,order)
    call sysc%post_event(lastchange_rebond)

  end subroutine set_bond_order

  ! Add a bond between cell atoms iat1 and iat2 + lvec with the given
  ! bond order by editing the connectivity in place. Posts a rebond
  ! event.
  module subroutine add_bond(sysc,iat1,iat2,lvec,order)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: iat1, iat2
    integer, intent(in) :: lvec(3)
    integer, intent(in) :: order

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return

    ! add the bond to the connectivity and signal a rebond
    call sys(isys)%c%add_bond(iat1,iat2,lvec,order)
    call sysc%post_event(lastchange_rebond)

  end subroutine add_bond

  ! For the atom identifier id corresponding to the given atom type,
  ! set the atomic position(s) in the system.
  module subroutine reread_geometry_from_file(sysc)
    use crystalseedmod, only: crystalseed, read_seeds_from_file, realloc_crystalseed
    use param, only: isformat_r_from_library, isformat_r_derived
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

    ! derived systems
    if (sysc%seed%isformat == isformat_r_derived) then
       call sys(isys)%c%struct_new(sysc%seed,crashfail=.true.)
       call sysc%post_event(lastchange_geometry)
       call sysc%undo_reset()
       return
    end if

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

    ! the structure was reloaded from file: start a fresh undo history
    call sysc%undo_reset()

  end subroutine reread_geometry_from_file

  ! Change the unit cell of the system to have cell lengths aa (bohr)
  ! and angles bb (degree). If forcewyc, force the system to keep
  ! symemtry.
  module subroutine move_cell(sysc,aa,bb,forcewyc,copybonding)
    use crystalmod, only: pointgroup_info
    use tools_math, only: m_x2c_from_cellpar, cellpar_from_metric
    class(sysconf), intent(inout) :: sysc
    real*8, intent(in) :: aa(3), bb(3)
    logical, intent(in) :: forcewyc
    logical, intent(in), optional :: copybonding

    integer :: isys
    real*8 :: aa_(3), bb_(3)
    logical :: copybonding_
    integer :: leqv, i, n
    real*8 :: g(3,3), gavg(3,3), da
    real*8, allocatable :: rotm(:,:,:)

    real*8, parameter :: tighteps = 1d-7
    real*8, parameter :: epsconv = 1d-7
    integer, parameter :: maxstep = 200

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    aa_ = aa
    bb_ = bb
    if (forcewyc .and. sys(isys)%c%spgavail) then
       ! calculate the point group operations
       call sys(isys)%c%pointgroup(leqv,rotm)

       ! iteratively find a set of cell parameters that satisfies the user's input
       n = 0
       da = 1d30
       do while (da > epsconv)
          n = n + 1
          if (n > maxstep) exit
          g = m_x2c_from_cellpar(aa_,bb_)
          g = matmul(transpose(g),g)

          ! calculate the symmetrizing matrix
          gavg = 0d0
          do i = 1, leqv
             gavg = gavg + matmul(transpose(rotm(:,:,i)),matmul(g,rotm(:,:,i)))
          end do
          gavg = gavg / leqv

          call cellpar_from_metric(gavg,aa_,bb_)

          if (abs(aa(1) - sys(isys)%c%aa(1)) > tighteps) then
             da = abs(aa(1) - aa_(1))
             aa_(1) = aa(1)
          elseif (abs(aa(2) - sys(isys)%c%aa(2)) > tighteps) then
             da = abs(aa(2) - aa_(2))
             aa_(2) = aa(2)
          elseif (abs(aa(3) - sys(isys)%c%aa(3)) > tighteps) then
             da = abs(aa(3) - aa_(3))
             aa_(3) = aa(3)
          elseif (abs(bb(1) - sys(isys)%c%bb(1)) > tighteps) then
             da = abs(bb(1) - bb_(1))
             bb_(1) = bb(1)
          elseif (abs(bb(2) - sys(isys)%c%bb(2)) > tighteps) then
             da = abs(bb(2) - bb_(2))
             bb_(2) = bb(2)
          elseif (abs(bb(3) - sys(isys)%c%bb(3)) > tighteps) then
             da = abs(bb(3) - bb_(3))
             bb_(3) = bb(3)
          end if
       end do

       ! do nothing if a cell could not be found
       if (da > epsconv) return
    end if

    ! move the unit cell
    copybonding_ = .false.
    if (present(copybonding)) copybonding_ = copybonding
    call sys(isys)%c%move_cell_all(aa_,bb_,copybonding=copybonding_)

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine move_cell

  !> Transform the unit cell to a standardized/reduced cell. mode is one
  !> of the celltransform_* constants (standard, primitive, primstd,
  !> niggli, delaunay). refine applies only to the standardized cells and
  !> refines the atomic positions to their ideal symmetric locations.
  module subroutine transform_cell(sysc,mode,refine,errmsg,keepcell)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: mode
    logical, intent(in) :: refine
    character(len=:), allocatable, intent(inout) :: errmsg
    logical, intent(in), optional :: keepcell

    integer :: isys
    real*8 :: x0(3,3)
    logical :: keepcell_

    errmsg = ""
    keepcell_ = .false.
    if (present(keepcell)) keepcell_ = keepcell

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    ! apply the transformation
    x0 = 0d0
    select case (mode)
    case (celltransform_standard)
       x0 = sys(isys)%c%cell_standard(.false.,.false.,refine,errmsg=errmsg,keepcell=keepcell_)
    case (celltransform_primitive)
       x0 = sys(isys)%c%cell_standard(.true.,.false.,refine,errmsg=errmsg)
    case (celltransform_primstd)
       x0 = sys(isys)%c%cell_standard(.true.,.true.,refine,errmsg=errmsg)
    case (celltransform_niggli)
       x0 = sys(isys)%c%cell_standard(.true.,.true.,.false.,errmsg=errmsg)
       if (.not.has_errmsg(errmsg)) x0 = sys(isys)%c%cell_niggli(errmsg=errmsg)
    case (celltransform_delaunay)
       x0 = sys(isys)%c%cell_standard(.true.,.true.,.false.,errmsg=errmsg)
       if (.not.has_errmsg(errmsg)) x0 = sys(isys)%c%cell_delaunay(errmsg=errmsg)
    case default
       return
    end select
    if (has_errmsg(errmsg)) return

    ! the geometry has changed (only if the cell actually changed)
    if (any(abs(x0) > 1d-5)) &
       call sysc%post_event(lastchange_geometry)

  end subroutine transform_cell

  !> Transform the unit cell using an arbitrary 3x3 transformation matrix
  !> x0 (lattice vectors of the new cell in the old setting, crystallographic
  !> coordinates) and origin shift t0 (fractional). If doinv, use the inverse
  !> of x0 as the transformation matrix.
  module subroutine transform_cell_matrix(sysc,x0,t0,doinv,errmsg)
    use tools_math, only: matinv
    class(sysconf), intent(inout) :: sysc
    real*8, intent(in) :: x0(3,3), t0(3)
    logical, intent(in) :: doinv
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: isys
    real*8 :: x0_(3,3)

    errmsg = ""

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    x0_ = x0
    if (doinv) call matinv(x0_,3)

    ! apply the transformation
    call sys(isys)%c%newcell(x0_,t0,errmsg=errmsg)
    if (has_errmsg(errmsg)) return

    ! the geometry has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine transform_cell_matrix

  !> Read-only passthrough to crystal%cell_nice_list: search for nice
  !> supercells of increasing size up to inice times the input cell and
  !> return the inscribed-sphere radii (rmax) and transformation matrices
  !> (mmax) for the nicest cell of each size.
  module subroutine cell_nice_list(sysc,inice,rmax,mmax)
    class(sysconf), intent(inout) :: sysc
    integer, intent(in) :: inice
    real*8, allocatable, intent(out) :: rmax(:)
    real*8, allocatable, intent(out) :: mmax(:,:,:)

    integer :: isys

    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    call sys(isys)%c%cell_nice_list(inice,rmax,mmax)

  end subroutine cell_nice_list

  !> Recalculate the crystal symmetry using the tolerance stored in
  !> sysc%symeps, then post a geometry-change event.
  module subroutine recalc_symmetry(sysc,errmsg)
    use global, only: symprec
    class(sysconf), intent(inout) :: sysc
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: isys
    real*8 :: osp

    errmsg = ""

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    ! recalculate the symmetry at the chosen tolerance
    osp = symprec
    symprec = sysc%symeps
    call sys(isys)%c%calcsym(.false.,errmsg)
    symprec = osp
    if (len_trim(errmsg) > 0) return

    ! the symmetry (and the non-equivalent atom list) has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine recalc_symmetry

  !> Clear the crystal symmetry (revert to P1), then post a
  !> geometry-change event.
  module subroutine clear_symmetry(sysc)
    class(sysconf), intent(inout) :: sysc

    integer :: isys

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    ! clear the symmetry
    call sys(isys)%clearsym()

    ! the symmetry (and the non-equivalent atom list) has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine clear_symmetry

  !> Delete the symmetry operations flagged in del and rebuild the crystal with
  !> the largest subgroup avoiding them, then post a geometry-change event.
  module subroutine reduce_symmetry(sysc,del,errmsg)
    class(sysconf), intent(inout) :: sysc
    logical, intent(in) :: del(:)
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: isys

    errmsg = ""

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    ! reduce the symmetry to the largest subgroup avoiding the deleted operations
    call sys(isys)%c%reduce_symmetry(del,errmsg)
    if (len_trim(errmsg) > 0) return

    ! the symmetry (and the non-equivalent atom list) has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine reduce_symmetry

  !> Refine the geometry to the ideal symmetry positions (idealize the
  !> cell parameters and atomic positions) using the tolerance stored
  !> in sysc%symeps, keeping the original cell choice (i.e. without
  !> switching to the conventional cell).
  module subroutine refine_symmetry(sysc,errmsg)
    use global, only: symprec
    class(sysconf), intent(inout) :: sysc
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: isys
    real*8 :: osp

    errmsg = ""

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    ! refine at the chosen tolerance, keeping the original cell
    ! (transform_cell posts the event)
    osp = symprec
    symprec = sysc%symeps
    call sysc%transform_cell(celltransform_standard,.true.,errmsg,keepcell=.true.)
    symprec = osp

  end subroutine refine_symmetry

  !> Re-assign atomic types so the asymmetric unit contains whole
  !> molecules, then post a geometry-change event. Only valid for a 3d
  !> molecular crystal.
  module subroutine wholemols_op(sysc,errmsg)
    class(sysconf), intent(inout) :: sysc
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: isys

    errmsg = ""

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return
    if (.not.sys(isys)%c%ismol3d) then
       errmsg = "WHOLEMOLS can only be applied to a molecular crystal"
       return
    end if

    ! apply wholemols
    call sys(isys)%c%wholemols()

    ! the non-equivalent atom list has changed
    call sysc%post_event(lastchange_geometry)

  end subroutine wholemols_op

  !> Scan the space group as a function of the symmetry
  !> tolerance. Returns, for each of a logarithmic series of
  !> tolerances (eps), the international symbol (sym) and the space
  !> group number (num). The system's symmetry is not modified.
  module subroutine spg_analysis(sysc,eps,sym,num)
    use global, only: symprec
    use spglib, only: SpglibDataset
    class(sysconf), intent(inout) :: sysc
    real*8, allocatable, intent(inout) :: eps(:)
    character(len=11), allocatable, intent(inout) :: sym(:)
    integer, allocatable, intent(inout) :: num(:)

    integer :: isys, i, n
    real*8 :: osp
    type(SpglibDataset) :: spg
    character(len=:), allocatable :: errmsg

    real*8, parameter :: spmin = 1d-10
    real*8, parameter :: factor = 10d0
    integer, parameter :: nstep = 11

    if (allocated(eps)) deallocate(eps)
    if (allocated(sym)) deallocate(sym)
    if (allocated(num)) deallocate(num)

    ! consistency checks
    isys = sysc%id
    if (.not.ok_system(isys,sys_init)) return
    if (sys(isys)%c%ismolecule) return

    ! scan the tolerance range, recording the space group at each step
    osp = symprec
    allocate(eps(nstep),sym(nstep),num(nstep))
    n = 0
    symprec = spmin / factor
    do i = 1, nstep
       symprec = symprec * factor
       call sys(isys)%c%spglib_wrap(spg,.false.,errmsg)
       if (len_trim(errmsg) == 0) then
          n = n + 1
          eps(n) = symprec
          sym(n) = spg%international_symbol
          num(n) = spg%spacegroup_number
       end if
    end do
    symprec = osp

    ! trim to the number of successful steps
    if (n == 0) then
       deallocate(eps,sym,num)
    elseif (n < nstep) then
       eps = eps(1:n)
       sym = sym(1:n)
       num = num(1:n)
    end if

  end subroutine spg_analysis

  !xx! private procedures

  ! Thread worker: run over all systems and initialize the ones that are not locked
  function initialization_thread_worker(arg) bind(c)
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

                ! post the first geometry change
                call sysc(i)%post_event(lastchange_geometry,keepfields=.true.)

                ! this system has been initialized
                sysc(i)%status = sys_init

                ! seed the undo/redo history with the initial geometry (must
                ! come after the status is set, since undo_capture only records
                ! systems that are sys_init)
                call sysc(i)%undo_reset()
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
