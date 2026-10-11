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

! Scene files (.scene): a JSON file that recreates a system as it is
! shown in the main view. It holds the structure in full (with the
! edits, the symmetry and the bonds), the bonding settings of the
! system, the Display, the background and the camera. The fields
! are not in the file: it has the commands that made them (the files
! they read are found again next to the scene file, or where they
! were), and they are made again when the system is built. Their
! critical points go to checkpoint files (<field file>.chk_cps, or
! next to the scene file for a field with no file), written when the
! scene is saved and read when it is loaded.
submodule (systems) scenefile
  implicit none

  ! version of the scene file format written; a reader refuses files
  ! with a later version
  integer, parameter :: scene_file_version = 1

  !xx! private procedures
  ! subroutine write_fields(json,p,isys,file)
  ! subroutine read_fields(json,p,file,pend)
  ! subroutine add_file_entry(json,r,path,dir)
  ! function find_file_entry(json,r,dir,changed) result(file)
  ! function replace_files(line,old,new) result(lnew)
  ! subroutine write_objects(json,p,isys)
  ! subroutine load_objects(isys,str,autoadded)

contains

  !> Write the system isys, as its main view shows it, to the scene
  !> file file. errmsg is non-empty if it could not be written.
  module subroutine scene_file_save(isys,file,errmsg)
    use json_module, only: json_value, json_core
    use config, only: getstring, istring_version
    integer, intent(in) :: isys
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg

    type(json_core), target :: json
    type(json_value), pointer :: p, o
    type(crystalseed) :: seed
    type(scene_camera) :: cam
    integer :: i, j, n
    integer, allocatable :: iv(:)
    logical :: ok, copysym

    errmsg = ""
    if (.not.ok_system(isys,sys_init)) then
       errmsg = "the system is not ready"
       return
    end if

    call json%initialize(compress_vectors=.true.)
    call json%create_object(p,'')
    call json%add(p,'critic2_scene',scene_file_version)
    call json%add(p,'program',"critic2 " // getstring(istring_version))
    call json%add(p,'source',sysc(isys)%fullname)

    ! the system settings
    call json%create_object(o,'system')
    call json%add(p,o)
    call json%add(o,'name',trim(sysc(isys)%seed%name))
    call json%add(o,'renamed',sysc(isys)%renamed)
    call json%add(o,'symeps',sysc(isys)%symeps)
    call json%add(o,'bondfactor',sysc(isys)%bondfactor)
    call json%add(o,'bonddelta',sysc(isys)%bonddelta)
    call json%add(o,'atmcov',sysc(isys)%atmcov)
    if (allocated(sysc(isys)%bondallowed)) then
       ! the pairs of atomic numbers that may not bond
       n = count(.not.sysc(isys)%bondallowed)
       allocate(iv(2*n))
       n = 0
       do j = 1, size(sysc(isys)%bondallowed,2)
          do i = 1, size(sysc(isys)%bondallowed,1)
             if (sysc(isys)%bondallowed(i,j)) cycle
             iv(n+1:n+2) = (/i,j/)
             n = n + 2
          end do
       end do
       call json%add(o,'bondforbidden',iv)
    end if

    ! the structure, as undo saves it: the same crystal is built from it
    copysym = (.not.sys(isys)%c%ismolecule .and. sys(isys)%c%spgavail)
    call sys(isys)%c%makeseed(seed,copysym=copysym,copybonding=.true.)
    seed%name = sysc(isys)%seed%name
    call seed%write_json(json,p)

    ! the fields and the objects
    call write_fields(json,p,isys,file)
    call write_objects(json,p,isys)

    ! the scene: Display, background, camera
    call sysc(isys)%sc%disp%write_json(json,p)
    call json%add(p,'background',real(sysc(isys)%sc%bgcolor,8))
    call sysc(isys)%sc%cam_get(cam)
    call cam%write_json(json,p)

    ! write the file
    call json%print(p,file)
    call json%check_for_errors(ok,errmsg)
    if (ok) errmsg = ""
    call json%destroy(p)

  end subroutine scene_file_save

  !> Read the scene file file and create the system it holds, with
  !> the settings of the file (those that need the built system wait
  !> in its scenefile, for scene_file_apply). The caller launches the
  !> initialization thread. errmsg is non-empty if no system was
  !> created.
  module subroutine scene_file_load(file,errmsg)
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg

    type(crystalseed), allocatable :: seed(:)
    type(scene_pending) :: pend

    call scene_file_read(file,seed,pend,errmsg)
    if (len_trim(errmsg) > 0) return
    call add_systems_from_seeds(1,seed,scenefile=pend)

  end subroutine scene_file_load

  !> Read the scene file file: the structure of the system it holds
  !> (seed(1)) and its settings (pend). errmsg is non-empty if the file
  !> could not be read.
  module subroutine scene_file_read(file,seed,pend,errmsg)
    use json_module, only: json_value, json_core
    use tools_io, only: string
    use utils, only: file_name_abs
    use param, only: isformat_r_derived
    use representations, only: reptype_NUM
    character*(*), intent(in) :: file
    type(crystalseed), allocatable, intent(inout) :: seed(:)
    type(scene_pending), intent(out) :: pend
    character(len=:), allocatable, intent(out) :: errmsg

    type(json_core) :: json
    type(json_value), pointer :: p, o
    integer, allocatable :: iv(:)
    real*8, allocatable :: rv(:)
    logical, allocatable :: lv(:)
    character(len=:), allocatable :: name
    integer :: i, iver
    logical :: found, ok

    errmsg = ""
    nullify(p)
    call json%initialize()
    call json%load(file=file,p=p)
    if (json%failed()) then
       call json%check_for_errors(ok,errmsg)
       errmsg = "not a valid JSON file: " // errmsg
       goto 999
    end if

    ! the format version
    call json%get(p,'critic2_scene',iver,found)
    if (.not.found) then
       errmsg = "not a critic2 scene file"
       goto 999
    elseif (iver > scene_file_version) then
       errmsg = "scene file version " // string(iver) // " is newer than this critic2 reads (" //&
          string(scene_file_version) // ")"
       goto 999
    end if

    ! the structure
    if (allocated(seed)) deallocate(seed)
    allocate(seed(1))
    call json%get(p,'structure',o,found)
    if (.not.found) then
       errmsg = "the scene file has no structure"
       goto 999
    end if
    call seed(1)%read_json(json,o,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! the system settings
    call json%get(p,'system',o,found)
    if (found) then
       call json%get(o,'renamed',pend%renamed,found,default=.false.)
       call json%get(o,'name',name,found)
       if (found .and. pend%renamed) seed(1)%name = name
       call json%get(o,'symeps',pend%symeps,found,default=-1d0)
       call json%get(o,'bondfactor',pend%bondfactor,found,default=-1d0)
       call json%get(o,'bonddelta',pend%bonddelta,found,default=-1d0)
       call json%get(o,'atmcov',rv,found)
       if (found) then
          if (size(rv) == maxzat0+1) pend%atmcov = rv
       end if
       call json%get(o,'bondforbidden',iv,found)
       if (found) then
          allocate(pend%bondallowed(maxzat0,maxzat0))
          pend%bondallowed = .true.
          do i = 1, size(iv)/2
             if (any(iv(2*i-1:2*i) < 1 .or. iv(2*i-1:2*i) > maxzat0)) cycle
             pend%bondallowed(iv(2*i-1),iv(2*i)) = .false.
          end do
       end if
    end if

    ! the system is the scene: it is known by the file name unless it
    ! was renamed, and it is an in-memory system (Save asks for a file,
    ! re-reading rebuilds it from this seed)
    if (.not.pend%renamed) seed(1)%name = file
    seed(1)%file = ""
    seed(1)%isformat = isformat_r_derived

    ! the fields
    call read_fields(json,p,file,pend)

    ! the objects, kept as JSON text until the system is built
    call json%get(p,'objects',o,found)
    if (found) call json%print_to_string(o,pend%objects)
    call json%get(p,'objects_autoadded',lv,found)
    if (found) then
       if (size(lv) == reptype_NUM) pend%autoadded = lv
    end if

    ! the scene settings
    call json%get(p,'display',o,found)
    if (found) then
       allocate(pend%disp)
       call pend%disp%read_json(json,o)
    end if
    call json%get(p,'background',rv,found)
    if (found) then
       if (size(rv) == 3) pend%bgcolor = real(rv,c_float)
    end if
    call json%get(p,'camera',o,found)
    if (found) then
       allocate(pend%cam)
       if (.not.pend%cam%read_json(json,o)) deallocate(pend%cam)
    end if

    pend%file = file_name_abs(file)

999 continue
    if (json%failed()) call json%clear_exceptions()
    if (associated(p)) call json%destroy(p)

  end subroutine scene_file_read

  !> Make the fields of the scene file waiting in the system isys, in
  !> the slots they had, with their names and options, and set the
  !> reference field. Run by the initialization thread (ti) after the
  !> system is built; problems are reported to the output console.
  module subroutine scene_file_load_fields(isys,ti)
    use tools_io, only: uout, string
    integer, intent(in) :: isys
    type(thread_info), intent(in), optional :: ti

    integer :: k, id, nleft, nleft0
    character(len=:), allocatable :: errmsg
    logical :: found
    logical, allocatable :: done(:)
    character(len=:), allocatable :: errlast(:)

    if (.not.allocated(sysc(isys)%scenefile)) return
    if (.not.allocated(sysc(isys)%scenefile%fields)) return
    if (.not.sys(isys)%isinit) return

    associate(fields => sysc(isys)%scenefile%fields)
      do k = 1, size(fields)
         if (allocated(fields(k)%warn)) &
            write (uout,'("!! Warning !! ",A)') fields(k)%warn
      end do

      ! Make the fields in the order they were made. A field that
      ! fails is tried again after the others, as long as some field
      ! is made in each pass: a ghost field may use fields that were
      ! made after it (in a slot that was freed and filled again).
      allocate(done(size(fields)))
      allocate(character(len=mlen) :: errlast(size(fields)))
      done = .false.
      errlast = ""
      nleft = size(fields)
      nleft0 = nleft + 1
      do while (nleft > 0 .and. nleft < nleft0)
         nleft0 = nleft
         do k = 1, size(fields)
            if (done(k)) cycle

            ! make the field (field 0 is always there)
            if (fields(k)%id == 0) then
               id = 0
            else
               call sys(isys)%load_field_string(fields(k)%load,.false.,id,errmsg,ti=ti,&
                  idwant=fields(k)%id)
               if (len_trim(errmsg) > 0 .or. .not.sys(isys)%goodfield(id)) then
                  errlast(k) = errmsg
                  cycle
               end if
            end if
            done(k) = .true.
            nleft = nleft - 1

            ! name and options
            if (allocated(fields(k)%name)) sys(isys)%f(id)%name = fields(k)%name
            if (allocated(fields(k)%options)) then
               call sys(isys)%set_field_options(id,fields(k)%options,errmsg)
               if (len_trim(errmsg) > 0) &
                  write (uout,'("!! Warning !! Could not set the options of field ",A," of the scene: ",A)') &
                  string(id), trim(errmsg)
            end if

            ! the critical points, after the options (the interpolation
            ! is in the checkpoint's fingerprint)
            if (allocated(fields(k)%cpfile)) then
               call sys(isys)%f(id)%read_chk_cps(fields(k)%cpfile,found,errmsg,ti=ti)
               if (.not.found) errmsg = "file not found"
               if (len_trim(errmsg) > 0) then
                  write (uout,'("!! Warning !! Could not read the critical points of field ",A," of the scene from ",A,": ",A)') &
                     string(id), fields(k)%cpfile, trim(errmsg)
               else
                  write (uout,'("* Read ",A," critical points for field ",A," from: ",A)') &
                     string(sys(isys)%f(id)%ncp), string(id), fields(k)%cpfile
               end if
            end if
         end do
      end do

      ! the fields that could not be made
      do k = 1, size(fields)
         if (done(k)) cycle
         write (uout,'("!! Warning !! Could not load field ",A," of the scene: ",A)') &
            string(fields(k)%id), trim(errlast(k))
         write (uout,'("!! Warning !! Command: LOAD ",A)') fields(k)%load
      end do
    end associate

    ! the reference field
    k = sysc(isys)%scenefile%iref
    if (k >= 0 .and. k /= sys(isys)%iref) then
       if (sys(isys)%goodfield(k)) call sys(isys)%set_reference(k,.false.)
    end if

  end subroutine scene_file_load_fields

  !> Apply the scene file settings waiting in the system isys to its
  !> scene, once the system is built (by the initialization thread,
  !> after the first geometry change, so the Display masks are stamped
  !> with the bonding revision they belong to). The camera is taken by
  !> the scene after its first camera reset.
  module subroutine scene_file_apply(isys)
    integer, intent(in) :: isys

    if (.not.ok_system(isys,sys_ready)) return
    if (.not.allocated(sysc(isys)%scenefile)) return

    associate(pend => sysc(isys)%scenefile, sc => sysc(isys)%sc)
      if (allocated(pend%objects)) call load_objects(isys,pend%objects,pend%autoadded)
      if (allocated(pend%disp)) then
         sc%disp = pend%disp
         if (sc%disp%atype < 1 .or. sc%disp%atype > atlisttype_NUM) then
            ! an unknown grouping: its masks mean nothing
            sc%disp%atype = atlisttype_species
            if (allocated(sc%disp%ashown)) deallocate(sc%disp%ashown)
         end if
         sc%disp%rev = sysc(isys)%rev_rebond
         call sc%disp%update(isys)
      end if
      if (allocated(pend%bgcolor)) sc%bgcolor = pend%bgcolor
      if (allocated(pend%cam)) sc%campending = pend%cam
      sc%forcebuildlists = .true.
    end associate
    deallocate(sysc(isys)%scenefile)

  end subroutine scene_file_apply

  !xx! private procedures

  !> Write the fields of system isys to the json object p, as the
  !> commands that made them, in the order they were made, and the
  !> reference field. The files they read are written as given, as
  !> absolute paths, and relative to the directory of the scene file
  !> (file) if they are in it, with their sizes. Fields that cannot be
  !> made again are left out, with a warning.
  subroutine write_fields(json,p,isys,file)
    use json_module, only: json_value, json_core
    use tools_io, only: uout, string
    use utils, only: file_name_abs, file_name_dir, file_name_root
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(in) :: p
    integer, intent(in) :: isys
    character*(*), intent(in) :: file

    type(json_value), pointer :: o, q, a, r
    integer :: i, j, k, n, iref
    integer, allocatable :: ord(:)
    character(len=:), allocatable :: dir, root, cpfile, errmsg
    character(len=mlen), allocatable :: cpdone(:)
    integer :: ncpdone

    associate(sy => sys(isys))
      ! the fields, in the order they were made (the order of their
      ! identifiers): a field may be made from the ones before it
      allocate(ord(sy%nf+1))
      n = 0
      do i = 0, sy%nf
         if (.not.sy%goodfield(i)) cycle
         j = n
         do while (j > 0)
            if (sy%f(ord(j))%uid <= sy%f(i)%uid .or. ord(j) == 0) exit
            ord(j+1) = ord(j)
            j = j - 1
         end do
         ord(j+1) = i
         n = n + 1
      end do

      dir = file_name_dir(file_name_abs(file))
      allocate(cpdone(n))
      ncpdone = 0
      iref = -1
      call json%create_array(o,'fields')
      call json%add(p,o)
      do k = 1, n
         i = ord(k)
         if (i > 0 .and. .not.allocated(sy%f(i)%loadcmd)) then
            write (uout,'("!! Warning !! Field ",A," (",A,") is not in the scene: it was not loaded from a file or a command")') &
               string(i), trim(sy%f(i)%name)
            cycle
         end if
         call json%create_object(q,'')
         call json%add(q,'id',i)
         call json%add(q,'name',trim(sy%f(i)%name))
         if (i > 0) then
            call json%add(q,'load',sy%f(i)%loadcmd)
            if (allocated(sy%f(i)%loadfiles)) then
               call json%create_array(a,'files')
               call json%add(q,a)
               do j = 1, size(sy%f(i)%loadfiles)
                  call json%create_object(r,'')
                  call add_file_entry(json,r,trim(sy%f(i)%loadfiles(j)),dir)
                  call json%add(a,r)
               end do
            end if
         end if
         if (allocated(sy%f(i)%optcmd)) call json%add(q,'options',sy%f(i)%optcmd)

         ! the critical points (beyond the nuclei), in a checkpoint:
         ! the field's own (<field file>.chk_cps), or one next to the
         ! scene file if the field has no file, that one cannot be
         ! written, or another field of the scene has it (same file)
         if (sy%f(i)%ncp > sy%c%nneq) then
            root = file_name_root(file) // "_field" // string(i)
            cpfile = sy%f(i)%chk_cps_file(root)
            errmsg = "in use"
            if (.not.any(cpdone(1:ncpdone) == cpfile)) call sy%f(i)%write_chk_cps(cpfile,errmsg)
            if (len_trim(errmsg) > 0 .and. cpfile /= root // ".chk_cps") then
               cpfile = root // ".chk_cps"
               call sy%f(i)%write_chk_cps(cpfile,errmsg)
            end if
            if (len_trim(errmsg) > 0) then
               write (uout,'("!! Warning !! The critical points of field ",A," are not in the scene: ",A)') &
                  string(i), trim(errmsg)
            else
               write (uout,'("* Critical points of field ",A," written to the checkpoint: ",A)') &
                  string(i), cpfile
               call json%create_object(r,'cpfile')
               call add_file_entry(json,r,cpfile,dir)
               call json%add(q,r)
               ncpdone = ncpdone + 1
               cpdone(ncpdone) = cpfile
            end if
         end if
         call json%add(o,q)
         if (i == sy%iref) iref = i
      end do
      if (iref >= 0) call json%add(p,'reference_field',iref)
    end associate

  end subroutine write_fields

  !> Write the objects of the main scene of system isys (not the
  !> transient ones) to the json object p, in the order of the object
  !> list, with the kinds added to the scene. A color bar names its
  !> isosurface object by its position in the file (from 1), and the
  !> isosurface by its position in that object (colorbar_link).
  subroutine write_objects(json,p,isys)
    use json_module, only: json_value, json_core
    use jsonxfer, only: json_xfer
    use representations, only: reptype_colorbar, reptype_isosurface, rep_slot_of, iso_slot_index
    use tools, only: mergesort
    type(json_core), intent(inout), target :: json
    type(json_value), pointer, intent(in) :: p
    integer, intent(in) :: isys

    type(json_xfer) :: io
    type(json_value), pointer :: a, q, r
    integer :: i, k, n, irep, iobj, islot
    logical :: found
    integer, allocatable :: iord(:), perm(:), ifile(:)

    io%json => json
    io%reading = .false.
    associate(sc => sysc(isys)%sc)
      ! the objects, in list order
      allocate(iord(sc%nrep),perm(sc%nrep))
      do i = 1, sc%nrep
         iord(i) = sc%rep(i)%iord
         perm(i) = i
      end do
      if (sc%nrep > 0) call mergesort(iord,perm,1,sc%nrep)

      ! position of each object in the file (0 = not in it)
      allocate(ifile(sc%nrep))
      ifile = 0
      n = 0
      do k = 1, sc%nrep
         i = perm(k)
         if (.not.sc%rep(i)%isinit .or. sc%rep(i)%owner /= 0) cycle
         n = n + 1
         ifile(i) = n
      end do

      call json%create_array(a,'objects')
      call json%add(p,a)
      do k = 1, sc%nrep
         i = perm(k)
         if (ifile(i) == 0) cycle
         call json%create_object(q,'')
         call io%x(q,'type',sc%rep(i)%type)
         call io%x(q,'flavor',sc%rep(i)%flavor)
         call sc%rep(i)%xfer_json(io,q,.false.)
         call sc%rep(i)%xfer_json(io,q,.true.)

         ! the isosurface of a color bar (none = automatic)
         if (sc%rep(i)%type == reptype_colorbar) then
            irep = rep_slot_of(sc%rep,sc%nrep,sc%rep(i)%colorbar%isouid)
            if (irep > 0) then
               if (ifile(irep) > 0 .and. sc%rep(irep)%type == reptype_isosurface) then
                  iobj = ifile(irep)
                  islot = iso_slot_index(sc%rep(irep)%iso,sc%rep(i)%colorbar%slotuid)
                  call io%obj(q,'colorbar_link',r,found)
                  call io%x(r,'object',iobj)
                  call io%x(r,'isosurface',islot)
               end if
            end if
         end if
         call json%add(a,q)
      end do
      call json%add(p,'objects_autoadded',sc%autoadded)
    end associate

  end subroutine write_objects

  !> Replace the objects of the main scene of system isys with those
  !> of the scene file, given as the JSON text str of its "objects"
  !> array, and set the kinds added to the scene (autoadded) if given.
  !> Each object is made with the defaults of its kind, its options
  !> are read, its styles made again for the system with them, and
  !> then its style rows read. Invalid entries are skipped. The color
  !> bars are linked to their isosurfaces after all the objects are
  !> made. Run by the initialization thread, before any window can
  !> point to the objects.
  subroutine load_objects(isys,str,autoadded)
    use json_module, only: json_value, json_core
    use jsonxfer, only: json_xfer
    use tools_io, only: uout
    use representations, only: representation, reptype_NUM, repflavor_NUM, repflavor_type,&
       reptype_colorbar, reptype_isosurface
    integer, intent(in) :: isys
    character*(*), intent(in) :: str
    logical, allocatable, intent(in) :: autoadded(:)

    type(json_core), target :: json
    type(json_xfer) :: io
    type(json_value), pointer :: a, q, r
    type(representation), allocatable :: aux(:)
    integer :: i, n, itype, iflavor, id, k, islot
    integer, allocatable :: idnew(:)
    logical :: ok, found
    character(len=:), allocatable :: errmsg

    nullify(a)
    call json%initialize()
    call json%deserialize(a,str)
    if (json%failed()) then
       call json%check_for_errors(ok,errmsg)
       write (uout,'("!! Warning !! Could not read the objects of the scene: ",A)') errmsg
       call json%clear_exceptions()
       if (associated(a)) call json%destroy(a)
       return
    end if
    n = json%count(a)

    io%json => json
    io%reading = .true.
    associate(sc => sysc(isys)%sc)
      ! the objects the scene was made with go, and their names and
      ! list positions with them
      do i = 1, sc%nrep
         if (sc%rep(i)%isinit .and. sc%rep(i)%owner == 0) call sc%rep(i)%end()
      end do
      sc%icount = 0

      ! room for all the objects, so that adding them does not grow the
      ! list (that would make the windows point to the objects again,
      ! from this thread)
      if (size(sc%rep,1) < sc%nrep + n) then
         allocate(aux(sc%nrep + n))
         aux(1:sc%nrep) = sc%rep(1:sc%nrep)
         call move_alloc(aux,sc%rep)
      end if

      ! the objects of the file, and the slot each is made in (0 = not
      ! made)
      allocate(idnew(n))
      idnew = 0
      do i = 1, n
         call json%get_child(a,i,q)
         if (.not.associated(q)) cycle
         itype = 0
         iflavor = -1
         call io%x(q,'type',itype,1,reptype_NUM)
         call io%x(q,'flavor',iflavor,0,repflavor_NUM)
         if (itype == 0 .or. iflavor < 0) cycle
         if (repflavor_type(iflavor) /= itype) cycle
         call sc%add_representation(itype,iflavor,id)
         if (.not.sc%rep(id)%isinit) cycle
         call sc%rep(id)%xfer_json(io,q,.false.)
         call sc%rep(id)%reset_all_styles(0)
         call sc%rep(id)%update()
         call sc%rep(id)%xfer_json(io,q,.true.)
         idnew(i) = id
      end do

      ! the isosurfaces of the color bars
      do i = 1, n
         if (idnew(i) == 0) cycle
         if (sc%rep(idnew(i))%type /= reptype_colorbar) cycle
         call json%get_child(a,i,q)
         call io%obj(q,'colorbar_link',r,found)
         if (.not.found) cycle
         k = 0
         islot = 0
         call io%x(r,'object',k,1,n)
         call io%x(r,'isosurface',islot,0,huge(1))
         if (k == 0) cycle
         if (idnew(k) == 0) cycle
         if (sc%rep(idnew(k))%type /= reptype_isosurface) cycle
         sc%rep(idnew(i))%colorbar%isouid = sc%rep(idnew(k))%uid
         if (islot >= 1 .and. islot <= sc%rep(idnew(k))%iso%niso) &
            sc%rep(idnew(i))%colorbar%slotuid = sc%rep(idnew(k))%iso%slot(islot)%uid
      end do
      if (allocated(autoadded)) sc%autoadded = autoadded
      sc%forcesort = .true.
      sc%forcebuildlists = .true.
    end associate
    if (json%failed()) call json%clear_exceptions()
    call json%destroy(a)

  end subroutine load_objects

  !> Read the fields of the scene file file from the json object p
  !> into pend. The files each field reads are looked for next to the
  !> scene file, then where they were, and the command is changed to
  !> read them from where they are found. Invalid entries are
  !> skipped.
  subroutine read_fields(json,p,file,pend)
    use json_module, only: json_value, json_core
    use tools_io, only: string
    use utils, only: file_name_dir
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(in) :: p
    character*(*), intent(in) :: file
    type(scene_pending), intent(inout) :: pend

    integer, parameter :: maxid = 1000 ! largest field slot accepted

    type(json_value), pointer :: o, q, a, r
    type(scene_field), allocatable :: fl(:)
    integer :: i, j, n, nfl, nfi, id
    character(len=:), allocatable :: dir, str, path, warn
    character(len=mlen), allocatable :: fold(:), fnew(:)
    logical :: found, changed

    call json%get(p,'fields',o,found)
    if (.not.found) return
    n = json%count(o)
    if (n <= 0) return
    dir = file_name_dir(file)

    allocate(fl(n))
    nfl = 0
    do i = 1, n
       call json%get_child(o,i,q)
       if (.not.associated(q)) cycle
       fl(nfl+1) = scene_field()

       ! slot and command; skip repeated or unknown slots
       call json%get(q,'id',id,found)
       if (.not.found) cycle
       if (id < 0 .or. id > maxid) cycle
       if (any(fl(1:nfl)%id == id)) cycle
       fl(nfl+1)%id = id
       if (id > 0) then
          call json%get(q,'load',str,found)
          if (.not.found) cycle
          if (len_trim(str) == 0) cycle
          fl(nfl+1)%load = str
       end if
       call json%get(q,'name',str,found)
       if (found) fl(nfl+1)%name = str
       call json%get(q,'options',str,found)
       if (found) then
          if (len_trim(str) > 0) fl(nfl+1)%options = str
       end if

       ! the files: next to the scene file, where they were, or as
       ! written in the command
       call json%get(q,'files',a,found)
       if (found .and. id > 0) then
          nfi = json%count(a)
          if (allocated(fold)) deallocate(fold)
          if (allocated(fnew)) deallocate(fnew)
          allocate(fold(max(nfi,0)),fnew(max(nfi,0)))
          warn = ""
          do j = 1, nfi
             call json%get_child(a,j,r)
             fold(j) = ""
             fnew(j) = ""
             if (.not.associated(r)) cycle
             call json%get(r,'path',path,found)
             if (.not.found) cycle
             fold(j) = path
             fnew(j) = find_file_entry(json,r,dir,changed)
             if (changed) then
                if (len(warn) > 0) warn = warn // ", "
                warn = warn // trim(fnew(j))
             end if
          end do
          if (len(warn) > 0) &
             fl(nfl+1)%warn = "Field " // string(id) // " of the scene: the file changed since the scene was saved: " // warn
          fl(nfl+1)%load = replace_files(fl(nfl+1)%load,fold,fnew)
       end if

       ! the checkpoint with the critical points
       call json%get(q,'cpfile',r,found)
       if (found) then
          path = find_file_entry(json,r,dir,changed)
          if (len_trim(path) > 0) fl(nfl+1)%cpfile = path
          if (changed) then
             str = "Field " // string(id) // " of the scene: the checkpoint changed since the scene was saved: " // path
             if (allocated(fl(nfl+1)%warn)) str = fl(nfl+1)%warn // "; " // str
             fl(nfl+1)%warn = str
          end if
       end if

       nfl = nfl + 1
    end do
    if (nfl > 0) pend%fields = fl(1:nfl)

    call json%get(p,'reference_field',id,found)
    if (found) then
       if (id >= 0 .and. id <= maxid) then
          if (any(fl(1:nfl)%id == id)) pend%iref = id
       end if
    end if

  end subroutine read_fields

  !> Add the file path to the JSON object r, as written (path), as an
  !> absolute path (abs), relative to the directory dir of the scene
  !> file if it is under it (rel), and its size in bytes (size; a real,
  !> json integers are 32-bit), for find_file_entry.
  subroutine add_file_entry(json,r,path,dir)
    use json_module, only: json_value, json_core
    use utils, only: file_name_abs
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(in) :: r
    character*(*), intent(in) :: path, dir

    character(len=:), allocatable :: ap
    integer*8 :: fsize
    logical :: ok

    ap = file_name_abs(path)
    call json%add(r,'path',path)
    call json%add(r,'abs',ap)
    if (len(dir) > 0 .and. len(ap) > len(dir)) then
       if (ap(1:len(dir)) == dir) call json%add(r,'rel',ap(len(dir)+1:))
    end if
    inquire(file=ap,exist=ok,size=fsize)
    if (ok .and. fsize >= 0) call json%add(r,'size',real(fsize,8))

  end subroutine add_file_entry

  !> The file of the JSON object r written by add_file_entry: next to
  !> the scene file (in directory dir), where it was, or as it was
  !> written, the first that exists (the last if none). changed is
  !> .true. if the file has a size different from the one it had when
  !> the scene was saved. Returns "" if r has no path.
  function find_file_entry(json,r,dir,changed) result(file)
    use json_module, only: json_value, json_core
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(in) :: r
    character*(*), intent(in) :: dir
    logical, intent(out) :: changed
    character(len=:), allocatable :: file

    character(len=:), allocatable :: path, rel, ap
    integer*8 :: fsize
    real*8 :: fsize0
    logical :: found, ok

    file = ""
    changed = .false.
    call json%get(r,'path',path,found)
    if (.not.found) return
    call json%get(r,'rel',rel,found)
    if (found .and. len_trim(rel) > 0) then
       file = dir // rel
       inquire(file=file,exist=ok)
       if (.not.ok) file = ""
    end if
    if (len(file) == 0) then
       call json%get(r,'abs',ap,found)
       if (found .and. len_trim(ap) > 0) then
          inquire(file=ap,exist=ok)
          if (ok) file = ap
       end if
    end if
    if (len(file) == 0) file = path

    call json%get(r,'size',fsize0,found)
    if (found) then
       inquire(file=file,exist=ok,size=fsize)
       changed = (ok .and. fsize >= 0 .and. real(fsize,8) /= fsize0)
    end if

  end function find_file_entry

  !> Replace the files old(:) in the command line with new(:): the
  !> first word equal to old(1) after the start of the line, then the
  !> first equal to old(2) after it, and so on (blank old(:) are
  !> skipped). The replaced words are quoted if they need it, and the
  !> rest of the line is unchanged.
  function replace_files(line,old,new) result(lnew)
    use tools_io, only: getword, quoteword
    character*(*), intent(in) :: line
    character(len=mlen), intent(in) :: old(:), new(:)
    character(len=:), allocatable :: lnew

    integer :: lp, lp0, i0, j
    character(len=:), allocatable :: word

    lnew = ""
    lp = 1
    lp0 = 1
    j = 1
    do while (j <= size(old))
       if (len_trim(old(j)) == 0) then
          j = j + 1
          cycle
       end if
       if (lp > len(line)) exit
       ! start of the next word
       i0 = lp - 1 + verify(line(lp:)," ")
       word = getword(line,lp)
       if (len_trim(word) == 0) exit
       if (word == old(j)) then
          if (old(j) /= new(j)) then
             lnew = lnew // line(lp0:i0-1) // quoteword(trim(new(j))) // " "
             lp0 = lp
          end if
          j = j + 1
       end if
    end do
    if (lp0 <= len(line)) lnew = lnew // line(lp0:)
    lnew = trim(lnew)

  end function replace_files

end submodule scenefile
