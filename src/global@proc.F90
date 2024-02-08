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

! Global variables, parameters and basic initialization subroutines.
! Contains global variables and parameters, and initialization
! procedures run at the beginning of the execution.
submodule (global) proc
  implicit none

contains

  !> Main driver routine for critic2
  module subroutine critic_main()
    use tricks, only: trick
    use molcalc, only: molcalc_driver
    use qtree, only: qtree_driver
    use stm, only: stm_driver
    use xdm, only: xdm_driver
    use hirshfeld, only: hirsh_nogrid
    use sigmahole, only: sigmahole_driver
    use bisect, only: basinplot, bundleplot, sphereintegrals, integrals
    use integration, only: intgrid_driver
    use flux, only: fluxprint
    use autocp, only: autocritic, cpreport
    use nci, only: nciplot
    use rhoplot, only: rhoplot_point, rhoplot_line, rhoplot_plane, rhoplot_cube,&
       rhoplot_grdvec
    use fieldmod, only: type_grid
    use struct_drivers, only: struct_crystal_input, struct_molcell,&
       struct_atomlabel, struct_sym, struct_sym, struct_charges, struct_write,&
       struct_powder, struct_rdf, struct_compare, struct_comparevc, struct_environ,&
       struct_econ, struct_edit,&
       struct_coord, struct_polyhedra, struct_packing, struct_vdw, struct_identify,&
       struct_makemols_neighcrys, struct_molreorder, struct_molmove,&
       struct_kpoints, struct_bz, struct_newcell, struct_amd
    use systemmod, only: sy
    use spglib, only: spg_list_spg
    use arithmetic, only: listvariables, listlibxc
    use tools_math, only: emd
    use tools_io, only: ferror, faterr, uin, ucopy, uout, getword, lgetword, getline,&
       equal, isinteger, isreal, string, usegui
    use param, only: param_init

    ! parsing
    integer :: lp, lpold
    character(len=:), allocatable :: subline, word, word2
    character(len=:), allocatable :: line, errmsg
    !
    integer :: id, idum
    integer :: i, nn, ismoli, ncom
    logical :: ok
    real*8 :: rdum
#ifdef HAVE_LIBXC
    logical :: doref, doname, doflags
#endif

    ! Start reading
    ncom = 1
    main: do while (getline(uin,line,ucopy=ucopy,nprompt=ncom))
       ncom = ncom + 1
       lp=1
       word = lgetword(line,lp)
       subline = line(lp:)

       ! crystal
       if (equal(word,'crystal') .or. equal(word,'molecule')) then
          if (usegui) then
             write (uout,'("!! The CRYSTAL and MOLECULE keywords cannot be used in the graphical interface.")')
             write (uout,'("!! Please use File->New or File->Open to add new structures.")')
             cycle main
          end if

          if (equal(word,'crystal')) then
             ismoli = 0
             if (iunit_isdef) iunit = iunit_bohr
          else
             ismoli = 1
             if (iunit_isdef) iunit = iunit_ang
          end if
          call struct_crystal_input(subline,ismoli,.true.,.not.quiet,s0=sy)

          ! newcell
       elseif (equal(word,'edit')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_edit(sy,.not.quiet)

          ! newcell
       elseif (equal(word,'newcell')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_newcell(sy,subline,.not.quiet)

          ! molcell
       elseif (equal(word,'molcell')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_molcell(sy,subline)

          ! sym/symm/nosym/nosymm
       elseif (equal(word,'symm').or.equal(word,'sym').or.equal(word,'nosym').or.equal(word,'nosymm')) then
          if (equal(word,'nosym').or.equal(word,'nosymm')) then
             doguess = 0
          else
             word = lgetword(line,lp)
             if (isinteger(idum,word)) then
                doguess = idum
             elseif (isreal(rdum,word)) then
                symprec = rdum
                call check_structure_defined(ok,silent=.true.)
                if (ok) then
                   if (.not.sy%c%ismolecule) &
                      call struct_sym(sy,"recalc",.not.quiet)
                end if
             else
                call check_structure_defined(ok)
                if (.not.ok) cycle
                call struct_sym(sy,subline,.not.quiet)
             end if
          end if

          ! q/qat, zpsp, nocore
       elseif (equal(word,'q') .or. equal(word,'qat') .or. equal(word,'zpsp') .or. equal(word,'nocore')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_charges(sy,line,ok)

          ! atomlabel
       elseif (equal(word,'atomlabel')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_atomlabel(sy,subline)

          ! write
       elseif (equal(word,'write')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_write(sy,subline,.false.)

          ! load
       elseif (equal(word,'load')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call sy%load_field_string(subline,.true.,id,errmsg)
          if (id < 0 .or. len_trim(errmsg) > 0) &
             call ferror('load',errmsg,faterr,line,syntax=.true.)

          ! unload
       elseif (equal(word,'unload')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          lpold = lp
          word = getword(line,lp)
          nn = sy%fieldname_to_idx(word)
          if (nn >= 0) then
             if (.not.sy%goodfield(nn)) then
                call ferror('critic2','wrong field in UNLOAD',faterr,line,syntax=.true.)
                cycle
             end if
             if (nn == 0) then
                call ferror('critic2','can not unload the promolecular density',faterr,line,syntax=.true.)
                cycle
             end if
             call sy%unload_field(nn)
          else
             lp = lpold
             word = lgetword(line,lp)
             if (equal(word,'all')) then
                do i = 1, ubound(sy%f,1)
                   call sy%unload_field(i)
                end do
             else
                call ferror('critic2','Unknown keyword in UNLOAD',faterr,line,syntax=.true.)
             endif
          end if

          ! powder
       elseif (equal(word,'powder')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_powder(sy,line(lp:))

          ! rdf
       elseif (equal(word,'rdf')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_rdf(sy,line(lp:))

          ! amd
       elseif (equal(word,'amd')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_amd(sy,line(lp:))

          ! compare
       elseif (equal(word,'compare')) then
          call struct_compare(sy,line(lp:))

          ! comparevc
       elseif (equal(word,'comparevc')) then
          call struct_comparevc(sy,line(lp:))

          ! setfield
       elseif (equal(word,'setfield')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          lpold = lp
          word = getword(line,lp)
          id = sy%fieldname_to_idx(word)
          if (id < 0) then
             id = sy%iref
             lp = lpold
          end if
          if (.not.sy%goodfield(id)) then
             call ferror('critic2','wrong field in setfield',faterr,line,syntax=.true.)
             cycle
          end if
          call sy%f(id)%set_options(line(lp:),errmsg)
          if (len_trim(errmsg) > 0) then
             call ferror('setfield',errmsg,faterr,line,syntax=.true.)
             cycle
          end if
          write (uout,'("+ Field number ",A)') string(id)
          call sy%f(id)%printinfo(.false.,.true.)
          write (uout,*)

          ! reference
       elseif (equal(word,'reference')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          word = getword(line,lp)
          id = sy%fieldname_to_idx(word)
          if (id < 0) then
             call ferror('critic2','unknown field in REFERENCE',faterr,line,syntax=.true.)
             cycle
          end if
          if (.not.sy%goodfield(id)) then
             call ferror('critic2','REFERENCE: field is not allocated',faterr,line,syntax=.true.)
             cycle
          end if
          call sy%set_reference(id,.false.)
          write (uout,'("* Field number ",A," is now REFERENCE."/)') string(id)
          call sy%report(.false.,.false.,.true.,.false.,.false.,.false.,.true.)
          call check_no_extra_word(ok)

          ! point
       elseif (equal(word,'point')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call rhoplot_point(subline)

          ! line
       elseif (equal(word,'line')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call rhoplot_line(subline)

          ! plane
       elseif (equal(word,'plane')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call rhoplot_plane(subline)

          ! cube
       elseif (equal(word,'cube')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call rhoplot_cube(subline)

          ! grdvec
       elseif (equal(word,'grdvec')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call rhoplot_grdvec()

          ! nciplot
       elseif (equal(word,'nciplot')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call nciplot()

          ! benchmark
       elseif (equal(word,'benchmark')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          ok = eval_next(nn,line,lp)
          if (.not. ok) nn = 10000
          call check_no_extra_word(ok)
          if (.not.ok) cycle
          call sy%f(sy%iref)%benchmark(nn)

          ! auto
       elseif (equal(word,'auto')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call autocritic(subline)

          ! cpreport
       elseif (equal(word,'cpreport')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call cpreport(subline)

          ! fluxprint
       elseif (equal(word,'fluxprint')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call fluxprint()

          ! integrable
       elseif (equal(word,'integrable')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call sy%new_integrable_string(subline,errmsg)
          if (len_trim(errmsg) == 0) then
             call sy%report(.false.,.false.,.true.,.false.,.false.,.false.,.false.)
          else
             call ferror('integrable',errmsg,faterr,line,syntax=.true.)
          end if

          ! pointprop
       elseif (equal(word,'pointprop')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call sy%new_pointprop_string(subline,errmsg)
          if (len_trim(errmsg) == 0) then
             call sy%report(.false.,.false.,.false.,.true.,.false.,.false.,.false.)
          else
             call ferror('pointprop',errmsg,faterr,line,syntax=.true.)
          end if

          ! basinplot
       elseif (equal(word,'basinplot')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call basinplot(subline)

          ! bundleplot
       elseif (equal(word,'bundleplot')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call bundleplot(subline)

          ! sphereintegrals
       elseif (equal(word,'sphereintegrals')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call sphereintegrals(subline)

          ! integrals
       elseif (equal(word,'integrals')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call integrals(subline)

          ! qtree
       elseif (equal(word,'qtree')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call qtree_driver(subline)

          ! yt/bader/hirshfeld/voronoi/isosurface
       elseif (equal(word,'yt').or.equal(word,'bader').or.&
          equal(word,'hirshfeld').or.equal(word,'voronoi').or.&
          equal(word,'isosurface')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle

          if (sy%f(sy%iref)%type == type_grid) then
             call intgrid_driver(line)
          elseif (equal(word,'hirshfeld')) then
             call hirsh_nogrid()
          else
             call ferror("critic2",word // " can only be used with grids",faterr,line,syntax=.true.)
             cycle
          end if

          ! xdm
       elseif (equal(word,'xdm')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call xdm_driver(subline)

          ! stm
       elseif (equal(word,'stm')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call stm_driver(line(lp:))

          ! molcalc
       elseif (equal(word,'molcalc')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call molcalc_driver(line(lp:))

          ! root
       elseif (equal(word,'root')) then
          if (len_trim(line(lp:)) < 1) then
             call ferror('critic2','need a string for root',faterr,line,syntax=.true.)
             cycle
          end if
          fileroot = line(lp:)

          ! ewald
       elseif (equal(word,'ewald')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call check_no_extra_word(ok)
          if (.not.ok) cycle
          if (sy%c%ismolecule) then
             call ferror("critic2","EWALD can not be used with molecules",faterr)
             cycle
          end if

          rdum = sy%c%ewald_energy()
          write (uout,'("* Ewald electrostatic energy (Hartree) = ",A/)') &
             string(rdum,'e',decimal=12)

          ! environ
       elseif (equal(word,'environ')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_environ(sy,line(lp:))

          ! econ
       elseif (equal(word,'econ')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_econ(sy)

          ! coord
       elseif (equal(word,'coord')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_coord(sy,line(lp:))

          ! polyhedra
       elseif (equal(word,'polyhedra')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_polyhedra(sy,line(lp:))

          ! packing
       elseif (equal(word,'packing')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_packing(sy,line(lp:))

          ! sigmahole
       elseif (equal(word,'sigmahole')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call sigmahole_driver(sy,line(lp:))

          ! vdw
       elseif (equal(word,'vdw')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_vdw(sy,line(lp:))

          ! identify
       elseif (equal(word,'identify')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_identify(sy,line,lp)

          ! makemols_neighcrys
       elseif (equal(word,'makemolsnc')) then
          call struct_makemols_neighcrys(line,lp)

          ! molreorder
       elseif (equal(word,'molreorder')) then
          call struct_molreorder(line,lp)

          ! molmove
       elseif (equal(word,'molmove')) then
          call struct_molmove(line,lp)

          ! kpoints
       elseif (equal(word,'kpoints')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_kpoints(sy,line(lp:))

          ! bz
       elseif (equal(word,'bz')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call struct_bz(sy)

          ! sum/min/max/mean/count
       elseif (equal(word,'sum').or.equal(word,'min').or.equal(word,'max').or.&
          equal(word,'mean').or.equal(word,'count')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          word2 = getword(line,lp)
          id = sy%fieldname_to_idx(word2)
          if (id < 0) id = sy%iref
          if (.not.sy%goodfield(id)) then
             call ferror('critic2','field is not allocated in sum/min/max/mean/count',faterr,line,syntax=.true.)
             cycle
          end if
          if (.not.sy%f(id)%type == type_grid) then
             call ferror('critic2','field is not a grid',faterr,line,syntax=.true.)
             cycle
          end if
          if (.not.sy%f(id)%grid%isinit.or..not.allocated(sy%f(id)%grid%f)) then
             call ferror('critic2','grid is not initialized',faterr,line,syntax=.true.)
             cycle
          end if
          if (equal(word,"sum")) then
             write (uout,'("SUM(",A,") = ",A)') string(id), string(sum(sy%f(id)%grid%f),'e',decimal=12)
             write (uout,'("INTEGRAL(",A,") = ",A/)') string(id), &
                string(sum(sy%f(id)%grid%f * sy%f(id)%c%omega / real(product(sy%f(id)%grid%n),8)),'e',decimal=12)
          elseif(equal(word,"min")) then
             write (uout,'("MIN(",A,") = ",A/)') string(id), string(minval(sy%f(id)%grid%f),'e',decimal=12)
          elseif(equal(word,"max")) then
             write (uout,'("MAX(",A,") = ",A/)') string(id), string(maxval(sy%f(id)%grid%f),'e',decimal=12)
          elseif(equal(word,"mean")) then
             write (uout,'("MEAN(",A,") = ",A/)') string(id),&
                string(sum(sy%f(id)%grid%f) / (sy%f(id)%grid%n(1)*sy%f(id)%grid%n(2)*sy%f(id)%grid%n(3)),'e',decimal=12)
          elseif(equal(word,"count")) then
             ok = eval_next(rdum,line,lp)
             if (.not.ok) rdum = 0d0
             write (uout,'("COUNT(",A," > ",A,") = ",A/)') string(id), &
                string(rdum,'e',decimal=7), string(count(sy%f(id)%grid%f > rdum))
          end if

          ! testrmt
       elseif (equal(word,'testrmt')) then
          call check_structure_defined(ok)
          if (.not.ok) cycle
          call check_no_extra_word(ok)
          if (.not.ok) cycle
          call sy%f(sy%iref)%testrmt(3,errmsg)
          if (len_trim(errmsg) > 0) &
             call ferror('testrmt',errmsg,faterr,line,syntax=.true.)

          ! clear
       elseif (equal(word,'clear')) then
          call critic_clearvariable(line(lp:))

          ! list
       elseif (equal(word,'list')) then
          call check_no_extra_word(ok)
          if (.not.ok) cycle
          call listvariables()
          call sy%report(.false.,.true.,.true.,.true.,.true.,.false.,.false.)

          ! libxc
       elseif (equal(word,'libxc')) then
#ifdef HAVE_LIBXC
          doref = .false.
          doname = .false.
          doflags = .false.
          do while (.true.)
             word = getword(line,lp)
             if (equal(word,'ref').or.equal(word,'refs')) then
                doref = .true.
             elseif (equal(word,'name').or.equal(word,'names')) then
                doname = .true.
             elseif (equal(word,'flags')) then
                doflags = .true.
             elseif (equal(word,'all')) then
                doflags = .true.
                doname = .true.
                doref = .true.
             elseif (len_trim(word) > 0) then
                call ferror('critic2','Unknown keyword in LIBXC',faterr,line,syntax=.true.)
                cycle main
             else
                exit
             end if
          end do
          call listlibxc(doref,doname,doflags)
#else
          call ferror('critic2','critic2 was not compiled with LIBXC support',faterr,line,syntax=.true.)
          cycle
#endif

          ! spg
       elseif (equal(word,'spg')) then
          call spg_list_spg()

          ! reset
       elseif (equal(word,'reset')) then
          if (usegui) then
             write (uout,'("!! The RESET keyword cannot be used in the graphical interface.")')
             cycle main
          end if
          call check_no_extra_word(ok)
          if (.not.ok) cycle
          call critic_clearvariable("all")
          call global_set_defaults()
          call sy%end()
          call sy%init()

          ! run/system
       elseif (equal(word,'run') .or. equal(word,'system')) then
          call system(line(lp:))

          ! echo
       elseif (equal(word,'echo')) then
          write (uout,'(A)') string(line(lp:))

          ! trick
       elseif (equal(word,'trick')) then
          call trick(line(lp:))

          ! end
       elseif (equal(word,'end').or.equal(word,'exit')) then
          call check_no_extra_word(ok)
          if (.not.ok) cycle
          exit
          ! pass unknown to setvariables
       else
          lp = 1
          call critic_setvariables(line, lp)
       endif
    enddo main

  contains

    subroutine check_no_extra_word(ok)
      character(len=:), allocatable :: aux2
      logical, intent(out) :: ok
      aux2 = getword(line,lp)
      ok = .true.
      if (len_trim(aux2) > 0) then
         call ferror('critic','Unknown extra keyword',faterr,line,syntax=.true.)
         ok = .false.
      end if
    end subroutine check_no_extra_word

    subroutine check_structure_defined(ok,silent)
      logical, intent(out) :: ok
      logical, intent(in), optional :: silent

      ok = associated(sy)
      if (.not.ok) goto 999
      ok = allocated(sy%c)
      if (.not.ok) goto 999
      ok = sy%c%isinit
      if (.not.ok) goto 999

      return
999   continue
      if (present(silent)) then
         if (silent) return
      end if
      call ferror('critic2','need CRYSTAL/MOLECULE before using this keyword',faterr,line,syntax=.true.)

    end subroutine check_structure_defined

  end subroutine critic_main

  !> Initialize basic variables at the beginning of the run.
  !> Also sets default values by calling global_set_defaults.
  module subroutine global_init(ghome,datadir)
    use tools_io, only: string, ferror, warning, uout
    use param, only: dirsep
    character*(*) :: ghome, datadir
    integer :: isenv
    logical :: lchk
    character(len=:), allocatable :: wfcstr, msgr1, msgr2, msg1, msg2, msg3
    integer, parameter :: maxlenpath = 1024

    wfcstr = dirsep // "wfc" // dirsep // "h__pbe.wfc"

    ! read the -r option
    msgr1 = ""
    msgr2 = ""
    if (len_trim(ghome) > 0) then
       critic_home = string(ghome)
       inquire(file=trim(critic_home) // wfcstr,exist=lchk)
       if (lchk) goto 99
       msgr1 = "(!) 0. Not found (-r option): " // trim(critic_home) // wfcstr

       critic_home = string(ghome) // dirsep // "dat"
       inquire(file=trim(critic_home) // wfcstr,exist=lchk)
       if (lchk) goto 99
       msgr2 = "(!) 0. Not found (-r option): " // trim(critic_home) // wfcstr
    endif

    ! read env variable CRITIC_HOME
    if (allocated(critic_home)) deallocate(critic_home)
    allocate(character(len=maxlenpath)::critic_home)
    call get_environment_variable("CRITIC_HOME",critic_home,status=isenv)
    if (isenv ==0) then
       critic_home = trim(critic_home) // dirsep // "dat"
       inquire(file=trim(critic_home) // wfcstr,exist=lchk)
       if (lchk) goto 99
       msg1 = "(!) 1. Not found (CRITIC_HOME): " // trim(critic_home) // wfcstr
    else
       msg1 = "(!) 1. CRITIC_HOME environment variable not set"
    end if

    ! then the install path
    critic_home = trim(adjustl(datadir))
    inquire(file=trim(critic_home) // wfcstr,exist=lchk)
    if (lchk) goto 99
    msg2 = "(!) 2. Not found (install path): " // trim(critic_home) // wfcstr

    ! then the current directory
    critic_home = "."
    inquire(file=trim(critic_home) // wfcstr,exist=lchk)
    if (lchk) goto 99
    msg3 = "(!) 3. Not found (pwd): " // trim(critic_home) // wfcstr

    ! argh!
    call ferror("grda_init","Could not find data files.",warning)
    if (len_trim(msgr1) > 0 .and. len_trim(msgr2) > 0) &
       write (uout,'(A/A)') msgr1, msgr2
    write (uout,'(A/A/A)') msg1, msg2, msg3
    write (uout,'("(!) The density files and the structure library will not be available.")')
    critic_home = "."

99  continue

    ! library file path
    clib_file = trim(critic_home) // dirsep // "lib" // dirsep // "crystal.dat"
    mlib_file = trim(critic_home) // dirsep // "lib" // dirsep // "molecule.dat"

    ! set all default values
    call global_set_defaults()

  end subroutine global_init

  !> Set the default values for all the global variables
  module subroutine global_set_defaults()
    use meshmod, only: mesh_type_franchini, mesh_level_good

    ! global flags
    precisecube = .true.

    ! bond factor
    bondfactor = bondfactor_def

    ! symmetry
    doguess = -1
    symprec = 1d-2

    ! units
    iunit = iunit_bohr
    iunit_isdef = .true.

    ! navigation
    NAV_stepper = NAV_stepper_bs
    NAV_step = 0.3d0
    NAV_maxerr = 1d-4
    NAV_gradeps = 1d-7
    prunedist = 0.1d0

    ! integration
    INT_radquad_type = INT_gauleg
    INT_radquad_nr = 50
    INT_radquad_abserr = 1d0
    INT_radquad_relerr = 1d-12
    INT_radquad_errprop = 3
    INT_radquad_errprop_default = .true.
    INT_iasprec = 1d-5

    ! molecular mesh
    mesh_type = mesh_type_franchini
    mesh_level = mesh_level_good

  end subroutine global_set_defaults

  !> Print out the initial critic2 banner
  module subroutine initial_banner()
    use tools_io, only: uout

    write (uout,'("                  _   _     _          ___     ")')
    write (uout,'("                 (_) | |   (_)        |__ \    ")')
    write (uout,'("     ___   _ __   _  | |_   _    ___     ) |   ")')
    write (uout,'("    / __| | ''__| | | | __| | |  / __|   / /    ")')
    write (uout,'("   | (__  | |    | | | |_  | | | (__   / /_    ")')
    write (uout,'("    \___| |_|    |_|  \__| |_|  \___| |____|   ")')
    write (uout,*)
    write (uout,'("* CRITIC2: analysis of real-space scalar fields in solids and molecules.")')
    write (uout,'("  (c) 1996-2019 A. Otero-de-la-Roza, A. Martin-Pendas, V. Lua~na")')
    write (uout,'("  Distributed under GNU GPL v.3 (see COPYING for details)")')
    write (uout,'("  Bugs, requests, and rants: aoterodelaroza@gmail.com ")')
    write (uout,'("  Website: https://aoterodelaroza.github.io/critic2/")')
    write (uout,'("  If you find this software useful, please cite:")')
    write (uout,'("  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 185 (2014) 1007-1018.")')
    write (uout,'("  A. Otero-de-la-Roza et al., Comput. Phys. Commun. 180 (2009) 157-166.")')
    write (uout,*)

  end subroutine initial_banner

  !> Print out the help message
  module subroutine help_me()
    use tools_io, only: uout

    write (uout,'("Critic2 requires a sequence of keywords to operate. Please read the manual")')
    write (uout,'("distributed with the program (https://aoterodelaroza.github.io/critic2/).")')
    write (uout,'("The command-line syntax is:")')
    write (uout,'("")')
    write (uout,'("     critic2 [-q] [-h] [-r /path/to/critic2] [-l] [inputfile [outputfile]] ")')
    write (uout,'("")')
    write (uout,'("If the input or the output file are not present, stdin and stdout are used ")')
    write (uout,'("instead. The additonal options are:")')
    write (uout,'("")')
    write (uout,'("   -h: print this message.")')
    write (uout,'("   -r: have critic2 find its data in /path/to/critic2 or /path/to/critic2/dat ")')
    write (uout,'("   -l: critic2 always generates files in the current working directory")')
    write (uout,'("   -q: quiet mode.")')
    write (uout,'("   -t: testing mode.")')
    write (uout,'("")')

  end subroutine help_me

  !> Print out the compilation details and the hardwired paths
  module subroutine config_write()
    use spglib, only: spg_get_major_version, spg_get_minor_version, spg_get_micro_version
#ifdef HAVE_LIBXC
    use xc_f90_lib_m
#endif
    use config, only: getstring, istring_package, istring_atarget,&
       istring_adate, istring_datadir, istring_version
    use param, only: dirsep
    use tools_io, only: uout, string

    logical :: lchk
    integer :: iver(3)

    write (uout,'("+ ",A," (development), version ",A,"")') getstring(istring_package), getstring(istring_version)
    write (uout,'("         host: ",A)') getstring(istring_atarget)
    write (uout,'("         date: ",A)') getstring(istring_adate)
    write (uout,'(" compiled dat: ",A)') getstring(istring_datadir)
    write (uout,'("      datadir: ",A)') trim(critic_home)
    inquire(file=trim(critic_home) // dirsep // "wfc" // dirsep // "h__pbe.wfc",exist=lchk)
    write (uout,'("...was found?: ",L)') lchk

    iver(1) = spg_get_major_version()
    iver(2) = spg_get_minor_version()
    iver(3) = spg_get_micro_version()
    write (uout,'("       spglib: ",A,".",A,".",A)') string(iver(1)), string(iver(2)), string(iver(3))

#ifdef HAVE_LIBXC
    call xc_f90_version(iver(1),iver(2),iver(3))
    write (uout,'("        libxc: ",A,".",A,".",A)') string(iver(1)), string(iver(2)), string(iver(3))
#else
    write (uout,'("        libxc: <unavailable>")')
#endif
    write (uout,*)

  end subroutine config_write

  !> Parse the command line and set a global variable
  module subroutine critic_setvariables(line,lp)
    use meshmod, only: mesh_type_becke, mesh_type_franchini, mesh_level_small,&
       mesh_level_normal, mesh_level_good, mesh_level_vgood, mesh_level_amazing
    use arithmetic, only: eval, setvariable
    use tools_io, only: lgetword, getword, equal, isinteger, isreal, ferror, &
       faterr, string, uout, isassignment, getword, zatguess
    use param, only: maxzat0, atmcov, atmvdw
    use systemmod, only: sy ! xxxx
    use param, only: icrd_crys, icrd_rcrys, icrd_cart ! xxxx
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    character(len=:), allocatable :: word, var
    logical :: ok
    real*8 :: rdum
    integer :: lp2, iz
    logical :: iscov
    character(len=:), allocatable :: errmsg
    real*8 :: x(3), xx(3), dmax ! xxxx
    integer :: j, nat, ierr, lvec(3), nat1, nat2 ! xxxx
    integer, allocatable :: eid(:), lvec2(:,:) ! xxxx
    real*8, allocatable :: dist(:), up2dsp(:,:), up2dcidx(:) ! xxxx

    word = lgetword(line,lp)
    if (equal(word,'bondfactor')) then
       ok = isreal(bondfactor,line,lp)
       if (.not.ok) &
          call ferror('critic_setvariables','Wrong bondfactor',faterr,line,syntax=.true.)
       bondfactor = min(bondfactor,2.0d0)
       call check_no_extra_word(ok)
    else if (equal(word,'ode_mode')) then
       do while (.true.)
          word = lgetword(line,lp)
          if (equal(word,'method')) then
             word = lgetword(line,lp)
             if (equal(word,'euler')) then
                NAV_stepper = NAV_stepper_euler
             elseif (equal(word,'heun')) then
                NAV_stepper = NAV_stepper_heun
             else if (equal(word,'bs')) then
                NAV_stepper = NAV_stepper_bs
             else if (equal(word,'rkck')) then
                NAV_stepper = NAV_stepper_rkck
             else if (equal(word,'dp')) then
                NAV_stepper = NAV_stepper_dp
             else
                call ferror('critic_setvariables','Wrong ODE_MODE stepper',faterr,line,syntax=.true.)
             end if
          else if (equal(word,'maxstep')) then
             ok = isreal(NAV_step,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong ODE_MODE/MAXSTEP',faterr,line,syntax=.true.)
             NAV_step = NAV_step / dunit0(iunit)
          else if (equal(word,'maxerr')) then
             ok = isreal(NAV_maxerr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong ODE_MODE/MAXERR',faterr,line,syntax=.true.)
          else if (equal(word,'gradeps')) then
             ok = isreal(NAV_gradeps,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong ODE_MODE/GRADEPS',faterr,line,syntax=.true.)
          elseif (len_trim(word) > 0) then
             call ferror('critic_setvariables','Unknown keyword in ODE_MODE',faterr,line,syntax=.true.)
          else
             exit
          end if
       end do
    else if (equal(word,'prune_distance')) then
       ok = isreal(prunedist,line,lp)
       if (.not.ok) call ferror('critic_setvariables','Wrong PRUNE_DISTANCE',faterr,line,syntax=.true.)
       prunedist = prunedist / dunit0(iunit)
    else if (equal (word,'int_radial')) then
       do while(.true.)
          word = lgetword(line,lp)
          if (equal(word,'type')) then
             word = lgetword(line,lp)
             if (equal(word,'gauleg')) then
                INT_radquad_type = INT_gauleg
             else if (equal(word,'qags')) then
                INT_radquad_type = INT_qags
             else if (equal(word,'qng')) then
                INT_radquad_type = INT_qng
             else if (equal(word,'qag')) then
                INT_radquad_type = INT_qag
             else
                call ferror('critic_setvariables','Wrong INT_RADIAL type',faterr,line,syntax=.true.)
             end if
          elseif (equal(word,'nr')) then
             ok = isinteger(INT_radquad_nr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL nr',faterr,line,syntax=.true.)
          elseif (equal(word,'abserr')) then
             ok = isreal(INT_radquad_abserr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL abserr',faterr,line,syntax=.true.)
          elseif (equal(word,'relerr')) then
             ok = isreal(INT_radquad_relerr,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL relerr',faterr,line,syntax=.true.)
          elseif (equal(word,'errprop')) then
             ok = isinteger(INT_radquad_errprop,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL errprop',faterr,line,syntax=.true.)
             INT_radquad_errprop_default = .false.
          elseif (equal(word,'prec')) then
             ok = isreal(INT_iasprec,line,lp)
             if (.not.ok) call ferror('critic_setvariables','Wrong INT_RADIAL prec',faterr,line,syntax=.true.)
             INT_iasprec = INT_iasprec / dunit0(iunit)
          elseif (len_trim(word) > 0) then
             call ferror('critic_setvariables','Unknown keyword in INT_RADIAL',faterr,line,syntax=.true.)
          else
             exit
          end if
       end do
    else if (equal (word,'meshtype')) then
       ! meshtype {becke|franchini} [small|normal|good|verygood|amazing]
       word = lgetword(line,lp)
       if (equal(word,'becke')) then
          mesh_type = mesh_type_becke
       elseif (equal(word,'franchini')) then
          mesh_type = mesh_type_franchini
       else
          call ferror('critic_setvariables','Unknown keyword in MESHTYPE',faterr,line,syntax=.true.)
       end if

       lp2 = lp
       word = lgetword(line,lp)
       if (equal(word,'small')) then
          MESH_level = mesh_level_small
       else if (equal(word,'normal')) then
          MESH_level = mesh_level_normal
       else if (equal(word,'good')) then
          MESH_level = mesh_level_good
       else if (equal(word,'verygood')) then
          MESH_level = mesh_level_vgood
       else if (equal(word,'amazing')) then
          MESH_level = mesh_level_amazing
       else
          lp = lp2
       end if
       call check_no_extra_word(ok)
    elseif (equal(word,'library')) then
       word = lgetword(line,lp)
       if (equal(word,'crystal')) then
          clib_file = string(line(lp:))
       elseif (equal(word,'molecule')) then
          mlib_file = string(line(lp:))
       else
          call ferror('critic_setvariables','Unknown keyword in LIBRARY',faterr,line,syntax=.true.)
       end if
    elseif (equal(word,'units')) then
       do while(.true.)
          word = lgetword(line,lp)
          if (equal(word,'bohr').or.equal(word,'au').or.equal(word,'a.u.')) then
             iunit = iunit_bohr
             iunit_isdef = .false.
          elseif (equal(word,'ang').or.equal(word,'angstrom')) then
             iunit = iunit_ang
             iunit_isdef = .false.
          elseif (len_trim(word) > 0) then
             call ferror('critic_setvariables','Unknown keyword in UNITS',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do
    elseif (equal(word,'standardcube')) then
       precisecube = .false.
       call check_no_extra_word(ok)
    elseif (equal(word,'precisecube')) then
       precisecube = .true.
       call check_no_extra_word(ok)
    elseif (equal(word,'radii')) then
       word = lgetword(line,lp)
       if (equal(word,'cov')) then
          iscov = .true.
       elseif (equal(word,'vdw')) then
          iscov = .false.
       else
          call list_radii()
          return
       end if

       do while(.true.)
          lp2 = lp
          iz = 0
          ok = isinteger(iz,line,lp)
          if (.not.ok) then
             word = lgetword(line,lp)
             if (len_trim(word) == 0) then
                lp = lp2
                exit
             else
                iz = zatguess(word)
             end if
          end if

          if (iz < 1 .or. iz > maxzat0) then
             call ferror('critic2','Syntax error or wrong expression in RADII',faterr,line,syntax=.true.)
             return
          end if
          if (eval_next_real(rdum,line,lp)) then
             rdum = rdum / dunit0(iunit)
          else
             call ferror('critic2','Syntax error or wrong expression in RADII',faterr,line,syntax=.true.)
             return
          end if

          if (iscov) then
             atmcov(iz) = rdum
          else
             atmvdw(iz) = rdum
          end if
       end do
    elseif (equal(word,'temp')) then ! xxxx
       call sy%c%build_env()

       dmax = 7d0
       x = (/0.5d0,0.5d0,0.51d0/)
       ! x = 0.5d0
       ! x = (/0.9d0,1.1d0,0.3d0/)
       ! x=(/0.15191d0,0.20933d0,0.00913d0/)

       call sy%c%env%list_near_atoms(x,icrd_crys,.true.,nat,ierr,eid,dist,lvec,up2n=1)
       do j = 1, nat
          xx = sy%c%env%at(eid(j))%x + sy%c%env%x2xr(real(lvec,8))
          write (*,*) eid(j), xx, dist(j)
       end do
       write (*,*) "xx1 ", nat
       write (*,*)
       nat1 = nat

       call sy%c%list_near_atoms(x,icrd_crys,.true.,nat,eid,dist,lvec2,up2n=1)
       do j = 1, nat
          write (*,*) eid(j), sy%c%x2xr(sy%c%atcel(eid(j))%x + lvec2(:,j)), dist(j)
          ! write (*,*) eid(j), (sy%c%atcel(eid(j))%r + sy%c%molx0) * 0.529177d0, dist(j)
       end do
       write (*,*) "xx2 ", nat
       nat2 = nat
       ! write (*,*) "xx ", nat1, nat2

    elseif (isassignment(var,word,line)) then
       rdum = eval(word,errmsg)
       if (len_trim(errmsg) > 0) then
          call ferror('critic2','Syntax error or wrong expression',faterr,line,syntax=.true.)
          return
       end if
       call setvariable(trim(var),rdum)
    else
       word = string(line)
       rdum = eval(word,errmsg)
       if (len_trim(errmsg) > 0) then
          call ferror('critic2','Syntax error or wrong expression',faterr,line,syntax=.true.)
          return
       end if
       if (.not.quiet) then
          write (uout,'(1p,G22.14/)') rdum
       else
          write (uout,'(1p,G22.14)') rdum
       endif
    end if

  contains
    subroutine check_no_extra_word(ok)
      character(len=:), allocatable :: aux2
      logical :: ok
      aux2 = getword(line,lp)
      ok = .true.
      if (len_trim(aux2) > 0) then
         call ferror('critic_setvariables','Unknown extra keyword',faterr,line,syntax=.true.)
         ok = .false.
      end if
    end subroutine check_no_extra_word
  end subroutine critic_setvariables

  !> Clear the value of a variable
  module subroutine critic_clearvariable(line)
    use arithmetic, only: clearallvariables, clearvariable
    use tools_io, only: getword, lower, equal
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word
    integer :: lp

    lp = 1
    do while(.true.)
       word = getword(line,lp)
       if (equal(lower(word),"all")) then
          call clearallvariables()
       elseif (len_trim(word) > 0) then
          call clearvariable(word)
       else
          exit
       endif
    end do

  end subroutine critic_clearvariable

  !> Evaluate next expression of word in line and return a real number.
  module function eval_next_real(res,line,lp0)
    use arithmetic, only: eval
    use tools_io, only: isexpression_or_word, string

    logical :: eval_next_real
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    real*8, intent(out) :: res

    integer :: lp
    character(len=:), allocatable :: word
    character(len=:), allocatable :: errmsg

    res = 0d0
    lp = lp0
    word = ""
    eval_next_real = isexpression_or_word(word,line,lp)
    if (eval_next_real) then
       res = eval(string(word),errmsg)
       if (len_trim(errmsg) > 0) eval_next_real = .false.
       if (eval_next_real) lp0 = lp
    endif

  end function eval_next_real

  !> Evaluate next expression of word in line and return an integer.
  module function eval_next_int(res,line,lp0)
    use arithmetic
    use tools_io

    logical :: eval_next_int
    character*(*), intent(in) :: line !< Input line
    integer, intent(inout) :: lp0 !< Pointer to position on input line, updated after reading.
    integer, intent(out) :: res

    character(len=:), allocatable :: word
    integer :: lp
    real*8 :: rdum
    character(len=:), allocatable :: errmsg

    real*8, parameter :: eps = 1d-20

    res = 0
    lp = lp0
    word = ""
    eval_next_int = isexpression_or_word(word,line,lp)
    if (eval_next_int) then
       rdum = eval(word,errmsg)
       if (len_trim(errmsg) > 0) eval_next_int = .false.
       if (abs(rdum - nint(rdum)) > eps) then
          eval_next_int = .false.
          return
       endif
       if (eval_next_int) then
          lp0 = lp
          res = nint(rdum)
       endif
    endif

  end function eval_next_int

  !> Write to standard output the list of atomic radii
  module subroutine list_radii()
    use tools_io, only: uout, string, nameguess
    use global, only: dunit0, iunit, iunitname0
    use param, only: atmcov, atmvdw
    integer :: i
    character*(2) :: name

    write (uout,'("* List of atomic radii (per atomic number)")')
    write (uout,'("# All radii in ",A)') iunitname0(iunit)
    write (uout,'("# Z at  rcov  rvdw")')
    do i = 1, maxzat0
       name = nameguess(i,.true.)
       write (uout,'(999(" ",A))') string(i,length=3), name, &
          string(atmcov(i) * dunit0(iunit),'f',length=5,decimal=2), &
          string(atmvdw(i) * dunit0(iunit),'f',length=5,decimal=2)
    end do
    write (uout,*)

  end subroutine list_radii

end submodule proc
