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
!                           _ _   _      ____  
!                  ___ _ __(_) |_(_) ___|___ \
!                 / __| '__| | __| |/ __| __) |
!                | (__| |  | | |_| | (__ / __/ 
!                 \___|_|  |_|\__|_|\___|_____|
!                                     
program critic
  use tricks
  use stm
  use xdm
  use ewald
  use hirshfeld
  use qtree
  use bisect
  use integration
  use flux
  use autocp
  use nci
  use rhoplot
  use fields
  use varbas
  use grd_atomic
  use struct
  use struct_basic
  use wfn_private
  use pi_private
  use spgs
  use global
  use config
  use graphics
  use arithmetic
  use types
  use tools
  use tools_io
  use param
  implicit none

  ! command-line arguments
  character(len=:), allocatable :: optv, ghome
  ! parsing
  integer :: lp, lpold
  character(len=:), allocatable :: subline, word
  character(len=:), allocatable :: line
  !
  integer :: level, plevel
  integer :: i, id, nn
  logical :: ll1, ok
  real*8 :: rdum, rdum1(3,3)

  ! initialize parameters
  call param_init()

  ! input/output, arguments (tools_io)
  call ioinit()
  call stdargs(optv,ghome,fileroot)

  ! set default values and initialize the rest of the modules
  call global_init(ghome)
  call cr%init()
  call fields_init()
  call spgs_init()
  call graphics_init()

  ! parse global control options
  quiet = (index(optv,"q") /= 0)
  if (index(optv,"h") /= 0) then
     call initial_banner()
     call help_me()
     goto 999
  endif

  ! header, interface, date
  if (.not.quiet) then
     call initial_banner()
     call config_write()
     call tictac('CRITIC2')
     write (uout,*)
     ucopy = uout
  else
     ucopy = -1
  endif

  ! Start reading
  do while (getline(uin,line,ucopy=ucopy))
     lp=1
     word = lgetword(line,lp)
     subline = line(lp:)

     ! crystal
     if (equal(word,'crystal') .or. equal(word,'molecule')) then
        ! there is a previous crystal structure, clean up...
        if (cr%isinit) call clean_structure()
        ! read the crystal enviornment
        call struct_crystal_input(cr,subline,equal(word,'molecule'),.true.,.true.)
        ! initialize the radial densities
        call grda_init(.true.,.true.,.true.)
        ! set the promolecular density as reference
        call set_reference(0)

     elseif (equal(word,'newcell')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before newcell',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_newcell(subline)
        ! clean the CP list
        call varbas_end()
        ! unload all fields and set the promolecular density as reference
        call fields_end()
        call fields_init()
        call set_reference(0)
        
     elseif (equal(word,'molcell')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before molcell',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_molcell(subline)
        
     ! clearsym/clearsymm
     elseif (equal(word,'clearsym') .or. equal(word,'clearsymm')) then
        call struct_clearsym() 
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        ! clean the CP list
        call varbas_end()
        call init_cplist(.true.)

     ! q/qat, zpsp, nocore
     elseif (equal(word,'q') .or. equal(word,'qat') &
        .or. equal(word,'zpsp') .or. equal(word,'nocore')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before q/qat/zpsp/nocore',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_charges(line)
        ll1 = equal(word,'zpsp') .or. equal(word,'nocore')
        call grda_init(ll1,.not.ll1,.true.)
           
     ! write
     elseif (equal(word,'write')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before write',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_write(subline)

     ! load
     elseif (equal(word,'load')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before load',faterr,line,syntax=.true.)
           cycle
        end if
        call fields_load(subline,id)
        if (refden == 0) call set_reference(id)

     ! unload
     elseif (equal(word,'unload')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before unload',faterr,line,syntax=.true.)
           cycle
        end if
        lpold = lp
        word = getword(line,lp)
        nn = fieldname_to_idx(word)
        if (nn >= 0) then
           if (.not.goodfield(nn)) then
              call ferror('critic2','wrong field in UNLOAD',faterr,line,syntax=.true.)
              cycle
           end if
           if (nn == 0) then
              call ferror('critic2','can not unload the promolecular density',faterr,line,syntax=.true.)
              cycle
           end if
           call fields_unload(nn)
        else
           lp = lpold
           word = lgetword(line,lp)
           if (equal(word,'all')) then
              do i = 1, ubound(f,1)
                 if (fused(i)) call fields_unload(i)
              end do
           else
              call ferror('critic2','Unknown keyword in UNLOAD',faterr,line,syntax=.true.)
           endif
        end if
           
     ! powder
     elseif (equal(word,'powder')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before powder',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_powder(line(lp:),cr)

     ! rdf
     elseif (equal(word,'rdf')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before rdf',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_rdf(line(lp:),cr)

     ! compare
     elseif (equal(word,'compare')) then
        call struct_compare(line(lp:))

     ! setfield
     elseif (equal(word,'setfield')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before setfield',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) id = refden
        if (.not.goodfield(id)) then
           call ferror('critic2','wrong field in setfield',faterr,line,syntax=.true.)
           cycle
        end if
        call setfield(f(id),id,line(lp:),.true.)

     ! reference
     elseif (equal(word,'reference')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before reference',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) then
           call ferror('critic2','unknown field in REFERENCE',faterr,line,syntax=.true.)
           cycle
        end if
        if (.not.goodfield(id)) then
           call ferror('critic2','REFERENCE: field is not allocated',faterr,line,syntax=.true.)
           cycle
        end if
        call set_reference(id)
        call check_no_extra_word(ok)

     ! point
     elseif (equal(word,'point')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before point',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_point(subline)

     ! line
     elseif (equal(word,'line')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before line',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_line(subline)

     ! plane
     elseif (equal(word,'plane')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before plane',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_plane(subline)

     ! cube
     elseif (equal(word,'cube')) then
        if (.not. cr%isinit) then 
           call ferror('critic2','need crystal before cube',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_cube(subline)

     ! grdvec
     elseif (equal(word,'grdvec')) then
        if (.not. cr%isinit) then
           call ferror('critic','need crystal before grdvec',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_grdvec()

     ! nciplot
     elseif (equal(word,'nciplot')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before nciplot',faterr,line,syntax=.true.)
           cycle
        end if
        call nciplot()

     ! benchmark 
     elseif (equal(word,'benchmark')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before benchmark',faterr,line,syntax=.true.)
           cycle
        end if
        ok = eval_next(nn,line,lp)
        if (.not. ok) nn = 10000
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call benchmark(nn)

     ! auto
     elseif (equal(word,'auto')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before auto',faterr,line,syntax=.true.)
           cycle
        end if
        call autocritic(subline)

     ! cpreport
     elseif (equal(word,'cpreport')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before cpreport',faterr,line,syntax=.true.)
           cycle
        end if
        call cpreport(subline)

     ! fluxprint
     elseif (equal(word,'fluxprint')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before fluxprint',faterr,line,syntax=.true.)
           cycle
        end if
        call fluxprint()

     ! integrable
     elseif (equal(word,'integrable')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before integrable',faterr,line,syntax=.true.)
           cycle
        end if
        call fields_integrable(subline)

     ! pointprop
     elseif (equal(word,'pointprop')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before pointprop',faterr,line,syntax=.true.)
           cycle
        end if
        call fields_pointprop(subline)

     ! basinplot
     elseif (equal(word,'basinplot')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before basinplot',faterr,line,syntax=.true.)
           cycle
        end if
        call basinplot(subline)

     ! bundleplot
     elseif (equal(word,'bundleplot')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before bundleplot',faterr,line,syntax=.true.)
           cycle
        end if
        call bundleplot(subline)

     ! sphereintegrals
     elseif (equal(word,'sphereintegrals')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before sphereintegrals',faterr,line,syntax=.true.)
           cycle
        end if
        call sphereintegrals(subline)

     ! integrals
     elseif (equal(word,'integrals')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before integrals',faterr,line,syntax=.true.)
           cycle
        end if
        call integrals(subline)

     ! qtree 
     elseif (equal(word,'qtree')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before qtree',faterr,line,syntax=.true.)
           cycle
        end if
        ok = eval_next(level,line,lp)
        if (.not.ok) level = 6
        ok = eval_next(plevel,line,lp)
        if (.not.ok) plevel = 0
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call qtree_integration(level,plevel)

     ! yt
     elseif (equal(word,'yt')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before yt',faterr,line,syntax=.true.)
           cycle
        end if
        call intgrid_driver(line)

     ! bader
     elseif (equal(word,'bader')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before bader',faterr,line,syntax=.true.)
           cycle
        end if
        call intgrid_driver(line)

     ! xdm
     elseif (equal(word,'xdm')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before xdm',faterr,line,syntax=.true.)
           cycle
        end if
        call xdm_driver(subline)

     ! stm
     elseif (equal(word,'stm')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before stm',faterr,line,syntax=.true.)
           cycle
        end if
        call stm_driver(line(lp:))

     ! sphfactor
     elseif (equal(word,'sphfactor')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before sphfactor',faterr,line,syntax=.true.)
           cycle
        end if
        call qtree_setsphfactor(subline)

     ! root
     elseif (equal(word,'root')) then
        if (len_trim(line(lp:)) < 1) then
           call ferror('critic2','need a string for root',faterr,line,syntax=.true.)
           cycle
        end if
        fileroot = line(lp:)

     ! hirshfeld
     elseif (equal(word,'hirshfeld')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before hirshfeld',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call hirsh_props_grid()

     ! ewald
     elseif (equal(word,'ewald')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before ewald',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        if (cr%ismolecule) then
           call ferror("critic2","EWALD can not be used with molecules",faterr)
           cycle
        end if

        call ewald_energy(rdum)
        write (uout,'("* Ewald electrostatic energy (Hartree) = ",A/)') &
           string(rdum,'e',decimal=12)

     ! ws
     elseif (equal(word,'ws')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before ws',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        if (cr%ismolecule) then
           call ferror("critic2","WS can not be used with molecules",faterr,syntax=.true.)
           cycle
        end if
           
        call cr%wigner((/0d0,0d0,0d0/),.true.)

     ! environ
     elseif (equal(word,'environ')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before environ',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_environ(line(lp:))

     ! packing
     elseif (equal(word,'packing')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before packing',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_packing(line(lp:))
        
     ! identify
     elseif (equal(word,'identify')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before identify',faterr,line,syntax=.true.)
           cycle
        end if
        call varbas_identify(line,lp)

     ! sum
     elseif (equal(word,'sum')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before sum',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) id = refden
        if (.not.goodfield(id)) then
           call ferror('critic2','SUM: field is not allocated',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        write (uout,'("SUM(",A,") = ",A/)') string(id), string(sum(f(id)%f),'e',decimal=12)

     ! min
     elseif (equal(word,'min')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before min',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) id = refden
        if (.not.goodfield(id)) then
           call ferror('critic2','MIN: field is not allocated',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        write (uout,'("MIN(",A,") = ",A/)') string(id), string(minval(f(id)%f),'e',decimal=12)

     ! max
     elseif (equal(word,'max')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before max',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) id = refden
        if (.not.goodfield(id)) then
           call ferror('critic2','MAX: field is not allocated',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        write (uout,'("MAX(",A,") = ",A/)') string(id), string(maxval(f(id)%f),'e',decimal=12)

     ! mean
     elseif (equal(word,'mean')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before mean',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) id = refden
        if (.not.goodfield(id)) then
           call ferror('critic2','MEAN: field is not allocated',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        write (uout,'("MEAN(",A,") = ",A/)') string(id),&
           string(sum(f(id)%f) / (f(id)%n(1)*f(id)%n(2)*f(id)%n(3)),'e',decimal=12)

     ! count
     elseif (equal(word,'count')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before count',faterr,line,syntax=.true.)
           cycle
        end if
        word = getword(line,lp)
        id = fieldname_to_idx(word)
        if (id < 0) id = refden
        if (.not.goodfield(id)) then
           call ferror('critic2','COUNT: field is not allocated',faterr,line,syntax=.true.)
           cycle
        end if

        ok = eval_next(rdum,line,lp)
        if (.not.ok) rdum = 0d0

        if (f(id)%type /= type_grid) then
           call ferror("count","count can only be used with grids",faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        write (uout,'("COUNT(",A," > ",A,") = ",A/)') string(id), &

           string(rdum,'e',decimal=7), string(count(f(id)%f > rdum))

     ! testrmt
     elseif (equal(word,'testrmt')) then
        if (.not. cr%isinit) then
           call ferror('critic2','need crystal before testrmt',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call testrmt(refden,3)

     ! clear
     elseif (equal(word,'clear')) then
        call critic_clearvariable(line(lp:))

     ! list
     elseif (equal(word,'list')) then
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call listvariables()
        call listfieldalias()

     ! reset
     elseif (equal(word,'reset')) then
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call critic_clearvariable("all")
        call global_set_defaults()
        call clean_structure()

     ! run/system
     elseif (equal(word,'run') .or. equal(word,'system')) then
        call system(line(lp:))
        
     ! echo
     elseif (equal(word,'echo')) then
        write (uout,'(A)') string(line(lp:))
        
     ! trick
     elseif (equal(word,'trick')) then
        call trick(line(lp:))
        
     ! ! temp, for testing
     ! elseif (equal(word,'temp')) then
     !    call cr%classify()
     !    
     ! end
     elseif (equal(word,'end')) then
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        exit
        ! pass unknown to setvariables
     else
        lp = 1
        call critic_setvariables(line, lp)
     endif
  enddo

  ! clean up
  call pi_end()
  call wfn_end()
  call grda_end()
  call varbas_end()
  call cr%end()
  call fields_end()

  if (.not.quiet) then
     write (uout,'("CRITIC2 ended succesfully (",A," WARNINGS, ",A," COMMENTS)"/)')&
        string(nwarns), string(ncomms)
     call tictac('CRITIC2')
  endif

999 continue

  ! pause at the end of the windows execution so I can see the output
#ifdef WIN
  read (*,*) 
#endif

contains
  !> Set field number id as reference
  subroutine set_reference(id)
    implicit none

    integer, intent(in) :: id

    ! header and change refden
    refden = id
    write (uout,'("* Field number ",A," is now REFERENCE."/)') string(refden)

    ! initialize CP list, defer the calculation of nuclei properties to the report
    call init_cplist(.true.)

    ! define second integrable property as the valence charge.
    nprops = max(2,nprops)
    integ_prop(2)%used = .true.
    integ_prop(2)%itype = itype_fval
    integ_prop(2)%fid = id
    integ_prop(2)%prop_name = "Pop"

    ! define third integrable property as the valence laplacian.
    nprops = max(3,nprops)
    integ_prop(3)%used = .true.
    integ_prop(3)%itype = itype_lapval
    integ_prop(3)%fid = id
    integ_prop(3)%prop_name = "Lap"

    ! report
    call fields_integrable_report()

    ! reset defaults for qtree
    if (f(refden)%type == type_grid) then
       gradient_mode = 1
       if (INT_radquad_errprop_default) INT_radquad_errprop = 2
    else
       gradient_mode = 2
       if (INT_radquad_errprop_default) INT_radquad_errprop = 3
    end if

  end subroutine set_reference

  !> Cleans the structure and associated data
  subroutine clean_structure()
    ! ...the structural data
    call cr%end()
    call cr%init()
    ! ...the fields associated to the previous structure
    call fields_end()
    call fields_init()
    ! ...the loaded radial atomic and core densities
    call grda_end()
    ! ...the CP list
    call varbas_end()
  end subroutine clean_structure

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

end program

