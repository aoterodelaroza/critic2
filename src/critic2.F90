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
!
!                           _ _   _      ____  
!                  ___ _ __(_) |_(_) ___|___ \
!                 / __| '__| | __| |/ __| __) |
!                | (__| |  | | |_| | (__ / __/ 
!                 \___|_|  |_|\__|_|\___|_____|
!                                     
program critic
  use tricks, only: trick
  use molcalc, only: molcalc_driver
  use stm, only: stm_driver
  use xdm, only: xdm_driver
  use hirshfeld, only: hirsh_props_grid
  use qtree, only: qtree_integration, qtree_setsphfactor
  use bisect, only: basinplot, bundleplot, sphereintegrals, integrals
  use integration, only: intgrid_driver
  use flux, only: fluxprint
  use autocp, only: autocritic, cpreport
  use nci, only: nciplot
  use rhoplot, only: rhoplot_point, rhoplot_line, rhoplot_plane, rhoplot_cube,&
     rhoplot_grdvec
  use fieldmod, only: type_grid
  use struct_drivers, only: struct_crystal_input, struct_newcell, struct_molcell,&
     struct_clearsym, struct_charges, struct_atomlabel, struct_write,&
     struct_powder, struct_rdf, struct_environ, struct_packing,&
     struct_compare, struct_identify
  use systemmod, only: systemmod_init, systemmod_end, sy
  use spgs, only: spgs_init
  use global, only: fileroot, quiet, global_init, initial_banner, config_write, &
     help_me, iunit, iunit_isdef, iunit_ang, iunit_bohr, eval_next, &
     critic_clearvariable, critic_setvariables, global_set_defaults
  use arithmetic, only: listvariables
  use grid1mod, only: grid1_clean_grids
  use config, only: datadir, version, atarget, adate, f77, fflags, fc, &
     fcflags, cc, cflags, ldflags, enable_debug, package
  use tools_io, only: uout, ucopy, uin, getline, lgetword, equal, faterr,&
     ferror, getword, string, nwarns, ncomms, ioinit, stdargs, tictac, &
     start_clock, print_clock
  use param, only: param_init
  implicit none

  ! command-line arguments
  character(len=:), allocatable :: optv, ghome
  ! parsing
  integer :: lp, lpold
  character(len=:), allocatable :: subline, word, word2
  character(len=:), allocatable :: line, errmsg
  !
  integer :: level, plevel, id
  integer :: i, nn, ismoli
  logical :: ok
  real*8 :: rdum

  ! initialize parameters
  call start_clock()
  call param_init()

  ! input/output, arguments (tools_io)
  call ioinit()
  call stdargs(optv,ghome,fileroot)

  ! set default values and initialize the rest of the modules
  call global_init(ghome,datadir)
  call spgs_init()
  call systemmod_init(1)

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
     call config_write(package,version,atarget,adate,f77,fflags,fc,&
        fcflags,cc,cflags,ldflags,enable_debug,datadir)
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
        if (equal(word,'crystal')) then
           ismoli = 0
           if (iunit_isdef) iunit = iunit_bohr
        else
           ismoli = 1
           if (iunit_isdef) iunit = iunit_ang
        end if
        call struct_crystal_input(subline,ismoli,.true.,.true.,s0=sy)

        ! newcell
     elseif (equal(word,'newcell')) then
        if (.not.sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE/MOLECULE before NEWCELL',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_newcell(sy,subline)

        ! molcell
     elseif (equal(word,'molcell')) then
        if (.not.sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before MOLCELL',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_molcell(sy,subline)

        ! clearsym/clearsymm
     elseif (equal(word,'clearsym') .or. equal(word,'clearsymm')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before CLEARSYM',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_clearsym(sy) 
        call check_no_extra_word(ok)
        if (.not.ok) cycle

        ! q/qat, zpsp, nocore
     elseif (equal(word,'q') .or. equal(word,'qat') &
        .or. equal(word,'zpsp') .or. equal(word,'nocore')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before Q/QAT/ZPSP/NOCORE',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_charges(sy,line,ok)

        ! atomlabel
     elseif (equal(word,'atomlabel')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before ATOMLABEL',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_atomlabel(sy,subline)

        ! write
     elseif (equal(word,'write')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before WRITE',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_write(sy,subline)

        ! load
     elseif (equal(word,'load')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before LOAD',faterr,line,syntax=.true.)
           cycle
        end if
        call sy%load_field_string(subline,id,errmsg)
        if (id < 0 .or. len_trim(errmsg) > 0) &
           call ferror('load',errmsg,faterr,line,syntax=.true.)

        ! unload
     elseif (equal(word,'unload')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before UNLOAD',faterr,line,syntax=.true.)
           cycle
        end if
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before POWDER',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_powder(sy,line(lp:))

        ! rdf
     elseif (equal(word,'rdf')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before RDF',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_rdf(sy,line(lp:))

        ! compare
     elseif (equal(word,'compare')) then
        call struct_compare(sy,line(lp:))

        ! setfield
     elseif (equal(word,'setfield')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before SETFIELD',faterr,line,syntax=.true.)
           cycle
        end if
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before REFERENCE',faterr,line,syntax=.true.)
           cycle
        end if
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before POINT',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_point(subline)

        ! line
     elseif (equal(word,'line')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before LINE',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_line(subline)

        ! plane
     elseif (equal(word,'plane')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before PLANE',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_plane(subline)

        ! cube
     elseif (equal(word,'cube')) then
        if (.not. sy%c%isinit) then 
           call ferror('critic2','need CRYSTAL/MOLECULE before CUBE',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_cube(subline)

        ! grdvec
     elseif (equal(word,'grdvec')) then
        if (.not. sy%c%isinit) then
           call ferror('critic','need CRYSTAL/MOLECULE before GRDVEC',faterr,line,syntax=.true.)
           cycle
        end if
        call rhoplot_grdvec()

        ! nciplot
     elseif (equal(word,'nciplot')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before NCIPLOT',faterr,line,syntax=.true.)
           cycle
        end if
        call nciplot()

        ! benchmark 
     elseif (equal(word,'benchmark')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before BENCHMARK',faterr,line,syntax=.true.)
           cycle
        end if
        ok = eval_next(nn,line,lp)
        if (.not. ok) nn = 10000
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call sy%f(sy%iref)%benchmark(nn)

        ! auto
     elseif (equal(word,'auto')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before AUTO',faterr,line,syntax=.true.)
           cycle
        end if
        call autocritic(subline)

        ! cpreport
     elseif (equal(word,'cpreport')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before CPREPORT',faterr,line,syntax=.true.)
           cycle
        end if
        call cpreport(subline)

        ! fluxprint
     elseif (equal(word,'fluxprint')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before FLUXPRINT',faterr,line,syntax=.true.)
           cycle
        end if
        call fluxprint()

        ! integrable
     elseif (equal(word,'integrable')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before INTEGRABLE',faterr,line,syntax=.true.)
           cycle
        end if
        call sy%new_integrable_string(subline,errmsg)
        if (len_trim(errmsg) == 0) then
           call sy%report(.false.,.false.,.true.,.false.,.false.,.false.,.false.)
        else
           call ferror('integrable',errmsg,faterr,line,syntax=.true.)
        end if

        ! pointprop
     elseif (equal(word,'pointprop')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before POINTPROP',faterr,line,syntax=.true.)
           cycle
        end if
        call sy%new_pointprop_string(subline,errmsg)
        if (len_trim(errmsg) == 0) then
           call sy%report(.false.,.false.,.false.,.true.,.false.,.false.,.false.)
        else
           call ferror('pointprop',errmsg,faterr,line,syntax=.true.)
        end if

        ! basinplot
     elseif (equal(word,'basinplot')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before BASINPLOT',faterr,line,syntax=.true.)
           cycle
        end if
        call basinplot(subline)

        ! bundleplot
     elseif (equal(word,'bundleplot')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before BUNDLEPLOT',faterr,line,syntax=.true.)
           cycle
        end if
        call bundleplot(subline)

        ! sphereintegrals
     elseif (equal(word,'sphereintegrals')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before SPHEREINTEGRALS',faterr,line,syntax=.true.)
           cycle
        end if
        call sphereintegrals(subline)
        
        ! integrals
     elseif (equal(word,'integrals')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before INTEGRALS',faterr,line,syntax=.true.)
           cycle
        end if
        call integrals(subline)
        
        ! qtree 
     elseif (equal(word,'qtree')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before QTREE',faterr,line,syntax=.true.)
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before YT',faterr,line,syntax=.true.)
           cycle
        end if
        call intgrid_driver(line)

        ! bader
     elseif (equal(word,'bader')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before BADER',faterr,line,syntax=.true.)
           cycle
        end if
        call intgrid_driver(line)

        ! xdm
     elseif (equal(word,'xdm')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before XDM',faterr,line,syntax=.true.)
           cycle
        end if
        call xdm_driver(subline)

        ! stm
     elseif (equal(word,'stm')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before STM',faterr,line,syntax=.true.)
           cycle
        end if
        call stm_driver(line(lp:))

        ! molcalc
     elseif (equal(word,'molcalc')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need MOLECULE before MOLCALC',faterr,line,syntax=.true.)
           cycle
        end if
        if (.not.sy%c%ismolecule) then
           call ferror("critic2","MOLCALC can not be used with crystals",faterr)
           cycle
        end if
        call molcalc_driver(line(lp:))

        ! sphfactor
     elseif (equal(word,'sphfactor')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before SPHFACTOR',faterr,line,syntax=.true.)
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before HIRSHFELD',faterr,line,syntax=.true.)
           cycle
        end if
        call check_no_extra_word(ok)
        if (.not.ok) cycle
        call hirsh_props_grid()

        ! ewald
     elseif (equal(word,'ewald')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before EWALD',faterr,line,syntax=.true.)
           cycle
        end if
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before environ',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_environ(sy,line(lp:))

        ! packing
     elseif (equal(word,'packing')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before packing',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_packing(sy,line(lp:))

        ! identify
     elseif (equal(word,'identify')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before identify',faterr,line,syntax=.true.)
           cycle
        end if
        call struct_identify(sy,line,lp)

        ! sum/min/max/mean/count
     elseif (equal(word,'sum').or.equal(word,'min').or.equal(word,'max').or.&
        equal(word,'mean').or.equal(word,'count')) then
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before sum/min/max/mean/count',faterr,line,syntax=.true.)
           cycle
        end if
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
           write (uout,'("SUM(",A,") = ",A/)') string(id), string(sum(sy%f(id)%grid%f),'e',decimal=12)
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
        if (.not. sy%c%isinit) then
           call ferror('critic2','need CRYSTAL/MOLECULE before testrmt',faterr,line,syntax=.true.)
           cycle
        end if
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

        ! reset
     elseif (equal(word,'reset')) then
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

        ! temp, for testing
     elseif (equal(word,'temp')) then
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

  call grid1_clean_grids()
  ! call systemmod_end() ! older ifort compilers have trouble deallocating sy

  if (.not.quiet) then
     write (uout,'("CRITIC2 ended succesfully (",A," WARNINGS, ",A," COMMENTS)"/)')&
        string(nwarns), string(ncomms)
     call print_clock()
     call tictac('CRITIC2')
  endif

999 continue

  ! pause at the end of the windows execution so I can see the output
#ifdef WIN
  read (*,*) 
#endif

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

end program critic

