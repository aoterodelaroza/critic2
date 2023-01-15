! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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
#ifdef HAVE_GUI
  use gui_main, only: gui_start
#endif
  use spglib, only: spg_build_hall_mapping
  use systemmod, only: systemmod_init
  use grid1mod, only: grid1_clean_grids
  use global, only: fileroot, quiet, testing, initial_banner, config_write, help_me,&
     global_init, critic_main
  use config, only: getstring, istring_datadir
  use tools_io, only: ncomms, nwarns, ucopy, uout, string, start_clock, stdargs,&
     history_init, history_end, tictac, print_clock, usegui
  use param, only: param_init
  implicit none

  ! command-line arguments
  character(len=:), allocatable :: optv, ghome

  ! initialize parameters
  call start_clock()
  call param_init()

  ! input/output, arguments (tools_io)
  call stdargs(optv,ghome,fileroot)
  call history_init()

  ! set default values and initialize the rest of the modules
  call global_init(ghome,getstring(istring_datadir))
  call systemmod_init(1)
  call spg_build_hall_mapping() ! for spglib

  ! parse global control options
  testing = (index(optv,"t") /= 0)
  quiet = (index(optv,"q") /= 0 .or.usegui)
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

  if (usegui) then
#ifdef HAVE_GUI
     ! the GUI
     call gui_start()
#endif
  else
     ! the text interface
     call critic_main()
  end if

  call grid1_clean_grids()
  call history_end()

  if (.not.quiet) then
     write (uout,'("CRITIC2 ended successfully (",A," WARNINGS, ",A," COMMENTS)"/)')&
        string(nwarns), string(ncomms)
     call print_clock()
     call tictac('CRITIC2')
  endif

999 continue

    ! pause at the end of the windows execution so I can see the output
#ifdef _WIN32
    read (*,*)
#endif

end program critic
