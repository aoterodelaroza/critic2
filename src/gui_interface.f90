! Copyright (c) 2017 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>, Robin Myhr <x@example.com>, Isaac
! Visintainer <x@example.com>, Richard Greaves <x@example.com>, Ángel
! Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
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

!> Interface for the critic2 GUI.
module gui_interface
  use iso_c_binding, only: c_ptr, c_null_ptr
  implicit none

  private

  public :: gui_initialize
  public :: draw_scene

  type(c_ptr) :: window = c_null_ptr

contains

  !> Initialize the critic2 GUI.
  subroutine gui_initialize(window_) bind(c)
    use iso_fortran_env, only: input_unit, output_unit
    use systemmod, only: systemmod_init
    use spgs, only: spgs_init
    use config, only: datadir, version, atarget, adate, f77, fflags, fc, &
       fcflags, cc, cflags, ldflags, enable_debug, package
    use global, only: global_init, config_write, initial_banner
    use tools_io, only: ioinit, ucopy, uout, start_clock, &
       tictac, interactive, uin
    use param, only: param_init

    type(c_ptr), intent(in), value :: window_

    ! save the window pointer
    window = window_
    call assert_non_null(window,"gui_initialize","window")

    ! initialize parameters
    call start_clock()
    call param_init()

    ! input/output, arguments (tools_io)
    call ioinit()
    uin = input_unit
    uout = output_unit
    interactive = .false.

    ! set default values and initialize the rest of the modules
    call global_init("",datadir)
    call spgs_init()
    call systemmod_init(1)

    ! banner and compilation info; do not copy input
    call initial_banner()
    call config_write(package,version,atarget,adate,f77,fflags,fc,&
       fcflags,cc,cflags,ldflags,enable_debug,datadir)
    call tictac('CRITIC2')
    write (uout,*)
    ucopy = -1

  end subroutine gui_initialize

  subroutine draw_scene() bind(c)
    use gui_glfw
    use gui_glu
    use gui_gl !xxxx! put the only back when done

    type(c_ptr) :: pQuadric
    integer(c_int) :: w, h

    ! call glClearColor(0.0,0.0,0.0,0.0)
    ! call glMatrixMode(GL_PROJECTION)
    ! call glLoadIdentity()
    ! call glOrtho(0.d0,1.d0,0.d0,1.d0,-1.d0,1.d0)
    ! call glClear(GL_COLOR_BUFFER_BIT)
    ! call glColor3f(1.0,1.0,1.0)
    ! call glBegin(GL_POLYGON)
    ! call glVertex3f(0.25,0.25,0.0)
    ! call glVertex3f(0.75,0.25,0.0)
    ! call glVertex3f(0.75,0.75,0.0)
    ! call glVertex3f(0.25,0.75,0.0)
    ! call glEnd()
    ! call glFlush()

    call glfwGetFramebufferSize(window,w,h)
    call glViewport(0,0,w,h)
    call glMatrixMode(GL_PROJECTION)
    call glLoadIdentity()
    call glOrtho(-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0)
    call glMatrixMode(GL_MODELVIEW)
    call glLoadIdentity()
    ! call glClearColor(1.0,1.0,1.0,0.0)
    call glClearColor(0.0,0.0,0.4,0.0)
    call glClear(GL_COLOR_BUFFER_BIT)
    ! pQuadric = gluNewQuadric()
    ! call assert_non_null(pQuadric,"draw_scene","pQuadric")
    ! call gluSphere(pQuadric,1.5d0,32,8)
    ! call glFlush()
    
  end subroutine draw_scene

  subroutine assert_non_null(ptr,routine,ptrname)
    use iso_c_binding, only: c_associated
    use tools_io, only: ferror, faterr
    character*(*), intent(in) :: routine, ptrname
    type(c_ptr), intent(in) :: ptr

    if (.not.c_associated(ptr)) &
       call ferror(trim(routine),"Unexpected NULL pointer - " // trim(ptrname),faterr)

  end subroutine assert_non_null

end module gui_interface
