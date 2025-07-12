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

! Routines for the scfplot window.
submodule (windows) scfplot
  use interfaces_cimgui
  implicit none
contains

  !> Update tasks for the SCF window, before the window is created.
  module subroutine update_scfplot(w)
    use gui_main, only: sysc, sys_init, ok_system
    class(window), intent(inout), target :: w

    integer :: isys
    logical :: doquit

    isys = w%isys
    doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = (sysc(isys)%collapse >= 0)
    if (doquit) call w%end()

  end subroutine update_scfplot

  !> Draw the SCF plot window.
  module subroutine draw_scfplot(w)
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use gui_main, only: nsys, sysc, sys_init, ok_system
    use utils, only: iw_text
    use types, only: realloc
    use tools_io, only: string
    class(window), intent(inout), target :: w

    real(c_double) :: xmin, xmax
    integer(c_int) :: nticks
    real*8 :: dy
    logical :: doquit
    integer :: i, isys, num
    type(ImVec2) :: sz
    type(ImVec4) :: auto
    character(len=:,kind=c_char), allocatable, target :: str1, str2

    integer, parameter :: maxxticks = 10

    isys = w%isys
    doquit = .not.ok_system(isys,sys_init)
    if (.not.doquit) doquit = (sysc(isys)%collapse >= 0)
    if (w%firstpass) w%plotn = 0

    if (.not.doquit) then
       ! get the property and prepare x and y
       if (w%firstpass) then
          if (allocated(w%plotx)) deallocate(w%plotx)
          if (allocated(w%ploty)) deallocate(w%ploty)
          num = count(sysc(1:nsys)%collapse == isys) + 1
          allocate(w%plotx(num),w%ploty(num))
          w%plotn = 0
          do i = 1, nsys
             if (sysc(i)%collapse == isys .and. sysc(i)%seed%energy /= huge(1d0)) then
                w%plotn = w%plotn + 1
                w%plotx(w%plotn) = w%plotn
                w%ploty(w%plotn) = sysc(i)%seed%energy
             end if
          end do
          if (sysc(isys)%seed%energy /= huge(1d0)) then
             w%plotn = w%plotn + 1
             w%plotx(w%plotn) = w%plotn
             w%ploty(w%plotn) = sysc(isys)%seed%energy
          end if
          call realloc(w%plotx,w%plotn)
          call realloc(w%ploty,w%plotn)
          w%ymax = maxval(w%ploty)
          w%ymin = minval(w%ploty)
       end if

       ! make the plot
       if (w%plotn > 0) then
          str1 = "##scfiterations" // string(isys) // c_null_char
          call igGetContentRegionAvail(sz)
          if (ipBeginPlot(c_loc(str1),sz,ImPlotFlags_None)) then
             str1 = "SCF Iteration" // c_null_char
             str2 = "Energy (Ha)" // c_null_char
             call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_AutoFit,ImPlotAxisFlags_AutoFit)

             str1 = "%d" // c_null_char
             nticks = size(w%plotx,1)
             xmin = w%plotx(1)
             xmax = w%plotx(nticks)
             if (nticks > maxxticks) then
                nticks = maxxticks + 1
                xmax = xmin + ceiling((xmax - xmin) / maxxticks) * maxxticks
             end if
             call ipSetupAxisTicks(ImAxis_X1,xmin,xmax,nticks)
             call ipSetupAxisFormat(ImAxis_X1,c_loc(str1))

             dy = (w%ymax - w%ymin)
             str1 = "%." // string(min(ceiling(max(abs(log10(dy)),0d0)) + 1,10)) // "f" // c_null_char
             call ipSetupAxisFormat(ImAxis_Y1,c_loc(str1))
             call ipGetPlotCurrentLimits(xmin,xmax,w%ymin,w%ymax) ! these need to be switched

             str1 = "Energy" // c_null_char
             auto%x = 0._c_float
             auto%y = 0._c_float
             auto%z = 0._c_float
             auto%w = -1._c_float
             call ipSetNextMarkerStyle(ImPlotMarker_Circle,-1._c_float,auto,-1._c_float,auto)

             call ipPlotLine(c_loc(str1),c_loc(w%plotx),c_loc(w%ploty),w%plotn,ImPlotLineFlags_None,0_c_int)
             call ipEndPlot()
          end if
       end if
    end if

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    if (doquit) call w%end()

  end subroutine draw_scfplot

end submodule scfplot
