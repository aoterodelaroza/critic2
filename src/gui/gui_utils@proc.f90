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

! Some utilities for building the GUI (e.g. wrappers around ImGui routines).
submodule (gui_utils) proc
  use iso_c_binding
  implicit none

contains

  !> Create a wrapped tooltip, maybe with a delay to show. ttshown
  !> activates the delay, and is the show flag for the delayed tooltip.
  module subroutine wrapped_tooltip(str,ttshown)
    use gui_interfaces_cimgui
    use gui_main, only: tooltip_wrap_factor, tooltip_delay
    character(len=*,kind=c_char), intent(in) :: str
    logical, intent(inout), optional :: ttshown

    character(len=:,kind=c_char), allocatable, target :: strloc

    if (present(ttshown)) then
       if (igIsItemHovered_delayed(ImGuiHoveredFlags_None,tooltip_delay,ttshown)) call show_tooltip()
    else
       if (igIsItemHovered(ImGuiHoveredFlags_None)) call show_tooltip()
    end if

  contains
    subroutine show_tooltip()
      strloc = trim(str) // c_null_char
      call igBeginTooltip()
      call igPushTextWrapPos(tooltip_wrap_factor * igGetFontSize())
      call igTextWrapped(c_loc(strloc))
      call igPopTextWrapPos()
      call igEndTooltip()
      ! call igSetTooltip(c_loc(strloc))
    end subroutine show_tooltip
  end subroutine wrapped_tooltip

  ! Returns true if the last item has been hovered for at least thr
  ! seconds. If already_shown (the tooltip has already been displayed),
  ! do not use the delay.
  module function igIsItemHovered_delayed(flags,thr,already_shown)
    use gui_main, only: g
    use gui_interfaces_cimgui, only: igIsItemHovered
    integer(c_int), value :: flags
    real(c_float), intent(in) :: thr
       logical, intent(in) :: already_shown
    logical(c_bool) :: igIsItemHovered_delayed

    igIsItemHovered_delayed = igIsItemHovered(flags) .and. (already_shown .or. g%HoveredIdTimer >= thr)

  end function igIsItemHovered_delayed

  ! Get a date/time string in a format adequate for the GUI
  module function get_time_string() result(output)
    use tools_io, only: string
    character(len=:), allocatable :: output

    integer :: values(8)

    call date_and_time(values=values)

    output = string(values(5),2,pad0=.true.) // ":" // string(values(6),2,pad0=.true.) //&
       ":" // string(values(7),2,pad0=.true.) // ", " // string(values(1)) // "/" // &
       string(values(2)) // "/" // string(values(3))

  end function get_time_string

end submodule proc
