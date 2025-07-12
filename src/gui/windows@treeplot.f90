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

! Routines for the view/edit geometry window.
submodule (windows) treeplot
  use interfaces_cimgui
  implicit none
contains

  !> Draw the tree plot window
  module subroutine draw_treeplot(w)
    use interfaces_glfw, only: glfwGetTime
    use gui_main, only: sysc, sys_empty, sys_init, sys
    use windows, only: win, iwin_tree
    use utils, only: iw_calcwidth, iw_combo_simple, iw_tooltip, iw_text
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_CLOSE_ALL_DIALOGS,&
       BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string
    use types, only: realloc
    use param, only: bohrtoa
    class(window), intent(inout), target :: w

    type(ImVec2) :: sz
    logical :: ok
    integer :: i, j, nshown
    real*8 :: valx, valy
    character(len=:,kind=c_char), allocatable, target :: str1, str2
    type(ImVec4) :: auto
    logical(c_bool) :: ch, forceupdate

    integer, save :: ic_plotx = 0, ic_ploty = 0
    logical, save :: ttshown = .false. ! tooltip flag

    integer, parameter :: ictrans(0:14) = (/ic_tree_id,ic_tree_e,ic_tree_v,ic_tree_vmol,ic_tree_nneq,ic_tree_ncel,&
       ic_tree_nmol,ic_tree_a,ic_tree_b,ic_tree_c,ic_tree_alpha,ic_tree_beta,ic_tree_gamma,ic_tree_emol,ic_tree_p/)

    ! initialize
    forceupdate = .false.
    if (w%firstpass) then
       ic_plotx = 2
       ic_ploty = 1
       forceupdate = .true.
    end if

    ! update the plot based on time signals between dependent windows
    !! if the tree info has been updated, also update the plot
    if (w%timelast_plot_update < win(iwin_tree)%timelast_tree_update) forceupdate = .true.

    ! x-choose and y-choose
    str1 = "ID"//c_null_char//"Energy (Ha)"//c_null_char//"Volume (Å³)"//c_null_char//&
       "Volume/Z (Å³)"//c_null_char//&
       "Number of symmetry-unique atoms"//c_null_char//"Number of cell atoms"//c_null_char//&
       "Number of molecules"//c_null_char//"a (Å)"//c_null_char//&
       "b (Å)"//c_null_char//"c (Å)"//c_null_char//"α (°)"//c_null_char//"β (°)"//c_null_char//&
       "γ (°)"//c_null_char//"Energy/Z (Ha)"//c_null_char//"Pressure (GPa)"//c_null_char//c_null_char
    call iw_text("x: ")
    call iw_combo_simple("##xselect",str1,ic_plotx,changed=ch,sameline=.true.)
    call iw_tooltip("Property to represent on the x-axis",ttshown)
    forceupdate = forceupdate .or. ch

    call iw_text("y: ")
    call iw_combo_simple("##yselect",str1,ic_ploty,changed=ch,sameline=.true.)
    call iw_tooltip("Property to represent on the y-axis",ttshown)
    forceupdate = forceupdate .or. ch

    ! update the plot data if necessary
    if (forceupdate) then
       w%plotn = 0
       if (allocated(w%plotx)) deallocate(w%plotx)
       if (allocated(w%ploty)) deallocate(w%ploty)
       nshown = size(win(iwin_tree)%iord,1)
       allocate(w%plotx(nshown),w%ploty(nshown))
       do j = 1, nshown
          i = win(iwin_tree)%iord(j)
          if (sysc(i)%status == sys_empty .or. sysc(i)%hidden) cycle

          ok = getvalue(valx,ictrans(ic_plotx),i)
          ok = ok .and. getvalue(valy,ictrans(ic_ploty),i)
          if (ok) then
             w%plotn = w%plotn + 1
             w%plotx(w%plotn) = valx
             w%ploty(w%plotn) = valy
          end if
       end do
       if (w%plotn > 0) then
          call realloc(w%plotx,w%plotn)
          call realloc(w%ploty,w%plotn)
       end if
       w%timelast_plot_update = glfwGetTime()
    end if

    ! make the plot
    if (w%plotn > 0) then
       str1 = "##treeplot" // c_null_char
       call igGetContentRegionAvail(sz)
       if (ipBeginPlot(c_loc(str1),sz,ImPlotFlags_None)) then
          call getname(str1,ictrans(ic_plotx))
          call getname(str2,ictrans(ic_ploty))
          call ipSetupAxes(c_loc(str1),c_loc(str2),ImPlotAxisFlags_AutoFit,ImPlotAxisFlags_AutoFit)

          auto%x = 0._c_float
          auto%y = 0._c_float
          auto%z = 0._c_float
          auto%w = -1._c_float
          call ipSetNextMarkerStyle(ImPlotMarker_Circle,-1._c_float,auto,-1._c_float,auto)

          call ipPlotLine(c_loc(str2),c_loc(w%plotx),c_loc(w%ploty),w%plotn,ImPlotLineFlags_None,0_c_int)
          call ipEndPlot()
       end if
    end if

    ! close
    if ((w%focused() .and. (is_bind_event(BIND_OK_FOCUSED_DIALOG) .or. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG) .or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)))) &
       call w%end()

  contains
    ! Get value ic for system i and return it in val. The function returns ok if it
    ! is a valid value.
    function getvalue(val,ic,i)
      real*8, intent(out) :: val
      integer, intent(in) :: ic, i
      logical :: getvalue

      val = 0d0
      getvalue = .false.
      if (ic == ic_tree_id) then
         val = real(i,8)
      elseif (ic == ic_tree_e) then
         if (sysc(i)%seed%energy == huge(1d0)) return
         val = sysc(i)%seed%energy
      elseif (ic == ic_tree_v) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%omega*bohrtoa**3
      elseif (ic == ic_tree_vmol) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%omega*bohrtoa**3/sys(i)%c%nmol
      elseif (ic == ic_tree_nneq) then
         if (sysc(i)%status /= sys_init) return
         val = sys(i)%c%nneq
      elseif (ic == ic_tree_ncel) then
         if (sysc(i)%status /= sys_init) return
         val = sys(i)%c%ncel
      elseif (ic == ic_tree_nmol) then
         if (sysc(i)%status /= sys_init) return
         val = sys(i)%c%nmol
      elseif (ic == ic_tree_a) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%aa(1)*bohrtoa
      elseif (ic == ic_tree_b) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%aa(2)*bohrtoa
      elseif (ic == ic_tree_c) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%aa(3)*bohrtoa
      elseif (ic == ic_tree_alpha) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%bb(1)
      elseif (ic == ic_tree_beta) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%bb(2)
      elseif (ic == ic_tree_gamma) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule) return
         val = sys(i)%c%bb(3)
      elseif (ic == ic_tree_emol) then
         if (sysc(i)%status /= sys_init .or. sysc(i)%seed%energy == huge(1d0)) return
         val = sysc(i)%seed%energy/sys(i)%c%nmol
      elseif (ic == ic_tree_p) then
         if (sysc(i)%status /= sys_init .or. sys(i)%c%ismolecule .or.&
             sysc(i)%seed%pressure == huge(1d0)) return
         val = sysc(i)%seed%pressure
      end if
      getvalue = .true.

    end function getvalue

    ! Get value ic for system i and return it in val. The function returns ok if it
    ! is a valid value.
    subroutine getname(str,ic)
      character(len=:,kind=c_char), allocatable, target :: str
      integer, intent(in) :: ic

      if (ic == ic_tree_id) then
         str = "ID" // c_null_char
      elseif (ic == ic_tree_e) then
         str = "Energy (Ha)" // c_null_char
      elseif (ic == ic_tree_v) then
         str = "Volume (Å³)" // c_null_char
      elseif (ic == ic_tree_vmol) then
         str = "Volume/Z (Å³)" // c_null_char
      elseif (ic == ic_tree_nneq) then
         str = "Number of symmetry-unique atoms" // c_null_char
      elseif (ic == ic_tree_ncel) then
         str = "Number of cell atoms" // c_null_char
      elseif (ic == ic_tree_nmol) then
         str = "Number of molecules" // c_null_char
      elseif (ic == ic_tree_a) then
         str = "a (Å)" // c_null_char
      elseif (ic == ic_tree_b) then
         str = "b (Å)" // c_null_char
      elseif (ic == ic_tree_c) then
         str = "c (Å)" // c_null_char
      elseif (ic == ic_tree_alpha) then
         str = "α (°)" // c_null_char
      elseif (ic == ic_tree_beta) then
         str = "β (°)" // c_null_char
      elseif (ic == ic_tree_gamma) then
         str = "γ (°)" // c_null_char
      elseif (ic == ic_tree_emol) then
         str = "Energy/Z (Ha)" // c_null_char
      elseif (ic == ic_tree_p) then
         str = "Pressure (GPa)" // c_null_char
      end if

    end subroutine getname

  end subroutine draw_treeplot

end submodule treeplot
