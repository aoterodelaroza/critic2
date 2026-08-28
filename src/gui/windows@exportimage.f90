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

! Routines for the export image window.
submodule (windows) exportimage
  use interfaces_cimgui
  implicit none

  ! settings remembered from the last successful export (this session)
  integer(c_int) :: lastnsample = 16
  integer(c_int) :: lastjpgquality = 90
  integer(c_int) :: lastnpixel = 0 ! 0 = follow the view's buffer size
  logical :: lastexportview = .true.
  logical :: lasttransparentbg = .true.

contains

  !> Draw the export image window
  module subroutine draw_exportimage(w)
    use windows, only: wintype_dialog, wpurp_dialog_saveimagefile
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, iw_checkbox,&
       iw_close_event, iw_setpos_bottomright, iw_inputtext
    use keybindings, only: is_bind_event, BIND_OK_FOCUSED_DIALOG
    use tools_io, only: string, uout
    use param, only: dirsep
    class(window), intent(inout), target :: w

    logical :: doquit, ok, okvalid, syschanged
    integer :: isys, iview, iaux, iwidth, iheight
    character(kind=c_char,len=:), allocatable, target :: str1, str2
    logical(c_bool) :: ldum

    logical, save :: ttshown = .false. ! tooltip flag

    ! this window renders its anchor view; resolve it first, the first-pass
    ! default for the render buffer size is taken from it
    doquit = .not.w%anchor(iview,isys,syschanged)

    ! initialize state
    if (w%firstpass) then
       w%okfile = okfile_default(isys,"image","png")
       w%okfilter = "PNG"
       w%errmsg = ""
       w%nsample = lastnsample
       w%jpgquality = lastjpgquality
       w%exportview = lastexportview
       w%transparentbg = lasttransparentbg
       if (lastnpixel > 0) then
          w%npixel = lastnpixel
       elseif (.not.doquit) then
          w%npixel = win(iview)%FBOside
       end if
    end if

    ! there is nothing to export without a scene
    if (.not.doquit) doquit = .not.associated(win(iview)%sc)

    ! file name: editable field plus browse button
    call iw_text("File Name",highlight=.true.)
    if (iw_inputtext("##exportimagefile",bufsize=1023,texta=w%okfile,width=36)) then
       ! update the image format from the typed extension, if recognized
       str1 = filter_from_extension()
       if (len_trim(str1) > 0) w%okfilter = str1
    end if
    call iw_tooltip("File the image will be written to (the extension selects the image format)",ttshown)
    if (iw_button("Browse...",sameline=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_saveimagefile,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the image file with a file browser",ttshown)
    call w%okfile_warn_overwrite()
    if (len_trim(filter_from_extension()) == 0) &
       call iw_text("Unknown image extension, writing " // trim(w%okfilter) // " format",&
          danger=.true.,wrap=.true.)

    ! render settings
    call iw_text("Render Settings",highlight=.true.)
    str1 = "Render buffer size (pixels)" // c_null_char
    str2 = "%d" // c_null_char
    call igPushItemWidth(iw_calcwidth(5,1))
    ldum = igDragInt(c_loc(str1),w%npixel,10._c_float,1_c_int,81920_c_int,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Size of the render buffer in pixels",ttshown)

    ! render buffer size presets, relative to the view size
    if (.not.doquit) then
       call iw_text("Presets:",alignframe=.true.)
       if (iw_button("1x",sameline=.true.)) w%npixel = win(iview)%FBOside
       call iw_tooltip("Same resolution as the current view",ttshown)
       if (iw_button("2x",sameline=.true.)) w%npixel = 2 * win(iview)%FBOside
       call iw_tooltip("Twice the resolution of the current view",ttshown)
       if (iw_button("4x",sameline=.true.)) w%npixel = 4 * win(iview)%FBOside
       call iw_tooltip("Four times the resolution of the current view",ttshown)
    end if

    str1 = "Number of samples for anti-aliasing" // c_null_char
    str2 = "%d" // c_null_char
    call igPushItemWidth(iw_calcwidth(2,1))
    ldum = igDragInt(c_loc(str1),w%nsample,1._c_float,1_c_int,32_c_int,c_loc(str2),&
       ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Number of samples for the anti-aliasing render",ttshown)

    ldum = iw_checkbox("Export the viewport only",w%exportview)
    call iw_tooltip("If checked, export the only the current view. Otherwise,&
       & export the whole render buffer",ttshown)

    ldum = iw_checkbox("Transparent background",w%transparentbg)
    call iw_tooltip("Make the background transparent in the exported image",ttshown)

    ! size of the written image
    if (.not.doquit) then
       call win(iview)%export_image_size(w%npixel,w%exportview,iwidth,iheight)
       call iw_text("Output image size:",highlight=.true.)
       call iw_text(string(iwidth) // " x " // string(iheight) // " pixels",sameline=.true.)
    end if

    ! image settings
    if (w%okfilter(1:3) == "JPE") then
       call iw_text("Image Settings",highlight=.true.)
       str1 = "JPEG Quality" // c_null_char
       str2 = "%d" // c_null_char
       call igPushItemWidth(iw_calcwidth(3,1))
       ldum = igDragInt(c_loc(str1),w%jpgquality,1._c_float,1_c_int,100_c_int,c_loc(str2),&
          ImGuiSliderFlags_AlwaysClamp)
       call igPopItemWidth()
       call iw_tooltip("Quality and weight of the JPEG file",ttshown)
    end if

    ! maybe the error message
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.,wrap=.true.)

    ! right-align and bottom-align for the rest of the contents
    call iw_setpos_bottomright(8,2)

    ! final buttons: OK
    okvalid = (len_trim(w%okfile) > 0) .and. .not.doquit
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .and. okvalid
    ok = ok .or. iw_button("OK",disabled=.not.okvalid)
    if (ok) then
       ! export the image
       call win(iview)%export_to_image(w%okfile,w%okfilter,w%nsample,w%npixel,&
          w%transparentbg,w%exportview,w%jpgquality,w%errmsg)

       ! on success, remember the settings, report, and quit
       if (len_trim(w%errmsg) == 0) then
          lastnsample = w%nsample
          lastjpgquality = w%jpgquality
          lastexportview = w%exportview
          lasttransparentbg = w%transparentbg
          lastnpixel = merge(0_c_int,w%npixel,w%npixel == win(iview)%FBOside)
          call okfile_save_dir(w%okfile)
          write (uout,'("Exported image file: ",A)') trim(w%okfile)
          doquit = .true.
       end if
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  contains
    ! image format corresponding to okfile's extension ("" = not recognized)
    function filter_from_extension() result(filter)
      use tools_io, only: equal, lower
      character(len=:), allocatable :: filter

      character(len=:), allocatable :: ext
      integer :: i0, i1

      filter = ""
      i0 = index(w%okfile,dirsep,back=.true.)
      i1 = index(w%okfile(i0+1:),'.',back=.true.)
      if (i1 == 0) return
      ext = lower(w%okfile(i0+i1+1:))
      if (equal(ext,'png')) then
         filter = "PNG"
      elseif (equal(ext,'bmp')) then
         filter = "BMP"
      elseif (equal(ext,'tga')) then
         filter = "TGA"
      elseif (equal(ext,'jpg') .or. equal(ext,'jpeg')) then
         filter = "JPEG"
      end if

    end function filter_from_extension
  end subroutine draw_exportimage

end submodule exportimage
