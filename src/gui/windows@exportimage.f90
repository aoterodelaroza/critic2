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

! Routines for the about window.
submodule (windows) exportimage
  use interfaces_cimgui
  implicit none
contains

  !> Draw the export image window
  module subroutine draw_exportimage(w)
    use interfaces_opengl3
    use interfaces_stb
    use gui_main, only: g
    use systems, only: sys_init, ok_system
    use windows, only: wintype_dialog, wpurp_dialog_saveimagefile
    use utils, only: iw_text, iw_button, iw_calcwidth, iw_tooltip, get_current_working_dir,&
       iw_checkbox
    use keybindings, only: is_bind_event, BIND_CLOSE_FOCUSED_DIALOG, BIND_OK_FOCUSED_DIALOG,&
       BIND_CLOSE_ALL_DIALOGS
    use tools_io, only: ferror, string
    use param, only: dirsep
    class(window), intent(inout), target :: w

    logical :: doquit, ok, goodsys, okvalid
    integer :: isys, width, height, iaux
    integer(c_int), target :: msFBO, endFBO ! framebuffer
    integer(c_int), target :: msFBOdepth, endFBOdepth ! framebuffer, depth buffer
    integer(c_int), target :: msFBOtex, endFBOtex ! framebuffer, texture
    integer(c_signed_char), allocatable, target :: data(:)
    integer(c_int) :: idum, origin(2)
    character(kind=c_char,len=:), allocatable, target :: str, str1, str2
    type(ImVec2) :: x0, x1, szavail
    logical(c_bool) :: ldum

    logical, save :: ttshown = .false. ! tooltip flag

    ! initialize state
    if (w%firstpass) then
       w%okfile = get_current_working_dir() // dirsep // "image.png"
       w%okfilter = "PNG"
       w%nsample = 16
       w%jpgquality = 90
       w%exportview = .true.
       w%npixel = win(w%idparent)%FBOside
       w%transparentbg = .true.
       w%errmsg = ""
    end if

    ! initialize
    doquit = .false.
    if (associated(win(w%idparent)%sc)) then
       isys = win(w%idparent)%sc%id
    else
       isys = win(w%idparent)%view_selected
       doquit = .true.
    end if

    ! Image file button
    call iw_text("Image File",highlight=.true.)
    if (iw_button("File",danger=.true.)) &
       iaux = stack_create_window(wintype_dialog,.true.,wpurp_dialog_saveimagefile,idparent=w%id,orraise=-1)
    call iw_tooltip("Choose the file to save the image to",ttshown)
    call iw_text(w%okfile,sameline=.true.)

    ! render settings
    call iw_text("Render Settings",highlight=.true.)
    str1 = "Render buffer size (pixels)" // c_null_char
    str2 = "%d" // c_null_char
    call igPushItemWidth(iw_calcwidth(5,1))
    ldum = igDragInt(c_loc(str1),w%npixel,10._c_float,1_c_int,81920_c_int,c_loc(str2),ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Size of the render buffer in pixels",ttshown)

    str1 = "Number of samples for anti-aliasing" // c_null_char
    str2 = "%d" // c_null_char
    call igPushItemWidth(iw_calcwidth(2,1))
    ldum = igDragInt(c_loc(str1),w%nsample,1._c_float,1_c_int,32_c_int,c_loc(str2),&
       ImGuiSliderFlags_AlwaysClamp)
    call igPopItemWidth()
    call iw_tooltip("Number of samples for the anti-aliasing render",ttshown)

    ldum = iw_checkbox("Export the viewport only",w%exportview)
    call iw_tooltip("Export the viewport only or the whole render buffer",ttshown)

    ldum = iw_checkbox("Transparent background",w%transparentbg)
    call iw_tooltip("Make the background transparent in the exported image",ttshown)

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
    if (len_trim(w%errmsg) > 0) call iw_text(w%errmsg,danger=.true.)

    ! right-align and bottom-align for the rest of the contents
    call igGetContentRegionAvail(szavail)
    call igSetCursorPosX(iw_calcwidth(8,2,from_end=.true.) - g%Style%ScrollbarSize)
    if (szavail%y > igGetTextLineHeightWithSpacing() + g%Style%WindowPadding%y) &
       call igSetCursorPosY(igGetCursorPosY() + szavail%y - igGetTextLineHeightWithSpacing() - g%Style%WindowPadding%y)

    ! final buttons: OK
    okvalid = (len_trim(w%okfile) > 0)
    ok = (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) .and. okvalid
    ok = ok .or. iw_button("OK",disabled=.not.okvalid)
    if (ok) then
       ! reset the error message
       w%errmsg = ""

       ! generate textures and buffers
       call glGenTextures(1, c_loc(msFBOtex))
       call glGenRenderbuffers(1, c_loc(msFBOdepth))
       call glGenFramebuffers(1, c_loc(msFBO))
       call glGenTextures(1, c_loc(endFBOtex))
       call glGenRenderbuffers(1, c_loc(endFBOdepth))
       call glGenFramebuffers(1, c_loc(endFBO))

       ! textures
       call glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, msFBOtex)
       call glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, w%nsample, GL_RGBA, w%npixel, w%npixel,&
          int(GL_TRUE,c_signed_char))
       call glBindTexture(GL_TEXTURE_2D, 0)
       call glBindTexture(GL_TEXTURE_2D, endFBOtex)
       call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w%npixel, w%npixel, 0, GL_RGBA, GL_UNSIGNED_BYTE, c_null_ptr)
       call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
       call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
       call glBindTexture(GL_TEXTURE_2D, 0)

       ! render buffer
       call glBindRenderbuffer(GL_RENDERBUFFER, msFBOdepth)
       call glRenderbufferStorageMultisample(GL_RENDERBUFFER, w%nsample, GL_DEPTH_COMPONENT, w%npixel, w%npixel)
       call glBindRenderbuffer(GL_RENDERBUFFER, 0)
       call glBindRenderbuffer(GL_RENDERBUFFER, endFBOdepth)
       call glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w%npixel, w%npixel)
       call glBindRenderbuffer(GL_RENDERBUFFER, 0)

       ! frame buffer
       call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
       call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, msFBOtex, 0)
       call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, msFBOdepth)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
          w%errmsg = "Failed creating the multi-sampled framebuffer (too large?)"
          goto 999
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
       call glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, endFBOtex, 0)
       call glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, endFBOdepth)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) then
          w%errmsg = "Failed creating the render framebuffer (too large?)"
          goto 999
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! render the scene to the multisampled framebuffer
       call glBindFramebuffer(GL_FRAMEBUFFER, msFBO)
       call glViewport(0_c_int,0_c_int,w%npixel,w%npixel)
       if (w%transparentbg) then
          call glClearColor(win(w%idparent)%sc%bgcolor(1),win(w%idparent)%sc%bgcolor(2),&
             win(w%idparent)%sc%bgcolor(3),0._c_float)
       else
          call glClearColor(win(w%idparent)%sc%bgcolor(1),win(w%idparent)%sc%bgcolor(2),&
             win(w%idparent)%sc%bgcolor(3),1._c_float)
       end if
       call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
       goodsys = ok_system(isys,sys_init)
       if (goodsys) call win(w%idparent)%sc%render()
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)

       ! blit the multisampled buffer to the normal colorbuffer
       call glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO)
       call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, endFBO)
       call glBlitFramebuffer(0, 0, w%npixel, w%npixel, 0, 0, w%npixel, w%npixel, GL_COLOR_BUFFER_BIT, GL_LINEAR)
       call glBindFramebuffer(GL_READ_FRAMEBUFFER, 0)
       call glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0)

       ! Read from the regular framebuffer into the data array
       call glBindFramebuffer(GL_FRAMEBUFFER, endFBO)
       if (.not.w%exportview) then
          ! whole texture
          width = w%npixel
          height = w%npixel
          allocate(data(4 * width * height))
          data = 0_c_signed_char
          call glReadPixels(0, 0, w%npixel, w%npixel, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
       else
          ! viewport only
          x0%x = win(w%idparent)%v_rmin%x
          x0%y = win(w%idparent)%v_rmin%y
          x1%x = win(w%idparent)%v_rmax%x
          x1%y = win(w%idparent)%v_rmax%y
          call win(w%idparent)%mousepos_to_texpos(x0)
          call win(w%idparent)%mousepos_to_texpos(x1)
          width = min(nint((x1%x - x0%x) / real(win(w%idparent)%FBOside,8) * w%npixel),w%npixel)
          height = min(nint((x1%y - x0%y) / real(win(w%idparent)%FBOside,8) * w%npixel),w%npixel)
          origin(1) = max(nint(x0%x / real(win(w%idparent)%FBOside,8) * w%npixel),0)
          origin(2) = max(nint(x0%y / real(win(w%idparent)%FBOside,8) * w%npixel),0)
          allocate(data(4 * width * height))
          data = 0_c_signed_char
          call glReadPixels(origin(1), origin(2), width, height, GL_RGBA, GL_UNSIGNED_BYTE, c_loc(data))
       end if
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)
       if (glCheckFramebufferStatus(GL_FRAMEBUFFER) /= GL_FRAMEBUFFER_COMPLETE) &
          w%errmsg = "Error rendering export image"

       ! write the file
       str = trim(w%okfile) // c_null_char
       if (w%okfilter(1:3) == "PNG") then
          idum = stbi_write_png(c_loc(str), width, height, 4, c_loc(data), 4*width)
       elseif (w%okfilter(1:3) == "BMP") then
          idum = stbi_write_bmp(c_loc(str), width, height, 4, c_loc(data))
       elseif (w%okfilter(1:3) == "TGA") then
          idum = stbi_write_tga(c_loc(str), width, height, 4, c_loc(data))
       elseif (w%okfilter(1:3) == "JPE") then
          idum = stbi_write_jpg(c_loc(str), width, height, 4, c_loc(data), w%jpgquality)
       end if
       if (idum == 0) &
          w%errmsg = "Error exporting image to file: "  // string(w%okfile)

999    continue

       ! delete the buffers
       call glBindFramebuffer(GL_FRAMEBUFFER, 0)
       call glDeleteTextures(1, c_loc(msFBOtex))
       call glDeleteRenderbuffers(1, c_loc(msFBOdepth))
       call glDeleteFramebuffers(1, c_loc(msFBO))
       call glDeleteTextures(1, c_loc(endFBOtex))
       call glDeleteRenderbuffers(1, c_loc(endFBOdepth))
       call glDeleteFramebuffers(1, c_loc(endFBO))

       ! quit if no error message
       if (len_trim(w%errmsg) == 0) doquit = .true.
    end if

    ! final buttons: cancel
    if (iw_button("Cancel",sameline=.true.)) doquit = .true.

    ! exit if focused and received the close keybinding
    if (w%focused() .and. is_bind_event(BIND_OK_FOCUSED_DIALOG)) doquit = .true.
    if ((w%focused() .and. is_bind_event(BIND_CLOSE_FOCUSED_DIALOG)).or.&
       is_bind_event(BIND_CLOSE_ALL_DIALOGS)) doquit = .true.

    ! quit = close the window
    if (doquit) call w%end()

  end subroutine draw_exportimage

end submodule exportimage
