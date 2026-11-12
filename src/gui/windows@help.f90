! Copyright (c) 2019-2026 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Routines for the help window: one topic per tool window of the GUI,
! read from <CRITIC_HOME>/helpdoc/<topic>.help, with the link to its
! section of the online manual and the references to cite. The Help
! button the window framework draws in every tool window with a topic
! opens this window at that topic (help_open).
!
! Help file format: a header of "name: value" lines (title, keyword,
! manual, cite) up to the first blank line, then the body: blank lines
! separate paragraphs, "# " starts a section heading, "- " a bullet
! (continuation lines indented), anything else is paragraph text.
submodule (windows) help
  use interfaces_cimgui
  implicit none

  ! the topics, in the order of the list; the key is the name of the
  ! file in helpdoc/ and the key window_help_topic gives
  integer, parameter :: nhelp = 23
  character(len=10), parameter :: help_key(nhelp) = (/&
     "index     ",&
     "menus     ","tree      ","console   ","view      ",&
     "geometry  ","builder   ","extract   ",&
     "vibrations","dynamics  ","rattle    ",&
     "voids     ","kpoints   ","kpath     ",&
     "point     ","line      ","planemap  ","isosurface","cube      ","nci       ","mo        ",&
     "cp        ","gpaths    "/)

  ! the group of each topic: a header in the list, as in the Tools menu
  integer, parameter :: nhelpgroup = 7
  character(len=18), parameter :: help_group(nhelpgroup) = (/&
     "Overview          ","Critic2 interface ","Structure         ","Dynamics          ",&
     "Structure analysis","Field plots       ","QTAIM             "/)
  integer, parameter :: help_igroup(nhelp) = (/1, 2,2,2,2, 3,3,3, 4,4,4, 5,5,5, 6,6,6,6,6,6,6, 7,7/)

  ! width of the topic list (characters)
  integer, parameter :: help_list_width = 27

contains

  !> Draw the help window: the list of topics on the left and the
  !> selected topic on the right (title, console keyword, body,
  !> references), with the manual link in the bottom row.
  module subroutine draw_help(w)
    use gui_main, only: fontsize
    use utils, only: iw_text, iw_button, iw_tooltip, iw_setpos_bottomright, iw_close_event,&
       iw_separatortext, iw_citation
    use tools_io, only: string
    use param, only: critic2_url
    class(window), intent(inout), target :: w

    integer :: i, ig, it, ib
    type(ImVec2) :: sz, szero
    character(kind=c_char,len=:), allocatable, target :: str

    logical, save :: ttshown = .false. ! tooltip flag

    ! load the topics on the first pass; they are kept when the window
    ! is closed and reopened (the selected topic too: help_open sets it)
    if (w%firstpass) then
       call help_load_all(w%hp)
       w%hp%itopic = max(1,min(w%hp%itopic,nhelp))
       w%hp%sec_refs = iw_section()
    end if
    szero%x = 0._c_float
    szero%y = 0._c_float

    ! left: the list of topics, grouped
    sz%x = help_list_width * fontsize%x
    sz%y = 0._c_float
    str = "##helptopics" // c_null_char
    if (igBeginChild_Str(c_loc(str),sz,.true._c_bool,ImGuiWindowFlags_None)) then
       do ig = 1, nhelpgroup
          if (ig > 1) call igSpacing()
          call iw_text(trim(help_group(ig)),highlight=.true.)
          call igIndent(0._c_float)
          do it = 1, nhelp
             if (help_igroup(it) /= ig) cycle
             str = w%hp%tp(it)%title // "##helptopic" // string(it) // c_null_char
             if (igSelectable_Bool(c_loc(str),logical(it == w%hp%itopic,c_bool),&
                ImGuiSelectableFlags_None,szero)) then
                w%hp%itopic = it
                w%hp%scrolltop = .true.
             end if
          end do
          call igUnindent(0._c_float)
       end do
    end if
    call igEndChild()

    ! right: the selected topic
    call igSameLine(0._c_float,-1._c_float)
    sz%x = 0._c_float
    sz%y = 0._c_float
    str = "##helpbody" // c_null_char
    if (igBeginChild_Str(c_loc(str),sz,.false._c_bool,ImGuiWindowFlags_None)) then
       if (w%hp%scrolltop) then
          call igSetScrollY_Float(0._c_float)
          w%hp%scrolltop = .false.
       end if
       associate (tp => w%hp%tp(w%hp%itopic))
         ! title and console keyword
         call iw_text(tp%title,highlight=.true.)
         if (len_trim(tp%keyword) > 0) then
            call iw_text("Console keyword: " // tp%keyword)
            call iw_tooltip("The keyword that does the same from the input console; the Help button "//&
               "of the console opens its manual section",ttshown)
         end if

         ! the file could not be read
         if (.not.tp%found) then
            call iw_text("The help file for this topic was not found: " // tp%file,danger=.true.,wrap=.true.)
            call iw_text("Make sure the CRITIC_HOME environment variable points to the critic2 "//&
               "directory or that critic2 is installed properly.",wrap=.true.)
         end if

         ! the body
         do ib = 1, size(tp%block,1)
            if (tp%block(ib)%kind == hb_heading) then
               call iw_separatortext(tp%block(ib)%text)
            elseif (tp%block(ib)%kind == hb_bullet) then
               ! the text item starts after the dash, so its wrapped
               ! lines hang from there
               call igIndent(0._c_float)
               call iw_text("-")
               call iw_text(tp%block(ib)%text,sameline=.true.,wrap=.true.)
               call igUnindent(0._c_float)
            else
               if (ib > 1) then
                  if (tp%block(ib-1)%kind /= hb_heading) call igSpacing()
               end if
               call iw_text(tp%block(ib)%text,wrap=.true.)
            end if
         end do

         ! the references
         if (size(tp%icite,1) > 0) then
            call iw_separatortext("References",w%hp%sec_refs)
            if (w%hp%sec_refs%open) then
               do i = 1, size(tp%icite,1)
                  if (i > 1) call igSpacing()
                  call iw_citation(tp%icite(i),"helpcite" // string(i),ttshown,indent=.true.)
               end do
            end if
         end if
       end associate
    end if
    call igEndChild()

    ! bottom row: the manual link and close
    call iw_setpos_bottomright(19,2)
    if (iw_button("Open in Manual")) then
       str = critic2_url // "manual/" // trim(w%hp%tp(w%hp%itopic)%manual) // c_null_char
       call openLink(c_loc(str))
    end if
    call iw_tooltip("Open the section of the online manual for this topic in the browser (the "//&
       "front page of the manual, if the topic has no section)",ttshown)
    if (iw_button("Close",sameline=.true.)) w%isopen = .false.

    ! exit if focused and received the close keybinding
    if (iw_close_event(w%focused())) w%isopen = .false.

  end subroutine draw_help

  !> Open the help window (or focus it if it is open) at the topic
  !> with the given key; an unknown key opens the overview.
  module subroutine help_open(key)
    character(len=*), intent(in) :: key

    integer :: it

    if (iwin_help < 1 .or. iwin_help > nwin) return
    it = max(1,findloc(help_key,key,dim=1))
    if (win(iwin_help)%isopen) then
       call igSetWindowFocus_Str(c_loc(win(iwin_help)%name))
    else
       call win(iwin_help)%init(wintype_help,.true.,iwin_help)
    end if
    win(iwin_help)%hp%itopic = it
    win(iwin_help)%hp%scrolltop = .true.

  end subroutine help_open

  !> The citation ids of the help topic with the given key (an empty
  !> array if the topic has none or the key is unknown): the Cite
  !> button of a tool window copies them. Loads the topics into the
  !> help window if it has not been opened yet.
  module function help_topic_cites(key) result(ids)
    character(len=*), intent(in) :: key
    integer, allocatable :: ids(:)

    integer :: it

    allocate(ids(0))
    if (iwin_help < 1 .or. iwin_help > nwin) return
    it = findloc(help_key,key,dim=1)
    if (it == 0) return
    call help_load_all(win(iwin_help)%hp)
    ids = win(iwin_help)%hp%tp(it)%icite

  end function help_topic_cites

  !> The key of the help topic of window w, or an empty string if the
  !> window has none (then it gets no Help button).
  module function window_help_topic(w) result(key)
    use representations, only: reptype_slice, reptype_isosurface, reptype_kpath
    class(window), intent(in) :: w
    character(len=:), allocatable :: key

    key = ""
    if (w%type == wintype_tree) then
       key = "tree"
    elseif (w%type == wintype_view) then
       key = "view"
    elseif (w%type == wintype_console_input .or. w%type == wintype_console_output) then
       key = "console"
    elseif (w%type == wintype_geometry) then
       key = "geometry"
    elseif (w%type == wintype_builder) then
       key = "builder"
    elseif (w%type == wintype_extract) then
       key = "extract"
    elseif (w%type == wintype_vibrations) then
       key = "vibrations"
    elseif (w%type == wintype_dynamics) then
       key = "dynamics"
    elseif (w%type == wintype_rattle) then
       key = "rattle"
    elseif (w%type == wintype_voids) then
       key = "voids"
    elseif (w%type == wintype_kpoints) then
       key = "kpoints"
    elseif (w%type == wintype_point) then
       key = "point"
    elseif (w%type == wintype_line) then
       key = "line"
    elseif (w%type == wintype_cube) then
       key = "cube"
    elseif (w%type == wintype_nci) then
       key = "nci"
    elseif (w%type == wintype_mo) then
       key = "mo"
    elseif (w%type == wintype_cp) then
       key = "cp"
    elseif (w%type == wintype_gpaths) then
       key = "gpaths"
    elseif (w%type == wintype_editrep) then
       if (associated(w%rep)) then
          if (w%rep%type == reptype_slice) then
             key = "planemap"
          elseif (w%rep%type == reptype_isosurface) then
             key = "isosurface"
          elseif (w%rep%type == reptype_kpath) then
             key = "kpath"
          end if
       end if
    end if

  end function window_help_topic

  !xx! private procedures

  !> Load every topic into the help state hp, once (the topics are
  !> kept while the window is closed).
  subroutine help_load_all(hp)
    type(help_state), intent(inout) :: hp

    integer :: i

    if (allocated(hp%tp)) return
    allocate(hp%tp(nhelp))
    do i = 1, nhelp
       call help_load_topic(trim(help_key(i)),hp%tp(i))
    end do

  end subroutine help_load_all

  !> Read the help file of topic key into tp. A missing or unreadable
  !> file leaves tp%found false, with the key as title.
  subroutine help_load_topic(key,tp)
    use global, only: critic_home
    use citations, only: cite_find
    use tools_io, only: fopen_read, fclose, getline_raw, lower, getword, ferror, warning
    use param, only: dirsep
    character(len=*), intent(in) :: key
    type(help_topic), intent(out) :: tp

    character(len=:), allocatable :: line, name, val, word, cur
    integer :: lu, idx, curkind, id, lp
    logical :: exist, inheader

    tp%found = .false.
    tp%file = trim(critic_home) // dirsep // "helpdoc" // dirsep // key // ".help"
    tp%title = key
    tp%keyword = ""
    tp%manual = ""
    allocate(tp%icite(0))
    allocate(tp%block(0))

    inquire(file=tp%file,exist=exist)
    if (.not.exist) return
    lu = fopen_read(tp%file,errstop=.false.)
    if (lu < 0) return
    tp%found = .true.

    inheader = .true.
    curkind = 0
    cur = ""
    do while (getline_raw(lu,line,.false.))
       ! a CRLF file (Windows checkout) keeps its CR: drop it
       if (len(line) > 0) then
          if (line(len(line):len(line)) == char(13)) line = line(1:len(line)-1)
       end if
       if (inheader) then
          ! the header: name: value, up to the first blank line
          if (len_trim(line) == 0) then
             inheader = .false.
             cycle
          end if
          idx = index(line,":")
          if (idx == 0) cycle
          name = lower(trim(adjustl(line(1:idx-1))))
          val = trim(adjustl(line(idx+1:)))
          if (name == "title") then
             tp%title = val
          elseif (name == "keyword") then
             tp%keyword = val
          elseif (name == "manual") then
             tp%manual = val
          elseif (name == "cite") then
             lp = 1
             do while (.true.)
                word = getword(val,lp)
                if (len_trim(word) == 0) exit
                id = cite_find(word)
                if (id > 0) then
                   tp%icite = (/tp%icite, id/)
                else
                   call ferror("help_load_topic","unknown citation key " // word // " in " // tp%file,warning)
                end if
             end do
          end if
       else
          ! the body
          if (len_trim(line) == 0) then
             call flush_block()
          elseif (index(line,"# ") == 1) then
             call flush_block()
             tp%block = (/tp%block, help_block(hb_heading,trim(adjustl(line(3:))))/)
          elseif (index(line,"- ") == 1) then
             call flush_block()
             curkind = hb_bullet
             cur = trim(adjustl(line(3:)))
          elseif (curkind == 0) then
             curkind = hb_para
             cur = trim(adjustl(line))
          else
             cur = cur // " " // trim(adjustl(line))
          end if
       end if
    end do
    call flush_block()
    call fclose(lu)

  contains
    ! the pending paragraph or bullet becomes a block
    subroutine flush_block()

      if (curkind /= 0 .and. len_trim(cur) > 0) tp%block = (/tp%block, help_block(curkind,cur)/)
      curkind = 0
      cur = ""

    end subroutine flush_block
  end subroutine help_load_topic

end submodule help
