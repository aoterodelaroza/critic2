! Copyright (c) 2007-2026 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Registry of literature citations: implementation.
submodule (citations) proc
  implicit none

contains

  !> Write a "Please cite:" block for the citations ids(:) to the
  !> standard output. The header line (default "Please cite:") is
  !> written with indent leading blanks (default 2) and each citation
  !> on its own line, two blanks further in. With comment, the block
  !> is written as comment lines instead ("# Please cite:" and
  !> "#   citation"), for output that is a commented table. Unknown
  !> ids are skipped; if none remains, nothing is written.
  module subroutine cite_print(ids,header,indent,comment)
    use tools_io, only: uout
    integer, intent(in) :: ids(:)
    character*(*), intent(in), optional :: header
    integer, intent(in), optional :: indent
    logical, intent(in), optional :: comment

    integer :: i, ind
    integer, allocatable :: ids_(:)
    logical :: comment_
    character(len=:), allocatable :: head, pre

    ind = 2
    if (present(indent)) ind = max(indent,0)
    comment_ = .false.
    if (present(comment)) comment_ = comment
    head = "Please cite:"
    if (present(header)) head = header

    ids_ = pack(ids,ids >= 1 .and. ids <= ncite)
    if (size(ids_,1) == 0) return
    if (comment_) then
       pre = "#"
       write (uout,'(A,A)') "# ", head
    else
       pre = repeat(" ",ind)
       write (uout,'(A,A)') pre, head
    end if
    do i = 1, size(ids_,1)
       write (uout,'(A,A)') pre // "  ", cite_line(ids_(i))
    end do

  end subroutine cite_print

  !> Return the index of the citation with the given key (case
  !> sensitive, trailing blanks ignored), or 0 if there is none.
  module function cite_find(key) result(id)
    character*(*), intent(in) :: key
    integer :: id

    id = findloc(citelist%key,key,dim=1)

  end function cite_find

  !> The one-line form of citation id: the reference text followed by
  !> the DOI in parentheses, if there is one.
  module function cite_line(id) result(str)
    integer, intent(in) :: id
    character(len=:), allocatable :: str

    str = ""
    if (id < 1 .or. id > ncite) return
    str = trim(citelist(id)%text)
    if (len_trim(citelist(id)%doi) > 0) &
       str = str // " (doi:" // trim(citelist(id)%doi) // ")"

  end function cite_line

  !> The https://doi.org/ URL of citation id, or an empty string if
  !> the citation has no DOI.
  module function cite_doi_url(id) result(str)
    integer, intent(in) :: id
    character(len=:), allocatable :: str

    str = ""
    if (id < 1 .or. id > ncite) return
    if (len_trim(citelist(id)%doi) == 0) return
    str = "https://doi.org/" // url_encode(trim(citelist(id)%doi))

  contains
    ! percent-encode the characters a DOI may contain that are not
    ! allowed as they are in a URL
    function url_encode(doi) result(str)
      character(len=*), intent(in) :: doi
      character(len=:), allocatable :: str

      integer :: i

      str = ""
      do i = 1, len(doi)
         if (doi(i:i) == "<") then
            str = str // "%3C"
         elseif (doi(i:i) == ">") then
            str = str // "%3E"
         elseif (doi(i:i) == '"') then
            str = str // "%22"
         elseif (doi(i:i) == " ") then
            str = str // "%20"
         else
            str = str // doi(i:i)
         end if
      end do

    end function url_encode

  end function cite_doi_url

end submodule proc
