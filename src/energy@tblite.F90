! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> tblite (GFN-FF / GFN2-xTB) backend for the energy calculator. The real
!> implementation is compiled only with -DHAVE_TBLITE; otherwise these
!> procedures are stubs that report the missing feature so the calculator can
!> fall back to the built-in force field.
submodule (energy) tblite
  implicit none

contains

  module subroutine calc_init_tblite(cl,c,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    character(len=:), allocatable, intent(out), optional :: errmsg

#ifdef HAVE_TBLITE
    ! real implementation to be added in the tblite integration step
    if (present(errmsg)) errmsg = "tblite backend not yet implemented"
#else
    if (present(errmsg)) errmsg = "critic2 was compiled without tblite support"
#endif

  end subroutine calc_init_tblite

  module subroutine calc_eval_tblite(cl,c,ene,grad,stress,errmsg)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
    real*8, intent(out) :: ene
    real*8, intent(out) :: grad(:,:)
    real*8, intent(out), optional :: stress(3,3)
    character(len=:), allocatable, intent(out), optional :: errmsg

    ene = 0d0
    grad = 0d0
    if (present(stress)) stress = 0d0
#ifdef HAVE_TBLITE
    if (present(errmsg)) errmsg = "tblite backend not yet implemented"
#else
    if (present(errmsg)) errmsg = "critic2 was compiled without tblite support"
#endif

  end subroutine calc_eval_tblite

  module subroutine calc_update_tblite(cl,c)
    use crystalmod, only: crystal
    class(calculator), intent(inout) :: cl
    class(crystal), intent(inout) :: c
  end subroutine calc_update_tblite

  module subroutine calc_free_tblite(cl)
    class(calculator), intent(inout) :: cl
  end subroutine calc_free_tblite

end submodule tblite
