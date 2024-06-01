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

! Implementation of the critic2 C/C++ library.
submodule (libcritic2) proc
  use crystalmod, only: xrpd_lambda_def, xrpd_fpol_def, xrpd_alpha_def, xrpd_th2ini_def,&
     xrpd_th2end_def, xrpd_besteps_def_safe, xrpd_besteps_def_quick, xrpd_maxfeval_def_safe,&
     xrpd_maxfeval_def_quick, xrpd_max_elong_def_safe, xrpd_max_elong_def_quick,&
     xrpd_max_ang_def_safe, xrpd_max_ang_def_quick
  use iso_c_binding
  implicit none

  ! whether critic2 has been initialized
  logical :: critic2_init = .false.

  !xx! private procedures
  ! subroutine initialize_critic2()

contains

  !> Create a crystal structure object from a file. The format of the
  !> file is detected from the extension. Allocate space for the
  !> crystal structure and return a pointer to it or NULL if there was
  !> an error. Requires destruction of the object after its use.
  module function c2_crystal_from_file(file) bind(c,name="c2_crystal_from_file")
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use c_interface_module, only: c_f_string_alloc
    type(c_ptr), value, intent(in) :: file
    type(c_ptr) :: c2_crystal_from_file

    character(len=:), allocatable :: fname
    type(crystalseed) :: seed
    character(len=:), allocatable :: errmsg
    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()

    ! read seed from file
    c2_crystal_from_file = c_null_ptr
    call c_f_string_alloc(file,fname)
    call seed%read_any_file(fname,-1,errmsg)
    if (len_trim(errmsg) > 0) return

    ! output the crystal
    allocate(crf)
    call crf%struct_new(seed,.false.)
    if (.not.crf%isinit) return
    c2_crystal_from_file = c_loc(crf)

  end function c2_crystal_from_file

  !> Create a crystal structure from the lattice parameters (lattice,
  !> in bohr), number of atoms (natom), atomic positions (position,
  !> fractional coords), and atomic numbers (zat). Allocate space for
  !> the crystal structure and return a pointer to it or NULL if there
  !> was an error. Requires destruction of the object after its use.
  module function c2_crystal_from_lattice(natom,lattice,position,zat) bind(c,name="c2_crystal_from_lattice")
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use tools_io, only: nameguess
    use param, only: maxzat0
    integer(c_int), intent(in), value :: natom
    real(c_double), intent(in) :: lattice(3,3)
    real(c_double), intent(in) :: position(3,natom)
    integer(c_int), intent(in) :: zat(natom)
    type(c_ptr) :: c2_crystal_from_lattice

    integer :: isused(maxzat0)
    type(crystalseed) :: seed
    integer :: i
    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()

    ! build seed from lattice info
    c2_crystal_from_lattice = c_null_ptr
    call seed%end()
    seed%m_x2c = lattice
    seed%useabr = 2

    ! atomic positions and types
    isused = 0
    seed%nat = natom
    seed%nspc = 0
    allocate(seed%x(3,seed%nat),seed%is(seed%nat))
    seed%x = position
    do i = 1, natom
       if (zat(i) <= 0 .or. zat(i) > maxzat0) return
       if (isused(zat(i)) == 0) then
          seed%nspc = seed%nspc + 1
          isused(zat(i)) = seed%nspc
       end if
       seed%is(i) = isused(zat(i))
    end do

    ! atomic species
    allocate(seed%spc(seed%nspc))
    do i = 1, maxzat0
       if(isused(i) > 0) then
          seed%spc(isused(i))%z = i
          seed%spc(isused(i))%name = nameguess(i)
       end if
    end do

    ! rest of info
    seed%isused = .true.
    seed%ismolecule = .false.

    ! output the crystal
    allocate(crf)
    call crf%struct_new(seed,.false.)
    if (.not.crf%isinit) return
    c2_crystal_from_lattice = c_loc(crf)

  end function c2_crystal_from_lattice

  !> Create a crystal structure from the cell lengths (cel, in bohr),
  !> cell angles (ang, degrees), number of atoms (natom), atomic
  !> positions (position, fractional coords), and atomic numbers
  !> (zat). Allocate space for the crystal structure and return a
  !> pointer to it or NULL if there was an error. Requires destruction
  !> of the object after its use.
  module function c2_crystal_from_cellpar(natom,cel,ang,position,zat) bind(c,name="c2_crystal_from_cellpar")
    use tools_math, only: m_x2c_from_cellpar
    integer(c_int), intent(in), value :: natom
    real(c_double), intent(in) :: cel(3)
    real(c_double), intent(in) :: ang(3)
    real(c_double), intent(in) :: position(3,natom)
    integer(c_int), intent(in) :: zat(natom)
    type(c_ptr) :: c2_crystal_from_cellpar

    real(c_double) :: lattice(3,3)

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()

    lattice = m_x2c_from_cellpar(cel,ang)
    c2_crystal_from_cellpar = c2_crystal_from_lattice(natom,lattice,position,zat)

  end function c2_crystal_from_cellpar

  !> Write the report about the input crystal structure to standard output.
  module subroutine c2_describe_crystal(cr) bind(c,name="c2_describe_crystal")
    use crystalmod, only: crystal
    type(c_ptr), value, intent(in) :: cr

    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! write the report
    call crf%report(.true.,.false.)

  end subroutine c2_describe_crystal

  ! Return the number of atoms from crystal cr.
  module function c2_crystal_get_natom(cr) bind(c,name="c2_crystal_get_natom")
    use crystalmod, only: crystal
    type(c_ptr), value, intent(in) :: cr
    integer(c_int) :: c2_crystal_get_natom

    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    c2_crystal_get_natom = crf%ncel

  end function c2_crystal_get_natom

  ! Get the number of atoms (nat), cell lengths (cel, bohr), cell
  ! angles (ang, degrees), lattice vectors (lattice, bohr),
  ! atomic positions (position, fractional coords), and atomic numbers
  ! (zat) from the crystal structure cr. The input variables must point
  ! to enough space to contain the output (cel[3], ang[3], lattice[3][3],
  ! position[natom][3], zat[natom]). If the pointer is NULL, do not
  ! return that variable.
  module subroutine c2_crystal_get_structure(cr,natom,cel,ang,lattice,position,zat) &
     bind(c,name="c2_crystal_get_structure")
    use crystalmod, only: crystal
    type(c_ptr), value, intent(in) :: cr
    integer(c_int) :: natom
    type(c_ptr), value :: cel
    type(c_ptr), value :: ang
    type(c_ptr), value :: lattice
    type(c_ptr), value :: position
    type(c_ptr), value :: zat

    type(crystal), pointer :: crf
    real(c_double), pointer :: fptr(:), fptr2(:,:)
    integer(c_int), pointer :: iptr(:)
    integer :: i

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! number of atoms
    natom = crf%ncel

    ! cell lengths, angles, lattice vectors
    if (c_associated(cel)) then
       call c_f_pointer(cel,fptr,(/3/))
       fptr(1:3) = crf%aa
    end if
    if (c_associated(ang)) then
       call c_f_pointer(ang,fptr,(/3/))
       fptr(1:3) = crf%bb
    end if
    if (c_associated(lattice)) then
       call c_f_pointer(lattice,fptr2,(/3,3/))
       fptr2(1:3,1:3) = crf%m_x2c
    end if

    ! atomic positions and atomic numbers
    if (c_associated(position)) then
       call c_f_pointer(zat,iptr,(/natom/))
       call c_f_pointer(position,fptr2,(/3,natom/))
       do i = 1, crf%ncel
          iptr(i) = crf%spc(crf%atcel(i)%is)%z
          fptr2(1:3,i) = crf%atcel(i)%x
       end do
    end if

  end subroutine c2_crystal_get_structure

  !> Write the crystal structure to a file. The format of the file is
  !> detected from the extension.
  module subroutine c2_write_crystal(cr,file) bind(c,name="c2_write_crystal")
    use crystalmod, only: crystal
    use c_interface_module, only: c_f_string_alloc
    type(c_ptr), value, intent(in) :: cr
    type(c_ptr), value, intent(in) :: file

    type(crystal), pointer :: crf
    character(len=:), allocatable :: fname

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! write the file
    call c_f_string_alloc(file,fname)
    call crf%write_simple_driver(fname)

  end subroutine c2_write_crystal

  !> Destroy the input crystal structure object and free the memory.
  module subroutine c2_destroy_crystal(cr) bind(c,name="c2_destroy_crystal")
    use crystalmod, only: crystal
    type(c_ptr), value, intent(in) :: cr

    type(crystal), pointer :: crf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! destroy the crystal object and free memory
    call crf%end()
    deallocate(crf)

  end subroutine c2_destroy_crystal

  !> Calculate the XRPD peak positions from the crystal structure.
  !> Limit the peaks to the 2theta range th2ini to th2end. lambda is
  !> the wavelength in angstrom. fpol is the polarization correction
  !> factor (0 = unpolarized, 0.95 = synchrotron). If th2ini, th2end,
  !> lambda, fpol < 0, use default values. Returns NULL on error.
  module function c2_peaks_from_crystal(cr,th2ini,th2end,lambda,fpol) bind(c,name="c2_peaks_from_crystal")
    use crystalmod, only: crystal, xrpd_peaklist
    type(c_ptr), value, intent(in) :: cr
    real(c_double), value :: th2ini, th2end, lambda, fpol
    type(c_ptr) :: c2_peaks_from_crystal

    type(crystal), pointer :: crf
    type(xrpd_peaklist), pointer :: pk
    character(len=:), allocatable :: errmsg

    ! consistency checks
    c2_peaks_from_crystal = c_null_ptr
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(cr)) return
    call c_f_pointer(cr,crf)
    if (.not.associated(crf)) return

    ! handle default values
    if (th2ini < 0._c_double) th2ini = xrpd_th2ini_def
    if (th2end < 0._c_double) th2end = xrpd_th2end_def
    if (lambda < 0._c_double) lambda = xrpd_lambda_def
    if (fpol < 0._c_double) fpol = xrpd_fpol_def

    ! calculate peak positions and return
    allocate(pk)
    call pk%from_crystal(crf,th2ini,th2end,lambda,fpol,.false.,.false.,errmsg)
    if (len_trim(errmsg) > 0) return
    c2_peaks_from_crystal = c_loc(pk)

  end function c2_peaks_from_crystal

  !> Read xrpd peaks from a file. Allocate space for the crystal
  !> structure and return a pointer to it or NULL if there was an
  !> error.
  module function c2_peaks_from_file(file) bind(c,name="c2_peaks_from_file")
    use crystalmod, only: xrpd_peaklist
    use c_interface_module, only: c_f_string_alloc
    type(c_ptr), value, intent(in) :: file
    type(c_ptr) :: c2_peaks_from_file

    character(len=:), allocatable :: fname, errmsg
    type(xrpd_peaklist), pointer :: pk

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()

    ! read seed from file
    c2_peaks_from_file = c_null_ptr
    allocate(pk)
    call c_f_string_alloc(file,fname)
    ! xxxx !
    ! call xrpd_peaks_from_file(pk,fname,errmsg)
    if (len_trim(errmsg) > 0) return
    c2_peaks_from_file = c_loc(pk)

  end function c2_peaks_from_file

  !> Destroy the input XRPD peaks structure object and free the memory.
  module subroutine c2_destroy_peaks(pk) bind(c,name="c2_destroy_peaks")
    use crystalmod, only: xrpd_peaklist
    type(c_ptr), value, intent(in) :: pk

    type(xrpd_peaklist), pointer :: pkf

    ! consistency checks
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(pk)) return
    call c_f_pointer(pk,pkf)
    if (.not.associated(pkf)) return

    ! destroy the crystal object and free memory
    deallocate(pkf)

  end subroutine c2_destroy_peaks

  !> Compare crystal c1 and set of XRPD peaks p2 using GPWDF. alpha =
  !> Gaussian triangle width. lambda = wavelength in angstrom. fpol =
  !> polarization correction factor (0 = unpolarized, 0.95 =
  !> synchrotron). If alpha, lambda, fpol < 0, use default values.
  !> Returns -1 on error.
  module function c2_compare_gpwdf(c1,p2,alpha,lambda,fpol) bind(c,name="c2_compare_gpwdf")
    use crystalmod, only: crystal, xrpd_peaklist, gaussian_compare
    type(c_ptr), value, intent(in) :: c1
    type(c_ptr), value, intent(in) :: p2
    real(c_double), value :: alpha, lambda, fpol
    real(c_double) :: c2_compare_gpwdf

    type(crystal), pointer :: cr
    type(xrpd_peaklist), pointer :: pk
    character(len=:), allocatable :: errmsg

    ! consistency checks
    c2_compare_gpwdf = -1._c_double
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(c1)) return
    call c_f_pointer(c1,cr)
    if (.not.associated(cr)) return
    if (.not.c_associated(p2)) return
    call c_f_pointer(p2,pk)
    if (.not.associated(pk)) return

    ! set default values
    if (alpha < 0d0) alpha = xrpd_alpha_def
    if (lambda < 0d0) lambda = xrpd_lambda_def
    if (fpol < 0d0) fpol = xrpd_fpol_def

    ! run the comparison
    call gaussian_compare(cr,pk,0,c2_compare_gpwdf,errmsg,verbose0=.false.,&
       alpha0=alpha,lambda0=lambda,fpol0=fpol)
    if (len_trim(errmsg) > 0) then
       c2_compare_gpwdf = -1._c_double
       return
    end if

  end function c2_compare_gpwdf

  !> Compare crystal c1 and set of XRPD peaks p2 using variable-cell
  !> (VC-)GPWDF. If global, use the global minimization for the best
  !> GPWDF; otherwise use only a local minimization. If verbose, print
  !> search progress to stdout. alpha = Gaussian triangle
  !> width. lambda = wavelength in angstrom. fpol = polarization
  !> correction factor (0 = unpolarized, 0.95 = synchrotron). maxfeval
  !> = maximum number of function evaluations. besteps = do not
  !> restart the feval count if a diff lower than the best diff is
  !> found within besteps. max_elong = maximum cell length elongation
  !> (%). max_ang = maximum cell angle deformation (degrees). Returns
  !> the VC-GPDWF score and the deformed c1 structure in crout.
  !> Returns -1 on error.
  module function c2_compare_vcgpwdf(c1,p2,crout,global,verbose,alpha,lambda,fpol,maxfeval,&
     besteps,max_elong,max_ang) bind(c,name="c2_compare_vcgpwdf")
    use crystalseedmod, only: crystalseed
    use crystalmod, only: crystal, xrpd_peaklist, gaussian_compare
    type(c_ptr), value, intent(in) :: c1
    type(c_ptr), value, intent(in) :: p2
    type(c_ptr) :: crout
    logical(c_bool), value :: global, verbose
    real(c_double), value :: alpha, lambda, fpol, besteps, max_elong, max_ang
    integer(c_int), value :: maxfeval
    real(c_double) :: c2_compare_vcgpwdf

    type(crystal), pointer :: cr, croutf
    type(xrpd_peaklist), pointer :: pk
    integer :: imode
    type(crystalseed) :: seed
    character(len=:), allocatable :: errmsg

    ! consistency checks
    c2_compare_vcgpwdf = -1._c_double
    if (.not.critic2_init) call initialize_critic2()
    if (.not.c_associated(c1)) return
    call c_f_pointer(c1,cr)
    if (.not.associated(cr)) return
    if (.not.c_associated(p2)) return
    call c_f_pointer(p2,pk)
    if (.not.associated(pk)) return

    ! set default values
    if (alpha < 0d0) alpha = xrpd_alpha_def
    if (lambda < 0d0) lambda = xrpd_lambda_def
    if (fpol < 0d0) fpol = xrpd_fpol_def
    if (maxfeval < 0) maxfeval = xrpd_maxfeval_def_safe
    if (besteps < 0d0) besteps = xrpd_besteps_def_safe
    if (max_elong < 0d0) max_elong = xrpd_max_elong_def_safe
    if (max_ang < 0d0) max_ang = xrpd_max_ang_def_safe
    if (global) then
       imode = 2
    else
       imode = 1
    end if

    ! run the comparison
    call gaussian_compare(cr,pk,imode,c2_compare_vcgpwdf,errmsg,seed,logical(verbose),&
       alpha,lambda,fpol,maxfeval,besteps,max_elong,max_ang)
    if (len_trim(errmsg) > 0) then
       c2_compare_vcgpwdf = -1._c_double
       return
    end if

    ! crystal in output
    allocate(croutf)
    call croutf%struct_new(seed,.false.)
    if (croutf%isinit) then
       crout = c_loc(croutf)
    else
       crout = c_null_ptr
    endif

  end function c2_compare_vcgpwdf

  !> Compare crystal c1 and set of XRPD peaks p2 using variable-cell
  !> (VC-)GPWDF, global minimization. Use "safe" default settings. If
  !> verbose, print search progress to stdout. lambda = wavelength in
  !> angstrom. fpol = polarization correction factor (0 = unpolarized,
  !> 0.95 = synchrotron). Returns the VC-GPDWF score and the deformed
  !> c1 structure in crout. Returns -1 on error.
  module function c2_compare_vcgpwdf_global_safe(c1,p2,crout,verbose,lambda,fpol) &
     bind(c,name="c2_compare_vcgpwdf_global_safe")
    type(c_ptr), value, intent(in) :: c1
    type(c_ptr), value, intent(in) :: p2
    type(c_ptr) :: crout
    logical(c_bool), value :: verbose
    real(c_double), value :: lambda, fpol
    real(c_double) :: c2_compare_vcgpwdf_global_safe

    c2_compare_vcgpwdf_global_safe = c2_compare_vcgpwdf(c1,p2,crout,.true._c_bool,verbose,&
       xrpd_alpha_def,lambda,fpol,xrpd_maxfeval_def_safe,xrpd_besteps_def_safe,xrpd_max_elong_def_safe,&
       xrpd_max_ang_def_safe)

  end function c2_compare_vcgpwdf_global_safe

  !> Compare crystal c1 and set of XRPD peaks p2 using variable-cell
  !> (VC-)GPWDF, global minimization. Use "quick" default settings. If
  !> verbose, print search progress to stdout. lambda = wavelength in
  !> angstrom. fpol = polarization correction factor (0 = unpolarized,
  !> 0.95 = synchrotron). Returns the VC-GPDWF score and the deformed
  !> c1 structure in crout. Returns -1 on error.
  module function c2_compare_vcgpwdf_global_quick(c1,p2,crout,verbose,lambda,fpol) &
     bind(c,name="c2_compare_vcgpwdf_global_quick")
    type(c_ptr), value, intent(in) :: c1
    type(c_ptr), value, intent(in) :: p2
    type(c_ptr) :: crout
    logical(c_bool), value :: verbose
    real(c_double), value :: lambda, fpol
    real(c_double) :: c2_compare_vcgpwdf_global_quick

    c2_compare_vcgpwdf_global_quick = c2_compare_vcgpwdf(c1,p2,crout,.true._c_bool,verbose,&
       xrpd_alpha_def,lambda,fpol,xrpd_maxfeval_def_quick,xrpd_besteps_def_quick,xrpd_max_elong_def_quick,&
       xrpd_max_ang_def_quick)

  end function c2_compare_vcgpwdf_global_quick

  !xx! private procedures

  !> Initialize critic2 when used as a library.
  subroutine initialize_critic2()
    use spglib, only: spg_build_hall_mapping
    use systemmod, only: systemmod_init
    use global, only: global_init, fileroot, initial_banner, config_write
    use tools_io, only: start_clock, stdargs, history_init
    use param, only: param_init
    use config, only: getstring, istring_datadir

    character(len=:), allocatable :: optv, ghome

    ! initialize
    call start_clock()
    call param_init()

    call stdargs(optv,ghome,fileroot)
    call history_init()

    call global_init(ghome,getstring(istring_datadir))
    call systemmod_init(1)
    call spg_build_hall_mapping() ! for spglib
    critic2_init = .true.

  end subroutine initialize_critic2

end submodule proc
