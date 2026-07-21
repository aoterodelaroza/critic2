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

! Icon textures for the GUI
submodule (icons) proc
  use iso_c_binding
  implicit none

  ! property and UI icon PNG file names, in icon ID order
  character(len=*), parameter :: iconfile(icon_NUM) = (/character(len=20) ::&
     "prop_fields.png","prop_vib.png","prop_occ.png",&
     "ui_close.png","ui_expand.png","ui_collapse.png",&
     "ui_atoms.png","ui_bonds.png","ui_labels.png",&
     "ui_cell.png","ui_polyhedra.png","ui_label_num.png",&
     "ui_label_wyck.png","ui_camera.png","ui_bgcolor.png",&
     "ui_applyall.png","ui_reset.png","ui_draw.png",&
     "ui_objects.png","ui_tools.png","ui_newview.png"/)

  ! format icon PNG file names, in isformat_r_* order
  character(len=*), parameter :: fmtfile(0:icon_fmt_MAX) = (/character(len=20) ::&
     "fmt_unknown.png","fmt_input.png","fmt_library.png","fmt_derived.png",&
     "fmt_cif.png","fmt_shelx.png","fmt_f21.png","fmt_cube.png",&
     "fmt_bincube.png","fmt_struct.png","fmt_abinit.png","fmt_elk.png",&
     "fmt_qein.png","fmt_qeout.png","fmt_crystal.png","fmt_fploout.png",&
     "fmt_xyz.png","fmt_wfn.png","fmt_wfx.png","fmt_fchk.png",&
     "fmt_molden.png","fmt_gaussian.png","fmt_siesta.png","fmt_xsf.png",&
     "fmt_gen.png","fmt_vasp.png","fmt_pwc.png","fmt_axsf.png",&
     "fmt_dat.png","fmt_pgout.png","fmt_orca.png","fmt_dmain.png",&
     "fmt_aimsin.png","fmt_aimsout.png","fmt_tinkerfrac.png","fmt_gjf.png",&
     "fmt_castepcell.png","fmt_castepgeom.png","fmt_mol2.png","fmt_pdb.png",&
     "fmt_zmat.png","fmt_sdf.png","fmt_magres.png","fmt_alamode.png",&
     "fmt_castepphonon.png","fmt_akaikkr.png","fmt_xband.png"/)

  ! tint colors for the format icons, by code family
  real(c_float), parameter :: rgba_fam_internal(4) = (/0.70_c_float,0.70_c_float,0.70_c_float,1.0_c_float/) ! grey
  real(c_float), parameter :: rgba_fam_crys(4) = (/1.00_c_float,0.65_c_float,0.25_c_float,1.0_c_float/)     ! orange
  real(c_float), parameter :: rgba_fam_molgeom(4) = (/0.40_c_float,0.70_c_float,1.00_c_float,1.0_c_float/)  ! light blue
  real(c_float), parameter :: rgba_fam_qchem(4) = (/0.75_c_float,0.55_c_float,1.00_c_float,1.0_c_float/)    ! violet
  real(c_float), parameter :: rgba_fam_pwdft(4) = (/1.00_c_float,0.45_c_float,0.45_c_float,1.0_c_float/)    ! red
  real(c_float), parameter :: rgba_fam_aedft(4) = (/0.25_c_float,0.85_c_float,0.80_c_float,1.0_c_float/)    ! teal
  real(c_float), parameter :: rgba_fam_grid(4) = (/0.45_c_float,0.90_c_float,0.35_c_float,1.0_c_float/)     ! green
  real(c_float), parameter :: rgba_fam_vib(4) = (/0.90_c_float,0.85_c_float,0.20_c_float,1.0_c_float/)      ! yellow

contains

  !> Load the icon PNG files and upload them to OpenGL textures. Must
  !> be called with a live OpenGL context. Missing or unreadable
  !> files leave the corresponding icon_tex entry at zero.
  module subroutine icons_init()
    use interfaces_stb, only: stbi_load, stbi_image_free
    use interfaces_opengl3
    use global, only: critic_home
    use tools_io, only: ferror, warning
    use param, only: dirsep
    integer :: i

    icon_tex = 0
    do i = 1, icon_NUM
       call load_icon(iconfile(i),icon_tex(i))
    end do
    icon_tex_fmt = 0
    do i = 0, icon_fmt_MAX
       call load_icon(fmtfile(i),icon_tex_fmt(i))
       rgba_icon_fmt(:,i) = format_tint(i)
    end do

  contains
    !> Load one icon PNG from dat/assets/icons and upload it to an
    !> OpenGL texture. On failure, warn and leave tex at zero.
    subroutine load_icon(filename,tex)
      character(len=*), intent(in) :: filename
      integer(c_int), intent(inout), target :: tex

      integer(c_int) :: iwidth, iheight, ichannel
      type(c_ptr) :: pixels
      character(kind=c_char,len=:), allocatable, target :: file

      file = trim(critic_home) // dirsep // "assets" // dirsep // "icons" //&
         dirsep // trim(filename) // c_null_char
      pixels = stbi_load(c_loc(file),iwidth,iheight,ichannel,4)
      if (.not.c_associated(pixels)) then
         call ferror('icons_init','could not load icon file: ' // trim(filename),warning)
         return
      end if
      call glGenTextures(1, c_loc(tex))
      call glBindTexture(GL_TEXTURE_2D, tex)
      call glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
      call glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, iwidth, iheight, 0, GL_RGBA,&
         GL_UNSIGNED_BYTE, pixels)
      call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
      call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
      call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
      call glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
      call glBindTexture(GL_TEXTURE_2D, 0)
      call stbi_image_free(pixels)

    end subroutine load_icon
  end subroutine icons_init

  !> Delete the icon textures.
  module subroutine icons_end()
    use interfaces_opengl3
    integer :: i

    do i = 1, icon_NUM
       if (icon_tex(i) /= 0) call glDeleteTextures(1, c_loc(icon_tex(i)))
    end do
    icon_tex = 0
    do i = 0, icon_fmt_MAX
       if (icon_tex_fmt(i) /= 0) call glDeleteTextures(1, c_loc(icon_tex_fmt(i)))
    end do
    icon_tex_fmt = 0

  end subroutine icons_end

  !> Return the tint color for the icon of the given structure reading
  !> format (isformat_r_*). Formats from the same code family share
  !> the same hue.
  function format_tint(isformat) result(rgba)
    use param
    integer, intent(in) :: isformat
    real(c_float) :: rgba(4)

    select case (isformat)
    case (isformat_r_cif,isformat_r_shelx,isformat_r_f21,isformat_r_dmain,&
       isformat_r_pdb,isformat_r_tinkerfrac)
       rgba = rgba_fam_crys
    case (isformat_r_xyz,isformat_r_gjf,isformat_r_zmat,isformat_r_mol2,&
       isformat_r_sdf)
       rgba = rgba_fam_molgeom
    case (isformat_r_wfn,isformat_r_wfx,isformat_r_fchk,isformat_r_molden,&
       isformat_r_gaussian,isformat_r_orca,isformat_r_dat,isformat_r_pgout)
       rgba = rgba_fam_qchem
    case (isformat_r_abinit,isformat_r_qein,isformat_r_qeout,isformat_r_crystal,&
       isformat_r_siesta,isformat_r_gen,isformat_r_vasp,isformat_r_pwc,&
       isformat_r_aimsin,isformat_r_aimsout,isformat_r_castepcell,isformat_r_castepgeom)
       rgba = rgba_fam_pwdft
    case (isformat_r_struct,isformat_r_elk,isformat_r_fploout,isformat_r_akaikkr,&
       isformat_r_xband)
       rgba = rgba_fam_aedft
    case (isformat_r_cube,isformat_r_bincube,isformat_r_xsf,isformat_r_axsf)
       rgba = rgba_fam_grid
    case (isformat_r_magres,isformat_r_alamode,isformat_r_castepphonon)
       rgba = rgba_fam_vib
    case default
       rgba = rgba_fam_internal
    end select

  end function format_tint

  !> Return a human-readable name for the given structure reading
  !> format (isformat_r_*).
  module function format_name(isformat) result(name)
    use param
    integer, intent(in) :: isformat
    character(len=:), allocatable :: name

    select case (isformat)
    case (isformat_r_from_input)
       name = "critic2 input"
    case (isformat_r_from_library)
       name = "critic2 structure library"
    case (isformat_r_derived)
       name = "derived from another system"
    case (isformat_r_cif)
       name = "CIF file"
    case (isformat_r_shelx)
       name = "SHELX file"
    case (isformat_r_f21)
       name = "DMACRYS .21 file"
    case (isformat_r_cube)
       name = "cube file"
    case (isformat_r_bincube)
       name = "binary cube file"
    case (isformat_r_struct)
       name = "WIEN2k struct file"
    case (isformat_r_abinit)
       name = "abinit DEN-style file"
    case (isformat_r_elk)
       name = "elk GEOMETRY.OUT file"
    case (isformat_r_qein)
       name = "Quantum ESPRESSO input"
    case (isformat_r_qeout)
       name = "Quantum ESPRESSO output"
    case (isformat_r_crystal)
       name = "CRYSTAL output"
    case (isformat_r_fploout)
       name = "FPLO output"
    case (isformat_r_xyz)
       name = "xyz file"
    case (isformat_r_wfn)
       name = "Gaussian wfn file"
    case (isformat_r_wfx)
       name = "Gaussian wfx file"
    case (isformat_r_fchk)
       name = "Gaussian fchk file"
    case (isformat_r_molden)
       name = "Molden-style file"
    case (isformat_r_gaussian)
       name = "Gaussian output"
    case (isformat_r_siesta)
       name = "SIESTA STRUCT file"
    case (isformat_r_xsf)
       name = "Xcrysden xsf file"
    case (isformat_r_gen)
       name = "DFTB+ gen file"
    case (isformat_r_vasp)
       name = "VASP file"
    case (isformat_r_pwc)
       name = "Quantum ESPRESSO pwc file"
    case (isformat_r_axsf)
       name = "Xcrysden axsf file"
    case (isformat_r_dat)
       name = "psi4 output"
    case (isformat_r_pgout)
       name = "postg output"
    case (isformat_r_orca)
       name = "ORCA output"
    case (isformat_r_dmain)
       name = "DMACRYS dmain file"
    case (isformat_r_aimsin)
       name = "FHIaims input"
    case (isformat_r_aimsout)
       name = "FHIaims output"
    case (isformat_r_tinkerfrac)
       name = "TINKER frac file"
    case (isformat_r_gjf)
       name = "Gaussian input"
    case (isformat_r_castepcell)
       name = "CASTEP cell file"
    case (isformat_r_castepgeom)
       name = "CASTEP geom file"
    case (isformat_r_mol2)
       name = "mol2 file"
    case (isformat_r_pdb)
       name = "pdb file"
    case (isformat_r_zmat)
       name = "Gaussian zmat file"
    case (isformat_r_sdf)
       name = "sdf file"
    case (isformat_r_magres)
       name = "magres file"
    case (isformat_r_alamode)
       name = "alamode file"
    case (isformat_r_castepphonon)
       name = "CASTEP phonon file"
    case (isformat_r_akaikkr)
       name = "AkaiKKR input"
    case (isformat_r_xband)
       name = "xband sysfile"
    case default
       name = "unknown"
    end select

  end function format_name

end submodule proc
