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

! Driver routines for structural calculations
submodule (struct_drivers) proc
  implicit none

  !xx! private routines
  ! subroutine struct_comparevc_vcpwdf(s,line)
  ! subroutine struct_comparevc_vcgpwdf(s,line)

contains

  !xx! top-level routines
  !> Parse the input of the crystal keyword (line) and return an
  !> initialized crystal c. mol0=1, interpret the structure as a
  !> molecule; mol0=0, a crystal, mol0=-1, guess. If allownofile
  !> allow the program to read the following input lines to parse
  !> the crystal/molecule environment.
  module subroutine struct_crystal_input(line,mol0,allownofile,verbose,s0,cr0,seed0)
    use systemmod, only: system
    use crystalmod, only: crystal
    use param, only: isformat_cif, isformat_shelx, isformat_f21,&
       isformat_cube, isformat_bincube, isformat_struct, isformat_abinit,&
       isformat_elk, isformat_fploout,&
       isformat_qein, isformat_qeout, isformat_crystal, isformat_xyz, isformat_gjf,&
       isformat_wfn, isformat_wfx, isformat_fchk, isformat_molden,&
       isformat_gaussian, isformat_siesta, isformat_xsf, isformat_gen,&
       isformat_vasp, isformat_pwc, isformat_axsf, isformat_dat,&
       isformat_pgout, isformat_orca, isformat_dmain, isformat_aimsin,&
       isformat_aimsout, isformat_tinkerfrac, isformat_castepcell,&
       isformat_castepgeom, isformat_mol2, isformat_pdb, isformat_zmat,&
       isformat_sdf, isformat_unknown
    use crystalseedmod, only: crystalseed, struct_detect_format,&
       struct_detect_ismol
    use global, only: doguess, iunit, dunit0, rborder_def, eval_next
    use tools_io, only: getword, equal, ferror, faterr, zatguess, lgetword,&
       string, uin, isinteger, isreal, lower
    use types, only: realloc
    character*(*), intent(in) :: line
    integer, intent(in) :: mol0
    logical, intent(in) :: allownofile
    logical, intent(in) :: verbose
    type(system), intent(inout), optional :: s0
    type(crystal), intent(inout), optional :: cr0
    type(crystalseed), intent(inout), optional :: seed0

    integer :: lp, lp2, istruct
    character(len=:), allocatable :: word, word2, lword, subline
    integer :: nn, isformat, id
    real*8 :: rborder, raux, xnudge
    logical :: docube, ok, ismol, mol, hastypes, readtypes
    type(crystalseed) :: seed
    character(len=:), allocatable :: errmsg

    ! read and parse
    lp=1
    word = getword(line,lp)

    ! detect the format for this file; see if we the user is forcing a particular format
    lword = lower(word)
    isformat = isformat_unknown
    if (equal(lword,'cif')) then
       isformat = isformat_cif
    elseif (equal(lword,'shelx')) then
       isformat = isformat_shelx
    elseif (equal(lword,'21')) then
       isformat = isformat_f21
    elseif (equal(lword,'cube')) then
       isformat = isformat_cube
    elseif (equal(lword,'bincube')) then
       isformat = isformat_bincube
    elseif (equal(lword,'wien')) then
       isformat = isformat_struct
    elseif (equal(lword,'abinit')) then
       isformat = isformat_abinit
    elseif (equal(lword,'elk')) then
       isformat = isformat_elk
    elseif (equal(lword,'qe_in')) then
       isformat = isformat_qein
    elseif (equal(lword,'qe_out')) then
       isformat = isformat_qeout
    elseif (equal(lword,'crystal')) then
       isformat = isformat_crystal
    elseif (equal(lword,'fplo')) then
       isformat = isformat_fploout
    elseif (equal(lword,'xyz')) then
       isformat = isformat_xyz
    elseif (equal(lword,'gjf').or.equal(lword,'com')) then
       isformat = isformat_gjf
    elseif (equal(lword,'zmat')) then
       isformat = isformat_zmat
    elseif (equal(lword,'wfn')) then
       isformat = isformat_wfn
    elseif (equal(lword,'wfx')) then
       isformat = isformat_wfx
    elseif (equal(lword,'fchk')) then
       isformat = isformat_fchk
    elseif (equal(lword,'molden')) then
       isformat = isformat_molden
    elseif (equal(lword,'gaussian')) then
       isformat = isformat_gaussian
    elseif (equal(lword,'siesta')) then
       isformat = isformat_siesta
    elseif (equal(lword,'cell')) then
       isformat = isformat_castepcell
    elseif (equal(lword,'geom')) then
       isformat = isformat_castepgeom
    elseif (equal(lword,'xsf')) then
       isformat = isformat_xsf
    elseif (equal(lword,'gen')) then
       isformat = isformat_gen
    elseif (equal(lword,'vasp')) then
       isformat = isformat_vasp
    elseif (equal(lword,'pwc')) then
       isformat = isformat_pwc
    elseif (equal(lword,'axsf')) then
       isformat = isformat_axsf
    elseif (equal(lword,'dat')) then
       isformat = isformat_dat
    elseif (equal(lword,'pgout')) then
       isformat = isformat_pgout
    elseif (equal(lword,'orca')) then
       isformat = isformat_orca
    elseif (equal(lword,'dmain')) then
       isformat = isformat_dmain
    elseif (equal(lword,'fhiaims_in')) then
       isformat = isformat_aimsin
    elseif (equal(lword,'fhiaims_out')) then
       isformat = isformat_aimsout
    elseif (equal(lword,'frac')) then
       isformat = isformat_tinkerfrac
    elseif (equal(lword,'mol2')) then
       isformat = isformat_mol2
    elseif (equal(lword,'pdb')) then
       isformat = isformat_pdb
    elseif (equal(lword,'sdf')) then
       isformat = isformat_sdf
    end if
    if (isformat /= isformat_unknown) then
       word = getword(line,lp)
    else
       call struct_detect_format(word,isformat)
    end if
    call struct_detect_ismol(word,isformat,ismol)
    subline = line(lp:)
    word2 = getword(line,lp)
    if (len_trim(word2) == 0) word2 = " "

    ! is this a crystal or a molecule?
    if (mol0 == 1) then
       mol = .true.
    elseif (mol0 == 0) then
       mol = .false.
    elseif (mol0 == -1) then
       mol = ismol
    else
       call ferror("struct_crystal_input","unknown mol0",faterr)
    end if
    docube = .false.
    rborder = rborder_def
    id = 0
    lp2=1
    if (ismol) then
       do while(.true.)
          if (equal(word2,'cubic').or.equal(word2,'cube')) then
             docube = .true.
          elseif (equal(word2,'id')) then
             word2 = lgetword(line,lp)
             if (.not.isinteger(id,word2)) then
                call ferror('struct_crystal_input','Wrong ID in sdf/mol input',faterr,line,syntax=.true.)
                return
             end if
          elseif (eval_next(raux,word2,lp2)) then
             rborder = raux / dunit0(iunit)
          elseif (len_trim(word2) > 1) then
             if (isformat == isformat_mol2) then
                exit
             else
                call ferror('struct_crystal_input','Unknown extra keyword',faterr,line,syntax=.true.)
                return
             end if
          else
             exit
          end if
          lp2 = 1
          word2 = lgetword(line,lp)
       end do
    end if

    ! build the seed
    errmsg = ""
    if (isformat == isformat_cif) then
       call seed%read_cif(word,word2,mol,errmsg)

    elseif (isformat == isformat_mol2) then
       call seed%read_mol2(word,rborder,docube,word2,errmsg)

    elseif (isformat == isformat_sdf) then
       call seed%read_sdf(word,rborder,docube,errmsg,id)

    elseif (isformat == isformat_shelx) then
       call seed%read_shelx(word,mol,errmsg)

    elseif (isformat == isformat_f21) then
       call seed%read_f21(word,mol,errmsg)

    elseif (isformat == isformat_cube) then
       call seed%read_cube(word,mol,errmsg)

    elseif (isformat == isformat_bincube) then
       call seed%read_bincube(word,mol,errmsg)

    elseif (isformat == isformat_struct) then
       call seed%read_wien(word,mol,errmsg)

    elseif (isformat == isformat_vasp) then
       call seed%read_vasp(word,mol,hastypes,errmsg)
       if (.not.hastypes) then
          if (index(word2,'POTCAR') > 0) then
             call seed%read_potcar(word2,errmsg)
             readtypes = .true.
          else
             seed%nspc = 0
             allocate(seed%spc(2))
             do while(.true.)
                nn = zatguess(word2)
                if (nn >= 0) then
                   seed%nspc = seed%nspc + 1
                   if (seed%nspc > size(seed%spc,1)) &
                      call realloc(seed%spc,2*seed%nspc)
                   seed%spc(seed%nspc)%name = string(word2)
                   seed%spc(seed%nspc)%z = zatguess(word2)
                else
                   if (len_trim(word2) > 0) then
                      call ferror('struct_crystal_input','Unknown atom type in CRYSTAL',faterr,line,syntax=.true.)
                      return
                   end if
                   exit
                end if
                word2 = getword(line,lp)
             end do
          end if
          if (seed%nspc == 0) &
             errmsg = "Atom types not found (use POTCAR or give atom types after structure file)"
          call realloc(seed%spc,seed%nspc)
       end if

    elseif (isformat == isformat_abinit) then
       call seed%read_abinit(word,mol,errmsg)

    elseif (isformat == isformat_elk) then
       call seed%read_elk(word,mol,errmsg)

    elseif (isformat == isformat_qeout) then
       lp2 = 1
       ok = isinteger(istruct,word2,lp2)
       if (.not.ok) istruct = 0
       call seed%read_qeout(word,mol,istruct,errmsg)

    elseif (isformat == isformat_crystal) then
       call seed%read_crystalout(word,mol,errmsg)

    elseif (isformat == isformat_fploout) then
       call seed%read_fploout(word,mol,errmsg)

    elseif (isformat == isformat_qein) then
       call seed%read_qein(word,mol,errmsg)

    elseif (isformat == isformat_xyz.or.isformat == isformat_wfn.or.&
       isformat == isformat_wfx.or.isformat == isformat_fchk.or.&
       isformat == isformat_molden.or.isformat == isformat_gaussian.or.&
       isformat == isformat_dat.or.isformat == isformat_pgout.or.&
       isformat == isformat_orca.or.isformat == isformat_gjf.or.&
       isformat == isformat_pdb.or.isformat == isformat_zmat) then
       call seed%read_mol(word,isformat,rborder,docube,errmsg)

    elseif (isformat == isformat_siesta) then
       call seed%read_siesta(word,mol,errmsg)

    elseif (isformat == isformat_castepcell) then
       call seed%read_castep_cell(word,mol,errmsg)

    elseif (isformat == isformat_castepgeom) then
       call seed%read_castep_geom(word,mol,errmsg)

    elseif (isformat == isformat_xsf) then
       call seed%read_xsf(word,rborder,docube,errmsg)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_gen) then
       call seed%read_dftbp(word,rborder,docube,errmsg)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_pwc) then
       call seed%read_pwc(word,mol,errmsg)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_axsf) then
       istruct = 1
       xnudge = 0d0

       lp2 = 1
       ok = isinteger(nn,word2,lp2)
       if (ok) then
          istruct = nn
          word2 = getword(line,lp)
          lp2 = 1
          ok = isreal(raux,word2,lp2)
          if (ok) xnudge = raux
       end if

       call seed%read_axsf(word,istruct,xnudge,rborder,docube,errmsg)
       if (mol0 /= -1) &
          seed%ismolecule = mol

    elseif (isformat == isformat_dmain) then
       call seed%read_dmain(word,mol,errmsg)

    elseif (isformat == isformat_aimsin) then
       call seed%read_aimsin(word,mol,rborder,docube,errmsg)

    elseif (isformat == isformat_aimsout) then
       call seed%read_aimsout(word,mol,rborder,docube,errmsg)

    elseif (isformat == isformat_tinkerfrac) then
       call seed%read_tinkerfrac(word,mol,errmsg)

    else if (equal(lower(word),'library')) then
       call seed%read_library(subline,ok,mol=mol)
       if (.not.ok) return

    else if (len_trim(word) < 1) then
       if (.not.allownofile) then
          call ferror('struct_crystal_input','Attempted to parse CRYSTAL/MOLECULE but allownofile=.true.',faterr,syntax=.true.)
          return
       end if

       if (.not.mol) then
          call seed%parse_crystal_env(uin,ok)
          if (.not.ok) return
       else
          call seed%parse_molecule_env(uin,ok)
          if (.not.ok) return
       endif
    else
       call ferror('struct_crystal_input','unrecognized file format',faterr,line,syntax=.true.)
       return
    end if
    if (len_trim(errmsg) > 0) then
       call ferror('struct_crystal_input',errmsg,faterr,line,syntax=.true.)
       return
    end if

    ! handle the doguess option
    if (.not.seed%ismolecule) then
       if (doguess == 0) then
          seed%findsym = 0
       elseif (doguess == 1 .and. seed%havesym == 0) then
          seed%findsym = 1
       end if
    end if

    if (present(s0)) then
       call s0%new_from_seed(seed)
       if (verbose) &
          call s0%report(.true.,.true.,.true.,.true.,.true.,.true.,.false.)
    end if
    if (present(cr0)) then
       call cr0%struct_new(seed,.true.)
       if (verbose) &
          call cr0%report(.true.,.true.)
    end if
    if (present(seed0)) then
       seed0 = seed
    end if

  end subroutine struct_crystal_input

  !> Clear the symmetry in the system.
  module subroutine struct_sym(s,line,verbose)
    use iso_c_binding, only: c_double
    use crystalseedmod, only: crystalseed
    use spglib, only: SpglibDataset, spg_standardize_cell
    use global, only: symprec
    use tools_math, only: matinv
    use tools_io, only: uout, lgetword, equal, isinteger, ferror, faterr, string
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line
    logical, intent(in) :: verbose

    character(len=:), allocatable :: word, errmsg
    integer :: lp, i
    real*8 :: osp
    type(SpglibDataset) :: spg
    real*8 :: x0(3,3)
    logical :: isempty

    real*8, parameter :: spmin = 1d-10
    real*8, parameter :: factor = 10d0

    ! header
    write (uout,'("* SYM: manipulation of the crystal symmetry.")')
    write (uout,*)

    ! no molecules here
    if (s%c%ismolecule) then
       call ferror('struct_sym','Symmetry not available in a MOLECULE',faterr,line,syntax=.true.)
       return
    end if

    ! parse the input
    isempty = .false.
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,'clear')) then
       write (uout,'("+ CLEAR the symmetry for this crystal structure.")')
       call s%clearsym()

    elseif (equal(word,'recalc')) then
       write (uout,'("+ RECALCULATE the symmetry operations.")')
       call s%c%calcsym(.false.,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror("struct_sym","spglib: "//errmsg,faterr)

    elseif (equal(word,'analysis')) then
       write (uout,'("+ ANALYSIS of the crystal symmetry.")')
       write (uout,*)
       osp = symprec

       write (uout,'("# Sym.Prec. Space group")')
       symprec = spmin / factor
       do i = 1, 11
          symprec = symprec * factor
          call s%c%spglib_wrap(spg,.false.,errmsg)
          if (len_trim(errmsg) > 0) exit
          write (uout,'(1p,"  ",A,A," ","(",A,")")') string(symprec,'e',10,2), string(spg%international_symbol), &
             string(spg%spacegroup_number)
       end do
       write (uout,*)
       symprec = osp
       return ! do not clear the fields or re-write the structure

    elseif (equal(word,'refine')) then
       x0 = s%c%cell_standard(.false.,.false.,.true.)

    elseif (equal(word,'wholemols')) then
       if (.not.s%c%ismol3d) then
          call ferror('struct_sym','Cannot apply WHOLEMOLS to a non-molecular crystal',faterr,line,syntax=.true.)
          return
       end if
       call s%c%wholemols()

    else
       isempty = .true.
       call s%c%struct_report_symmetry()
    end if

    ! clear the fields and report the new structure
    if (.not.isempty) then
       call s%reset_fields()
       if (verbose) &
          call s%report(.true.,.true.,.true.,.true.,.true.,.true.,.false.)
    end if

  end subroutine struct_sym

  !> Change the charges and pseudopotential charges in the system.
  module subroutine struct_charges(s,line,oksyn)
    use systemmod, only: system
    use grid1mod, only: grid1_register_core
    use global, only: eval_next
    use tools_io, only: lgetword, equal, ferror, faterr, getword, zatguess
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line
    logical, intent(out) :: oksyn

    character(len=:), allocatable :: word
    integer :: lp, nn, i
    logical :: ok
    real*8 :: xx

    oksyn = .false.
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,'q') .or. equal(word,'qat')) then
       do while (.true.)
          word = getword(line,lp)
          if (len_trim(word) < 1) exit
          nn = zatguess(word)
          if (nn == -1) then
             call ferror('struct_charges','Unknown atomic symbol in Q/QAT',faterr,line,syntax=.true.)
             return
          end if
          ok = eval_next(xx,line,lp)
          if (.not.ok) then
             call ferror('struct_charges','Incorrect Q/QAT syntax',faterr,line,syntax=.true.)
             return
          end if
          do i = 1, s%c%nspc
             if (s%c%spc(i)%z == nn) &
                s%c%spc(i)%qat = xx
          end do
       end do
    endif
    oksyn = .true.

    ! report the charges
    call s%c%report(.false.,.true.)
    call s%report(.false.,.false.,.false.,.false.,.false.,.true.,.false.)

  end subroutine struct_charges

  ! Write the crystal structure to a file. If usexyznames, use the
  ! atom names when writing an xyz file.
  module subroutine struct_write(s,line,usexyznames)
    use systemmod, only: system
    use global, only: eval_next, dunit0, iunit
    use tools_io, only: getword, equal, lower, lgetword, ferror, faterr, uout, &
       string
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line
    logical, intent(in) :: usexyznames

    character(len=:), allocatable :: word, wext, file, wroot
    integer :: lp, ix(3), lp2, iaux, nmer, idx
    logical :: doborder, molmotif, dosym, docell, domolcell, ok, isqe
    logical :: onemotif, environ, lnmer, doexternal
    real*8 :: rsph, xsph(3), rcub, xcub(3), renv, rk

    lp = 1
    file = getword(line,lp)
    wext = lower(file(index(file,'.',.true.)+1:))
    wroot = file(:index(file,'.',.true.)-1)

    if (equal(wext,'xyz').or.equal(wext,'gjf').or.equal(wext,'cml').or.&
        equal(wext,'obj').or.equal(wext,'ply').or.equal(wext,'off')) then
       ! xyz, gjf, cml, obj, ply, off
       doborder = .false.
       lnmer = .false.
       onemotif = .false.
       molmotif = .false.
       environ = .false.
       docell = .false.
       domolcell = .false.
       nmer = 1
       ix = 1
       rsph = -1d0
       xsph = 0d0
       rcub = -1d0
       xcub = 0d0
       renv = 0d0
       do while(.true.)
          word = lgetword(line,lp)
          lp2 = 1
          if (equal(word,'border')) then
             doborder = .true.
          elseif (equal(word,'onemotif')) then
             onemotif = .true.
          elseif (equal(word,'molmotif')) then
             molmotif = .true.
          elseif (equal(word,'environ')) then
             environ = .true.
             ok = eval_next(renv,line,lp)
             if (.not.ok) &
                call ferror('struct_write','incorrect ENVIRON',faterr,line,syntax=.true.)
             renv = renv / dunit0(iunit)
          elseif (equal(word,'nmer')) then
             lnmer = .true.
             ok = eval_next(nmer,line,lp)
             if (.not.ok) &
                call ferror('struct_write','incorrect NMER',faterr,line,syntax=.true.)
          elseif (equal(word,'cell')) then
             docell = .true.
          elseif (equal(word,'molcell')) then
             domolcell = .true.
          elseif (equal(word,'sphere')) then
             ok = eval_next(rsph,line,lp)
             ok = ok .and. eval_next(xsph(1),line,lp)
             if (ok) then
                ok = ok .and. eval_next(xsph(2),line,lp)
                ok = ok .and. eval_next(xsph(3),line,lp)
                if (.not.ok) then
                   call ferror('struct_write','incorrect WRITE SPHERE syntax',faterr,line,syntax=.true.)
                   return
                end if
             else
                xsph = 0d0
             end if
             ! convert units and coordinates
             rsph = rsph / dunit0(iunit)
             if (s%c%ismolecule) then
                xsph = xsph / dunit0(iunit) - s%c%molx0
                xsph = s%c%c2x(xsph)
             endif
          elseif (equal(word,'cube')) then
             ok = eval_next(rcub,line,lp)
             ok = ok .and. eval_next(xcub(1),line,lp)
             if (ok) then
                ok = ok .and. eval_next(xcub(2),line,lp)
                ok = ok .and. eval_next(xcub(3),line,lp)
                if (.not.ok) then
                   call ferror('struct_write','incorrect WRITE CUBE syntax',faterr,line,syntax=.true.)
                   return
                end if
             else
                xcub = 0d0
             end if
             ! convert units and coordinates
             rcub = rcub / dunit0(iunit)
             if (s%c%ismolecule) then
                xcub = xcub / dunit0(iunit) - s%c%molx0
                xcub = s%c%c2x(xcub)
             endif
          elseif (eval_next(iaux,word,lp2)) then
             ix(1) = iaux
             ok = eval_next(ix(2),line,lp)
             ok = ok .and. eval_next(ix(3),line,lp)
             if (.not.ok) then
                call ferror('struct_write','incorrect WRITE syntax',faterr,line,syntax=.true.)
                return
             end if
          elseif (len_trim(word) > 1) then
             call ferror('struct_write','Unknown extra keyword',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do

       if (.not.lnmer) then
          write (uout,'("* WRITE ",A," file: ",A)') trim(wext), string(file)
       else
          write (uout,'("* WRITE ",A," files: ",A,"_*.",A)') trim(wext), &
             string(wroot), trim(wext)
       end if

       if (equal(wext,'xyz').or.equal(wext,'gjf').or.equal(wext,'cml')) then
          call s%c%write_mol(file,wext,ix,doborder,onemotif,molmotif,&
             environ,renv,lnmer,nmer,rsph,xsph,rcub,xcub,usexyznames)
       else
          call s%c%write_3dmodel(file,wext,ix,doborder,onemotif,molmotif,&
             docell,domolcell,rsph,xsph,rcub,xcub)
       end if
    elseif (equal(wext,'gau')) then
       ! gaussian periodic boundary conditions
       write (uout,'("* WRITE Gaussian file: ",A)') string(file)
       call s%c%write_gaussian(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'in')) then
       ! quantum espresso (scf.in) and FHI aims (geometry.in)
       idx = index(wroot,'.',.true.)
       isqe = .false.
       if (idx > 0) then
          wext = lower(wroot(idx+1:))
          isqe = equal(wext,'scf')
       endif

       ok = eval_next(rk,line,lp)

       ! espresso
       if (isqe) then
          write (uout,'("* WRITE espresso file: ",A)') string(file)
          if (ok) then
             call s%c%write_espresso(file,rk)
          else
             call s%c%write_espresso(file)
          end if
       else
          write (uout,'("* WRITE FHIaims file: ",A)') string(file)
          if (ok) then
             call s%c%write_fhi(file,.true.,rk)
          else
             call s%c%write_fhi(file,.true.)
          end if
       end if

    elseif (equal(wext,'pwi')) then
       ! quantum espresso input
       ok = eval_next(rk,line,lp)
       write (uout,'("* WRITE espresso file: ",A)') string(file)
       if (ok) then
          call s%c%write_espresso(file,rk)
       else
          call s%c%write_espresso(file)
       end if
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'cell')) then
       ! CASTEP cell
       ok = eval_next(rk,line,lp)
       write (uout,'("* WRITE CASTEP cell file: ",A)') string(file)
       if (ok) then
          call s%c%write_castep_cell(file,rk)
       else
          call s%c%write_castep_cell(file)
       end if
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'poscar') .or. equal(wext,'contcar')) then
       ! vasp
       write (uout,'("* WRITE VASP file: ",A)') string(file)
       call s%c%write_vasp(file,.true.)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'abin')) then
       ! abinit
       write (uout,'("* WRITE abinit file: ",A)') string(file)
       call s%c%write_abinit(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'elk')) then
       ! elk
       write (uout,'("* WRITE elk file: ",A)') string(file)
       call s%c%write_elk(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'tess')) then
       ! tessel
       write (uout,'("* WRITE tess file: ",A)') string(file)
       call s%c%write_tessel(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'incritic').or.equal(wext,'cri')) then
       ! critic2
       write (uout,'("* WRITE critic2 input file: ",A)') string(file)
       call s%c%write_critic(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'cif') .or. equal(wext,'d12') .or. equal(wext,'res')) then
       ! cif and d12
       dosym = .true.
       doexternal = .true.
       do while(.true.)
          word = lgetword(line,lp)
          if (equal(word,'nosym').or.equal(word,'nosymm')) then
             dosym = .false.
          elseif (equal(word,'noexternal')) then
             doexternal = .false.
          elseif (len_trim(word) > 1) then
             call ferror('struct_write','Unknown extra keyword',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do

       if (equal(wext,'cif')) then
          write (uout,'("* WRITE cif file: ",A)') string(file)
          call s%c%write_cif(file,dosym)
       elseif (equal(wext,'d12')) then
          write (uout,'("* WRITE crystal file: ",A)') string(file)
          if (doexternal) &
             write (uout,'(/"+ WRITE fort.34 file: ",A)') file(:index(file,'.',.true.)-1) // ".fort.34"
          call s%c%write_d12(file,dosym,doexternal)
       elseif (equal(wext,'res')) then
          write (uout,'("* WRITE res file: ",A)') string(file)
          if (dosym) then
             call s%c%write_res(file,1)
          else
             call s%c%write_res(file,0)
          end if
       end if
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'m')) then
       ! escher
       write (uout,'("* WRITE escher file: ",A)') string(file)
       call s%c%write_escher(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'db')) then
       ! db
       write (uout,'("* WRITE db file: ",A)') string(file)
       call s%c%write_db(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'gin')) then
       ! gulp
       write (uout,'("* WRITE gulp file: ",A)') string(file)
       call s%c%write_gulp(file)
    elseif (equal(wext,'lammps')) then
       write (uout,'("* WRITE lammps file: ",A)') string(file)
       call s%c%write_lammps(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'fdf')) then
       write (uout,'("* WRITE fdf file: ",A)') string(file)
       call s%c%write_siesta_fdf(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'struct_in')) then
       write (uout,'("* WRITE STRUCT_IN file: ",A)') string(file)
       call s%c%write_siesta_in(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'hsd')) then
       write (uout,'("* WRITE hsd file: ",A)') string(file)
       call s%c%write_dftbp_hsd(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'gen')) then
       write (uout,'("* WRITE gen file: ",A)') string(file)
       call s%c%write_dftbp_gen(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'pyscf')) then
       ! pyscf (python script)
       write (uout,'("* WRITE pyscf file: ",A)') string(file)
       call s%c%write_pyscf(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'frac')) then
       ! abinit
       write (uout,'("* WRITE TINKER frac file: ",A)') string(file)
       call s%c%write_tinkerfrac(file)
       ok = check_no_extra_word()
       if (.not.ok) return
    else
       call ferror('struct_write','unrecognized file format',faterr,line,syntax=.true.)
       return
    end if
    write (uout,*)

  contains
    function check_no_extra_word()
      character(len=:), allocatable :: aux2
      logical :: check_no_extra_word

      check_no_extra_word = .true.
      do while (.true.)
         aux2 = getword(line,lp)
         if (equal(aux2,'primitive')) cycle
         if (len_trim(aux2) > 0) then
            call ferror('struct_write','Unknown extra keyword',faterr,line,syntax=.true.)
            check_no_extra_word = .false.
         end if
         exit
      end do

    end function check_no_extra_word
  end subroutine struct_write

  !> Relabel atoms based on user's input
  module subroutine struct_atomlabel(s,line)
    use systemmod, only: system
    use global, only: iunitname0, dunit0, iunit
    use tools_io, only: tab, string, nameguess, lower, ioj_center, uout, equal
    use types, only: species, realloc
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line

    character(len=:), allocatable :: templ, aux, aux2

    integer :: i, j, inum, idx, iz
    integer :: nspc
    integer, allocatable :: is(:)
    type(species), allocatable :: spc(:)
    logical :: found

    nspc = 0
    allocate(spc(1))
    allocate(is(s%c%nneq))

    ! clean up the input label and build the template
    templ = trim(adjustl(line))
    aux = ""
    do i = 1, len(templ)
       if (templ(i:i) == "'" .or. templ(i:i) == '"' .or. templ(i:i) == tab) cycle
       aux2 = trim(aux) // templ(i:i)
       aux = aux2
    end do
    templ = aux

    do i = 1, s%c%nneq
       iz = s%c%spc(s%c%at(i)%is)%z
       aux = trim(templ)
       do while(.true.)
          if (index(aux,"%aid") > 0) then
             ! the absolute index for this atom
             idx = index(aux,"%aid")
             aux2 = aux(1:idx-1) // string(i) // aux(idx+4:)
             aux = trim(aux2)
          elseif (index(aux,"%id") > 0) then
             ! the index by counting atoms only of this type
             inum = 0
             do j = 1, i
                if (s%c%spc(s%c%at(j)%is)%z == iz) inum = inum + 1
             end do
             idx = index(aux,"%id")
             aux2 = aux(1:idx-1) // string(inum) // aux(idx+3:)
             aux = trim(aux2)
          elseif (index(aux,"%S") > 0) then
             ! the atomic symbol
             idx = index(aux,"%S")
             aux2 = aux(1:idx-1) // string(nameguess(iz,.true.)) // aux(idx+2:)
             aux = trim(aux2)
          elseif (index(aux,"%s") > 0) then
             ! the atomic symbol, lowercase
             idx = index(aux,"%s")
             aux2 = aux(1:idx-1) // string(lower(nameguess(iz,.true.))) // aux(idx+2:)
             aux = trim(aux2)
          elseif (index(aux,"%l") > 0) then
             ! the original atom label
             idx = index(aux,"%l")
             aux2 = aux(1:idx-1) // string(s%c%spc(s%c%at(i)%is)%name) // aux(idx+2:)
             aux = trim(aux2)
          else
             exit
          endif
       end do

       found = .false.
       do j = 1, nspc
          if (equal(aux,spc(j)%name) .and. s%c%spc(s%c%at(i)%is)%z == spc(j)%z) then
             found = .true.
             is(i) = j
             exit
          end if
       end do

       if (.not.found) then
          nspc = nspc + 1
          if (nspc > size(spc,1)) &
             call realloc(spc,2*nspc)
          spc(nspc)%name = trim(aux)
          spc(nspc)%z = s%c%spc(s%c%at(i)%is)%z
          spc(nspc)%qat = s%c%spc(s%c%at(i)%is)%qat
          is(i) = nspc
       end if
    end do
    call realloc(spc,nspc)

    s%c%nspc = nspc
    if (allocated(s%c%spc)) deallocate(s%c%spc)
    allocate(s%c%spc(s%c%nspc))
    s%c%spc = spc
    do i = 1, s%c%nneq
       s%c%at(i)%is = is(i)
    end do
    if (allocated(s%c%atcel)) then
       do i = 1, s%c%ncel
          s%c%atcel(i)%is = s%c%at(s%c%atcel(i)%idx)%is
       end do
    end if
    deallocate(is)

    ! Write the list of atomic coordinates
    if (.not.s%c%ismolecule) then
       write (uout,'("+ List of non-equivalent atoms in the unit cell (cryst. coords.): ")')
       write (uout,'("# ",7(A," "))') string("nat",3,ioj_center), &
          string("x",14,ioj_center), string("y",14,ioj_center),&
          string("z",14,ioj_center), string("name",10,ioj_center), &
          string("mult",4,ioj_center), string("Z",4,ioj_center)
       do i=1, s%c%nneq
          write (uout,'("  ",7(A," "))') &
             string(i,3,ioj_center),&
             string(s%c%at(i)%x(1),'f',length=14,decimal=10,justify=3),&
             string(s%c%at(i)%x(2),'f',length=14,decimal=10,justify=3),&
             string(s%c%at(i)%x(3),'f',length=14,decimal=10,justify=3),&
             string(s%c%spc(s%c%at(i)%is)%name,10,ioj_center), &
             string(s%c%at(i)%mult,4,ioj_center), &
             string(s%c%spc(s%c%at(i)%is)%z,4,ioj_center)
       enddo
       write (uout,*)
    else
       write (uout,'("+ List of atoms in Cartesian coordinates (",A,"): ")') iunitname0(iunit)
       write (uout,'("# ",6(A," "))') string("at",3,ioj_center), &
          string("x",16,ioj_center), string("y",16,ioj_center),&
          string("z",16,ioj_center), string("name",10,ioj_center),&
          string("Z",4,ioj_center)
       do i=1,s%c%ncel
          write (uout,'("  ",6(A," "))') &
             string(i,3,ioj_center),&
             (string((s%c%atcel(i)%r(j)+s%c%molx0(j))*dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3),&
             string(s%c%spc(s%c%atcel(i)%is)%name,10,ioj_center),&
             string(s%c%spc(s%c%atcel(i)%is)%z,4,ioj_center)
       enddo
       write (uout,*)
    end if

  end subroutine struct_atomlabel

  !> Calculate the powder diffraction pattern for the current
  !structure.
  module subroutine struct_powder(s,line)
    use crystalmod, only: xrpd_lambda_def, xrpd_fpol_def, xrpd_sigma_def
    use systemmod, only: system
    use global, only: fileroot, eval_next
    use tools_io, only: ferror, faterr, uout, lgetword, equal, getword, &
       fopen_write, string, ioj_center, string, fclose
    use param, only: pi
    type(system), intent(in) :: s
    character*(*), intent(in) :: line

    integer :: i, lp, lu, np, npts
    real*8 :: th2ini, th2end, lambda, fpol, sigma
    real*8, allocatable :: t(:), ih(:), th2p(:), ip(:)
    integer, allocatable :: hvecp(:,:)
    character(len=:), allocatable :: root, word
    logical :: ok, ishard

    if (s%c%ismolecule) then
       call ferror("struct_powder","POWDER cannot be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! default values
    th2ini = 5d0
    th2end = 90d0
    lambda = xrpd_lambda_def
    fpol = xrpd_fpol_def
    npts = 10001
    sigma = xrpd_sigma_def
    root = trim(fileroot) // "_xrd"
    ishard = .false.

    ! header
    write (uout,'("* POWDER: powder diffraction pattern")')
    write (uout,*)

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"th2ini")) then
          ok = eval_next(th2ini,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect TH2INI',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"hard")) then
          ishard = .true.
       elseif (equal(word,"soft")) then
          ishard = .false.
       elseif (equal(word,"th2end")) then
          ok = eval_next(th2end,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect TH2END',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"l").or.equal(word,"lambda")) then
          ok = eval_next(lambda,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect LAMBDA',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"fpol")) then
          ok = eval_next(fpol,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect FPOL',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"npts")) then
          ok = eval_next(npts,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect NPTS',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"sigma")) then
          ok = eval_next(sigma,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect SIGMA',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"root")) then
          root = getword(line,lp)
       elseif (len_trim(word) > 0) then
          call ferror('struct_powder','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    call s%c%powder(0,th2ini,th2end,lambda,fpol,npts=npts,sigma=sigma,ishard=ishard,&
       t=t,ih=ih,th2p=th2p,ip=ip,hvecp=hvecp)
    np = size(th2p)

    ! write the data file
    lu = fopen_write(trim(root) // ".dat")
    write (lu,'("# ",A,A)') string("2*theta",13,ioj_center), string("Intensity",13,ioj_center)
    do i = 1, npts
       write (lu,'(A," ",A)') string(t(i),"f",15,7,ioj_center), &
          string(ih(i),"f",15,7,ioj_center)
    end do
    call fclose(lu)
    deallocate(t,ih)

    ! write the gnuplot input file
    lu = fopen_write(trim(root) // ".gnu")
    write (lu,'("set terminal postscript eps color enhanced ""Helvetica"" 25")')
    write (lu,'("set output """,A,".eps""")') trim(root)
    write (lu,*)
    do i = 1, np
       write (lu,'("set label ",A," """,A,A,A,""" at ",A,",",A," center rotate by 90 font ""Helvetica,12""")') &
          string(i), string(hvecp(1,i)), string(hvecp(2,i)), string(hvecp(3,i)), &
          string(th2p(i)*180/pi,"f"), string(min(ip(i)+4,102d0),"f")
    end do
    write (lu,*)
    write (lu,'("set xlabel ""2{/Symbol Q} (degrees)""")')
    write (lu,'("set ylabel ""Intensity (arb. units)""")')
    write (lu,'("set xrange [",A,":",A,"]")') string(th2ini,"f"), string(th2end,"f")
    write (lu,'("set style data lines")')
    write (lu,'("set grid")')
    write (lu,'("unset key")')
    write (lu,'("plot """,A,".dat"" w lines")') string(root)
    write (lu,*)
    call fclose(lu)

    ! list of peaks and intensity in output
    write(uout,'("#i  h  k  l       2*theta        Intensity ")')
    do i = 1, np
       write (uout,'(99(A," "))') string(i,3), string(hvecp(1,i),2), string(hvecp(2,i),2),&
          string(hvecp(3,i),2), string(th2p(i)*180/pi,"f",15,7,7), string(ip(i),"f",15,7,7)
    end do
    write (uout,*)

    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    if (allocated(th2p)) deallocate(th2p)
    if (allocated(ip)) deallocate(ip)
    if (allocated(hvecp)) deallocate(hvecp)

  end subroutine struct_powder

  !> Calculate the radial distribution function for the current
  !> structure.
  module subroutine struct_rdf(s,line)
    use systemmod, only: system
    use global, only: fileroot, eval_next, dunit0, iunit, iunitname0, iunit_bohr
    use tools_io, only: faterr, ferror, uout, lgetword, equal, fopen_write,&
       ioj_center, getword, string, fclose, isinteger
    use types, only: realloc
    type(system), intent(inout), target :: s
    character*(*), intent(in) :: line

    real*8 :: rini, rend
    character(len=:), allocatable :: root, word
    logical :: ok, ishard
    integer :: npts
    real*8, allocatable :: t(:), ih(:)
    integer, allocatable :: ipairs(:,:)
    integer :: lp, lu, i, npairs
    real*8 :: sigma

    ! default values
    rini = 0d0
    rend = 25d0
    root = trim(fileroot) // "_rdf"
    npts = 10001
    sigma = 0.05d0
    npairs = 0
    allocate(ipairs(2,10))
    ishard = .false.

    ! header
    write (uout,'("* RDF: radial distribution function")')
    write (uout,*)

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"rend")) then
          ok = eval_next(rend,line,lp)
          if (.not.ok) then
             call ferror('struct_rdf','Incorrect REND',faterr,line,syntax=.true.)
             return
          end if
          rend = rend / dunit0(iunit)
       elseif (equal(word,"rini")) then
          ok = eval_next(rini,line,lp)
          if (.not.ok) then
             call ferror('struct_rdf','Incorrect RINI',faterr,line,syntax=.true.)
             return
          end if
          rini = rini / dunit0(iunit)
       elseif (equal(word,"sigma")) then
          ok = eval_next(sigma,line,lp)
          if (.not.ok) then
             call ferror('struct_rdf','Incorrect SIGMA',faterr,line,syntax=.true.)
             return
          end if
          sigma = sigma / dunit0(iunit)
       elseif (equal(word,"hard")) then
          ishard = .true.
       elseif (equal(word,"soft")) then
          ishard = .false.
       elseif (equal(word,"npts")) then
          ok = eval_next(npts,line,lp)
          if (.not.ok) then
             call ferror('struct_rdf','Incorrect NPTS',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"root")) then
          root = getword(line,lp)
       elseif (equal(word,"pair")) then
          npairs = npairs + 1
          if (npairs > size(ipairs,2)) then
             call realloc(ipairs,2,2*npairs)
          end if
          ipairs(1,npairs) = s%c%identify_spc(getword(line,lp))
          ipairs(2,npairs) = s%c%identify_spc(getword(line,lp))
          if (ipairs(1,npairs) == 0 .or. ipairs(2,npairs) == 0) then
             call ferror('struct_rdf','Unknown species in RDF/PAIR option',faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror('struct_rdf','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    call s%c%rdf(rini,rend,sigma,ishard,npts,t,ih,npairs,ipairs)

    ! write the data file
    lu = fopen_write(trim(root) // ".dat")
    write (lu,'("# ",A,A)') string("r (" // iunitname0(iunit) // ")",13,ioj_center),&
       string("RDF(r)",13,ioj_center)
    do i = 1, npts
       write (lu,'(A," ",A)') string(t(i) * dunit0(iunit),"f",15,7,ioj_center), &
          string(ih(i),"f",15,7,ioj_center)
    end do
    call fclose(lu)
    deallocate(t,ih)

    ! write the gnuplot input file
    lu = fopen_write(trim(root) // ".gnu")
    write (lu,'("set terminal postscript eps color enhanced ""Helvetica"" 25")')
    write (lu,'("set output """,A,".eps""")') trim(root)
    write (lu,*)
    if (iunit == iunit_bohr) then
       write (lu,'("set xlabel ""r (bohr)""")')
    else
       write (lu,'("set xlabel ""r (ang)""")')
    end if
    write (lu,'("set ylabel ""RDF(r)""")')
    write (lu,'("set xrange [",A,":",A,"]")') string(rini * dunit0(iunit),"f"), &
       string(rend * dunit0(iunit),"f")
    write (lu,'("set style data lines")')
    write (lu,'("set grid")')
    write (lu,'("unset key")')
    write (lu,'("plot """,A,".dat"" w lines")') string(root)
    write (lu,*)
    call fclose(lu)

    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)

  end subroutine struct_rdf

  !> Calculate the average minimum distances (AMD).
  !> Widdowson et al., "Average Minimum Distances of Periodic Point Sets-Foundational Invariants for Mapping Periodic Crystals"
  !> Match. Commun. Math. Comput. Chem., 87 (2022) 529, doi:10.46793/match.87-3.529W
  !> Reference implementation: https://github.com/dwiddo/average-minimum-distance
  module subroutine struct_amd(s,line)
    use global, only: iunitname0, dunit0, iunit
    use tools_io, only: uout, isinteger, ferror, faterr, string, ioj_left, ioj_right, warning
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line

    integer :: i, lp
    integer :: imax
    real*8, allocatable :: amd(:)
    logical :: ok

    integer, parameter :: imax_default = 100 ! default maximum nearest neighbor ordinal

    ! initialize
    imax = imax_default

    ! parse the rest of the options
    lp = 1
    ok = isinteger(imax,line,lp)
    if (.not.ok) imax = imax_default

    ! header
    write (uout,'("* AMD: average minimum distances")')
    write (uout,'("# Please cite:")')
    write (uout,'("#   Widdowson et al., Match. Commun. Math. Comput. Chem., 87 (2022) 529 (doi:10.46793/match.87-3.529W)")')
    write (uout,*)

    ! check input and allocate
    if (imax < 1) &
       call ferror('struct_amd','The size of the AMD vector must be positive',faterr)
    if (s%c%ismolecule .and. imax > s%c%ncel-1) then
       imax = s%c%ncel - 1
       call ferror('struct_amd','The AMD in a molecule can only have nn up to the number of atoms - 1',warning)
    end if

    ! calcualte the amd
    allocate(amd(imax))
    call s%c%amd(imax,amd)
    amd = amd * dunit0(iunit)

    ! print results
    write (uout,'("+ AMD vector in ",A)') iunitname0(iunit)
    write (uout,'("# nn = nearest neighbor ordinal (1 = first nn, 2 = second nn, etc.")')
    write (uout,'("# avgd = average k-nearest neighbor distance")')
    write (uout,'("# nn avgd")')
    do i = 1, imax
       write (uout,'(A," ",A)') string(i,5,ioj_left), string(amd(i),'f',10,6,ioj_right)
    end do
    write (uout,*)

  end subroutine struct_amd

  !> Compare two or more structures using the powder diffraction
  !> patterns or the radial distribution functions. Uses the
  !> similarity based on cross-correlation functions proposed in
  !>   de Gelder et al., J. Comput. Chem., 22 (2001) 273.
  module subroutine struct_compare(s,line)
    use systemmod, only: system
    use crystalmod, only: crystal, xrpd_lambda_def, xrpd_peaklist, xrpd_alpha_def,&
       xrpd_sigma_def, crosscorr_gaussian, xrpd_th2ini_def, xrpd_th2end_def, xrpd_fpol_def
    use crystalseedmod, only: struct_detect_format, struct_detect_ismol, crystalseed
    use global, only: doguess, eval_next, dunit0, iunit, iunitname0
    use tools_math, only: crosscorr_triangle, rmsd_walker, umeyama_graph_matching,&
       ullmann_graph_matching, emd, synthetic_powder
    use tools_io, only: getword, equal, faterr, ferror, uout, string, ioj_center,&
       ioj_left, string, lower, lgetword, fopen_read, fclose, isreal,&
       file_read_xy
    use types, only: realloc
    use param, only: isformat_unknown, maxzat, pi
    type(system), intent(in) :: s
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word, lword, tname, difstr, errmsg
    integer :: doguess0
    integer :: lp, i, j, k, l, n, iz, is, idx, ncount
    integer :: ns, imol, isformat, mcon, nlist
    integer :: amd_norm ! 0 = norm-inf (default), 1 = norm-1, 2 = norm-2
    type(crystal), allocatable :: c(:)
    real*8 :: tini, tend, nor, h, xend, sigma, epsreduce, diffmin, diff2, xnorm2, lambda
    character*15 :: auxdif(5)
    type(xrpd_peaklist), allocatable :: xp(:)
    real*8, allocatable :: t(:), ih(:), iha(:,:), th2a(:,:)
    real*8, allocatable :: dref(:,:), ddg(:,:), ddh(:,:)
    integer, allocatable :: zcount1(:), zcount2(:), isperm(:), list(:,:), na(:), irepeat(:)
    real*8, allocatable :: diff(:,:), xnorm(:), x1(:,:), x2(:,:), ip(:), th2p(:)
    integer, allocatable :: iz1(:), ncon1(:), idcon1(:,:), iz2(:), ncon2(:), idcon2(:,:)
    integer, allocatable :: singleatom(:)
    logical :: ok, noh
    logical :: ismol, laux, lzc
    character*1024, allocatable :: fname(:)
    type(crystalseed) :: seed
    integer, allocatable :: fname_type(:)

    integer, parameter :: msg_counter = 50

    integer, parameter :: npts = 10001
    real*8, parameter :: rend0 = 25d0
    integer, parameter :: amd_imax = 100 ! length of the AMD vector
    real*8, parameter :: emd_discardp = 0.1d0 ! discard signals below this intentsity (highest I is 100)

    integer, parameter :: fname_current = 0
    integer, parameter :: fname_structure = 1
    integer, parameter :: fname_peaks = 2

    integer :: imethod
    integer, parameter :: imethod_default = 0
    integer, parameter :: imethod_gpdwf = 1
    integer, parameter :: imethod_powder = 2
    integer, parameter :: imethod_rdf = 3
    integer, parameter :: imethod_sorted = 4
    integer, parameter :: imethod_umeyama = 5
    integer, parameter :: imethod_ullmann = 6
    integer, parameter :: imethod_amd = 7
    integer, parameter :: imethod_emd = 8

    ! initialize
    imethod = imethod_default
    lp = 1
    ns = 0
    xend = -1d0
    allocate(fname(10),fname_type(10))
    sigma = xrpd_sigma_def
    epsreduce = -1d0
    imol = -1
    amd_norm = 0
    noh = .false.
    lambda = xrpd_lambda_def

    ! read the input options
    do while(.true.)
       word = getword(line,lp)
       lword = lower(word)
       if (equal(lword,'xend')) then
          ok = eval_next(xend,line,lp)
          if (.not.ok) then
             call ferror('struct_compare','incorrect TH2END',faterr,syntax=.true.)
             return
          end if
       elseif (equal(lword,'sigma')) then
          ok = eval_next(sigma,line,lp)
          if (.not.ok) then
             call ferror('struct_compare','incorrect SIGMA',faterr,syntax=.true.)
             return
          end if
       elseif (equal(lword,'lambda')) then
          ok = eval_next(lambda,line,lp)
          if (.not.ok) then
             call ferror('struct_compare','incorrect LAMBDA',faterr,syntax=.true.)
             return
          end if
       elseif (equal(lword,'reduce')) then
          ok = eval_next(epsreduce,line,lp)
          if (.not.ok) then
             call ferror('struct_compare','incorrect REDUCE',faterr,syntax=.true.)
             return
          end if
       elseif (equal(lword,'gpwdf')) then
          imethod = imethod_gpdwf
       elseif (equal(lword,'powder')) then
          imethod = imethod_powder
       elseif (equal(lword,'rdf')) then
          imethod = imethod_rdf
       elseif (equal(lword,'amd')) then
          imethod = imethod_amd
       elseif (equal(lword,'ullmann')) then
          imethod = imethod_ullmann
       elseif (equal(lword,'umeyama')) then
          imethod = imethod_umeyama
       elseif (equal(lword,'sorted')) then
          imethod = imethod_sorted
       elseif (equal(lword,'emd')) then
          imethod = imethod_emd
       elseif (equal(lword,'molecule')) then
          imol = 1
       elseif (equal(lword,'crystal')) then
          imol = 0
       elseif (equal(lword,'noh')) then
          noh = .true.
       elseif (equal(lword,'norm')) then
          word = lgetword(line,lp)
          if (equal(word,'1')) then
             amd_norm = 1
          elseif (equal(word,'2')) then
             amd_norm = 2
          elseif (equal(word,'inf')) then
             amd_norm = 0
          else
             call ferror('struct_compare','incorrect NORM',faterr,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          ns = ns + 1
          if (ns > size(fname)) then
             call realloc(fname,2*ns)
             call realloc(fname_type,2*ns)
          end if
          fname(ns) = word
          if (equal(fname(ns),".")) then
             fname_type(ns) = fname_current
          else
             idx = index(word,'.',.true.)
             if (equal(word(idx+1:),"peaks")) then
                fname_type(ns) = fname_peaks
             else
                fname_type(ns) = fname_structure
             end if
          end if
       else
          exit
       end if
    end do
    if (ns < 2) &
       call ferror('struct_compare','At least 2 structures are needed for the comparison',faterr)

    ! no symmetry for new structures
    doguess0 = doguess
    doguess = 0

    ! determine whether to use crystal or molecule comparison
    if (imol == 0) then
       ismol = .false.
    elseif (imol == 1) then
       ismol = .true.
    else
       ismol = .true.
       do i = 1, ns
          if (fname_type(i) == fname_current) then
             if (.not.s%c%isinit) &
                call ferror('struct_compare','Current structure (".") is not initialized.',faterr)
             ismol = ismol .and. s%c%ismolecule
          elseif (fname_type(i) == fname_peaks) then
             ismol = .false.
             inquire(file=fname(i),exist=laux)
             if (.not.laux) &
                call ferror("struct_compare","file not found: " // string(fname(i)),faterr)
          elseif (fname_type(i) == fname_structure) then
             call struct_detect_format(fname(i),isformat)
             call struct_detect_ismol(fname(i),isformat,laux)
             ismol = ismol .and. laux
             if (isformat == isformat_unknown) &
                call ferror("struct_compare","unknown file format: " // string(fname(i)),faterr)
             inquire(file=fname(i),exist=laux)
             if (.not.laux) &
                call ferror("struct_compare","file not found: " // string(fname(i)),faterr)
          end if
       end do
       if (ismol) then
          imol = 1
       else
          imol = 0
       end if
    end if

    ! Read the structures and write header
    write (uout,'("* COMPARE: compare structures")')
    if (ismol) then
       tname = "Molecule"
    else
       tname = "Crystal"
    end if
    allocate(c(ns))
    do i = 1, ns
       if (fname_type(i) == fname_current) then
          write (uout,'("  ",A," ",A,": <current>")') string(tname), string(i,2)
          c(i) = s%c
       elseif (fname_type(i) == fname_peaks) then
          write (uout,'("  ",A," ",A,": ",A)') string(tname), string(i,2), string(fname(i))
          if (ismol) &
             call ferror("struct_compare","reading peaks files only allowed for crystals",faterr)
          if (imethod /= imethod_emd .and. imethod /= imethod_powder .and. imethod /= imethod_default) &
             call ferror("struct_compare","reading peaks files only allowed for EMD and POWDER",faterr)
       else
          write (uout,'("  ",A," ",A,": ",A)') string(tname), string(i,2), string(fname(i))
          call struct_crystal_input(fname(i),imol,.false.,.false.,cr0=c(i))
          if (.not.c(i)%isinit) &
             call ferror("struct_compare","could not load crystal structure" // string(fname(i)),faterr)
       end if
    end do

    ! decide on the method, if it is not already set to a valid value
    if (ismol) then
       ! get the Z count and the sequence for the first molecule
       allocate(zcount1(maxzat),zcount2(maxzat))
       zcount1 = 0
       do i = 1, c(1)%ncel
          iz = c(1)%spc(c(1)%atcel(i)%is)%z
          zcount1(iz) = zcount1(iz) + 1
       end do

       ! get the Z count for the other molecules
       lzc = .true.
       do is = 2, ns
          zcount2 = 0
          do i = 1, c(is)%ncel
             iz = c(is)%spc(c(is)%atcel(i)%is)%z
             zcount2(iz) = zcount2(iz) + 1
          end do
          lzc = lzc .and. all(zcount2 == zcount1)
          if (.not.lzc) exit
       end do
       deallocate(zcount1,zcount2)

       ! rdf if the atom counts are not the same, otherwise ullmann is default
       if (.not.lzc) then
          imethod = imethod_rdf
       elseif (imethod == imethod_default .or. imethod == imethod_powder) then
          imethod = imethod_ullmann
       end if
    else
       if (imethod /= imethod_rdf .and. imethod /= imethod_amd .and.&
           imethod /= imethod_emd .and. imethod /= imethod_powder)&
           imethod = imethod_gpdwf
    end if

    ! strip the hydrogens
    if (noh) then
       write (uout,'("# Removing hydrogens from the structures.")')
       do i = 1, ns
          if (fname_type(i) == fname_peaks) cycle
          call c(i)%makeseed(seed,.true.)
          call seed%strip_hydrogens()
          call c(i)%struct_new(seed,.true.)
       end do
    end if

    ! rest of the header and default variables
    if (imethod == imethod_gpdwf) then
       write (uout,'("# Using cross-correlated Gaussian powder diffraction patterns (GPWDF).")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   Otero-de-la-Roza, J. Appl. Cryst. 57 (2024) 1401-1414 (doi:10.1107/S1600576724007489)")')
       write (uout,'("# Two structures are exactly equal if DIFF = 0.")')
       if (xend < 0d0) xend = xrpd_th2end_def
       difstr = "DIFF"
    elseif (imethod == imethod_powder) then
       write (uout,'("# Using cross-correlated POWDER diffraction patterns.")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   de Gelder et al., J. Comput. Chem., 22 (2001) 273")')
       write (uout,'("#       (10.1002/1096-987X(200102)22:3%3C273::AID-JCC1001%3E3.0.CO;2-0)")')
       write (uout,'("# Two structures are exactly equal if DIFF = 0.")')
       if (xend < 0d0) xend = xrpd_th2end_def
       difstr = "DIFF"
    elseif (imethod == imethod_emd) then
       write (uout,'("# Using discrete powder diffraction patterns and the earth mover''s distance (EMD).")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   Rubner et al., Int. J. Comput. Vis. 40.2 (2000) 99-121 (10.1023/A:1026543900054)")')
       write (uout,'("# Two structures are exactly equal if DIFF = 0.")')
       if (xend < 0d0) xend = xrpd_th2end_def
       difstr = "DIFF"
    elseif (imethod == imethod_rdf) then
       write (uout,'("# Using cross-correlated radial distribution functions (RDF).")')
       write (uout,'("# Two structures are exactly equal if DIFF = 0.")')
       if (xend < 0d0) xend = rend0
       difstr = "DIFF"
    elseif (imethod == imethod_amd) then
       write (uout,'("# Using average minimum distances (AMD) in bohr.")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   Widdowson et al., Match. Commun. Math. Comput. Chem., 87 (2022) 529")')
       write (uout,'("#       (doi:10.46793/match.87-3.529W)")')
       write (uout,'("# Two structures are exactly equal if DIFF = 0.")')
       if (xend < 0d0) xend = rend0
       difstr = "DIFF"
    else
       if (imethod == imethod_sorted) then
          write (uout,'("# Assuming the atom sequence is the same in all molecules.")')
       elseif (imethod == imethod_umeyama) then
          write (uout,'("# Using Umeyama''s graph matching algorithm.")')
       elseif (imethod == imethod_ullmann) then
          write (uout,'("# Using Ullmann''s graph matching algorithm.")')
       end if
       write (uout,'("# RMS of the atomic positions in ",A)') iunitname0(iunit)
       difstr = "RMS"
    end if
    if (epsreduce > 0d0) then
       write (uout,'("# REDUCE: Only repeated/unique structures will be listed, at DIFF level: ",A)') &
          string(epsreduce,'e',10,4)
       allocate(irepeat(ns))
       irepeat = 0
    end if

    ! allocate space for difference/rms values
    allocate(diff(ns,ns))
    diff = 1d0

    if (imethod == imethod_gpdwf) then
       ! crystals: GPWDF

       ! calculate all the peak lists
       allocate(xp(ns))
       ncount = 0
       do i = 1, ns
          errmsg = ""
          ncount = ncount + 1
          if (mod(ncount-1,msg_counter) == 0) &
          write (uout,'("  ... calculating pattern ",A," of ",A,".")') string(ncount), string(ns)
          if (fname_type(i) == fname_current) then
             call xp(i)%from_crystal(s%c,xrpd_th2ini_def,xend,lambda,xrpd_fpol_def,.false.,.false.,errmsg)
          elseif (fname_type(i) == fname_peaks) then
             call xp(i)%from_peaks_file(fname(i),errmsg)
          elseif (fname_type(i) == fname_structure) then
             call xp(i)%from_crystal(c(i),xrpd_th2ini_def,xend,lambda,xrpd_fpol_def,.false.,.false.,errmsg)
          end if
          if (len(errmsg) > 0) then
             call ferror("struct_compare","error calculating pattern" // string(i) // " from file " //&
                string(fname(i)) // ": " // errmsg,faterr)
          end if
       end do
       write (uout,'("  ... finished calculating patterns")')

       ! self-correlation
       allocate(xnorm(ns))
       do i = 1, ns
          call crosscorr_gaussian(xp(i),xp(i),xrpd_alpha_def,sigma,&
             xnorm(i),errmsg,.false.)
          if (len(errmsg) > 0) &
             call ferror("struct_compare","error calculating Gaussian crosscorrelation",faterr)
       end do
       xnorm = sqrt(abs(xnorm))

       ! calculate the overlap between diffraction patterns
       diff = 0d0
       ncount = 0
       do i = 1, ns
          ncount = ncount + 1
          if (mod(i-1,msg_counter) == 0) &
             write (uout,'("  ... comparing pattern ",A," of ",A,".")') string(i), string(ns)
          if (epsreduce > 0d0) then
             if (irepeat(i) > 0) cycle
          end if
          !$omp parallel do firstprivate(errmsg) private(xnorm2)
          do j = i+1, ns
             if (epsreduce > 0d0) then
                if (irepeat(j) > 0) cycle
             end if
             call crosscorr_gaussian(xp(i),xp(j),xrpd_alpha_def,sigma,xnorm2,errmsg,.false.)
             !$omp critical (errmsg_)
             if (len(errmsg) > 0) &
                call ferror("struct_compare","error calculating Gaussian crosscorrelation",faterr)
             !$omp end critical (errmsg_)
             diff(i,j) = max(1d0 - xnorm2 / xnorm(i) / xnorm(j),0d0)
             diff(j,i) = diff(i,j)

             if (epsreduce > 0d0) then
                if (diff(i,j) < epsreduce) irepeat(j) = i
             end if
          end do
          !$omp end parallel do
       end do
       deallocate(xnorm,xp)
       write (uout,'("  ... finished comparing patterns")')

    elseif (imethod == imethod_powder .or. imethod == imethod_rdf) then
       ! crystals: POWDER and RDF
       allocate(iha(npts,ns),singleatom(ns),t(npts))
       singleatom = -1

       ncount = 0
       !$omp parallel do private(tini,tend,nor,n) firstprivate(t,ih,th2p,ip)
       do i = 1, ns
          !$omp critical (ncount_)
          ncount = ncount + 1
          if (mod(ncount-1,msg_counter) == 0) &
             write (uout,'("  ... calculating pattern ",A," of ",A,".")') string(ncount), string(ns)
          !$omp end critical (ncount_)

          if (fname_type(i) == fname_peaks) then
             !$omp critical (fileread)
             call file_read_xy(fname(i),n,th2p,ip)
             call synthetic_powder(xrpd_th2ini_def,xend,npts,sigma,th2p,ip,t,ih)
             !$omp end critical (fileread)
             iha(:,i) = ih

             ! normalize the integral of abs(ih)
             tini = ih(1)**2
             tend = ih(npts)**2
             nor = (2d0 * sum(ih(2:npts-1)**2) + tini + tend) * (xend - xrpd_th2ini_def) / 2d0 / real(npts-1,8)
             iha(:,i) = ih / sqrt(nor)
          else
             ! calculate the powder diffraction pattern
             if (imethod == imethod_powder) then
                call c(i)%powder(0,xrpd_th2ini_def,xend,lambda,xrpd_fpol_def,npts=npts,sigma=sigma,ishard=.false.,&
                   t=t,ih=ih)

                ! normalize the integral of abs(ih)
                tini = ih(1)**2
                tend = ih(npts)**2
                nor = (2d0 * sum(ih(2:npts-1)**2) + tini + tend) * (xend - xrpd_th2ini_def) / 2d0 / real(npts-1,8)
                iha(:,i) = ih / sqrt(nor)
             else
                if (c(i)%ncel == 1) then
                   singleatom(i) = c(i)%spc(c(i)%atcel(1)%is)%z
                else
                   call c(i)%rdf(0d0,xend,sigma,.false.,npts,t,ih)
                   iha(:,i) = ih
                end if
             end if
          end if
       end do
       !$omp end parallel do
       if (allocated(t)) deallocate(t)
       if (allocated(ih)) deallocate(ih)
       write (uout,'("  ... finished calculating patterns")')

       ! self-correlation
       allocate(xnorm(ns))
       h =  (xend-xrpd_th2ini_def) / real(npts-1,8)
       do i = 1, ns
          if (singleatom(i) < 0) then
             xnorm(i) = crosscorr_triangle(h,iha(:,i),iha(:,i),1d0)
          end if
       end do
       xnorm = sqrt(abs(xnorm))

       ! calculate the overlap between diffraction patterns
       diff = 0d0
       ncount = 0
       do i = 1, ns
          ncount = ncount + 1
          if (mod(i-1,msg_counter) == 0) &
             write (uout,'("  ... comparing pattern ",A," of ",A,".")') string(i), string(ns)
          if (epsreduce > 0d0) then
             if (irepeat(i) > 0) cycle
          end if
          !$omp parallel do
          do j = i+1, ns
             if (epsreduce > 0d0) then
                if (irepeat(j) > 0) cycle
             end if
             if (singleatom(i) > 0 .and. singleatom(j) > 0) then
                if (singleatom(i) == singleatom(j)) then
                   diff(i,j) = 0d0
                else
                   diff(i,j) = 1d0
                end if
             else if (singleatom(i) > 0 .or. singleatom(j) > 0) then
                diff(i,j) = 1d0
             else
                diff(i,j) = max(1d0 - crosscorr_triangle(h,iha(:,i),iha(:,j),1d0) / xnorm(i) / xnorm(j),0d0)
             end if
             diff(j,i) = diff(i,j)

             if (epsreduce > 0d0) then
                if (diff(i,j) < epsreduce) irepeat(j) = i
             end if
          end do
          !$omp end parallel do
       end do
       deallocate(xnorm)
       write (uout,'("  ... finished comparing patterns")')
    elseif (imethod == imethod_emd) then
       ! crystals: EMD
       allocate(na(ns),th2a(100,ns),iha(100,ns))
       do i = 1, ns
          if (fname_type(i) == fname_peaks) then
             call file_read_xy(fname(i),n,th2p,ip)
          else
             n = 0
             ! calculate the powder diffraction pattern
             call c(i)%powder(0,xrpd_th2ini_def,xend,lambda,xrpd_fpol_def,npts=npts,sigma=sigma,ishard=.false.,&
                th2p=th2p,ip=ip,discardp=emd_discardp)

             ! normalize
             n = size(ip,1)
             ip(1:n) = ip(1:n) / sum(ip(1:n))
             th2p(1:n) = th2p(1:n) * 180d0 / pi
          end if

          ! write it down
          if (n > 0) then
             if (n > size(na,1)) then
                call realloc(na,2*n)
                call realloc(th2a,2*n,ns)
                call realloc(iha,2*n,ns)
             end if
             na(i) = n
             th2a(1:n,i) = th2p(1:n)
             iha(1:n,i) = ip(1:n)
          end if
       end do

       ! calculate the distance between diffraction patterns
       diff = 0d0
       do i = 1, ns
          do j = i+1, ns
             diff(i,j) = emd(1,na(i),th2a(1:na(i),i),iha(1:na(i),i),na(j),th2a(1:na(j),j),iha(1:na(j),j),1)
             diff(j,i) = diff(i,j)
             if (epsreduce > 0d0) then
                if (diff(i,j) < epsreduce) irepeat(j) = i
             end if
          end do
       end do

    elseif (imethod == imethod_amd) then
       ! crystals: AMD
       allocate(iha(amd_imax,ns))
       do i = 1, ns
          ! calculate the AMD
          call c(i)%amd(amd_imax,iha(:,i))
       end do

       ! compare
       diff = 0d0
       allocate(ip(amd_imax))
       do i = 1, ns
          do j = i+1, ns
             ip = abs(iha(:,i)-iha(:,j))
             if (amd_norm == 0) then
                diff(i,j) = maxval(ip)
             elseif (amd_norm == 1) then
                diff(i,j) = sum(ip)
             else
                diff(i,j) = norm2(ip)
             end if
             diff(j,i) = diff(i,j)
             if (epsreduce > 0d0) then
                if (diff(i,j) < epsreduce) irepeat(j) = i
             end if
          end do
       end do
       deallocate(ip)
    else
       ! molecules, ullmann, sorted, and umeyama. We have already made sure
       ! that the atoms counts are the same, so now we need the permutation
       ! prepare...
       n = c(1)%ncel
       allocate(isperm(n),x1(3,n),x2(3,n))
       if (imethod == imethod_umeyama) then
          allocate(dref(n,n),ddg(n,n),ddh(n,n))
       elseif (imethod == imethod_ullmann) then
          mcon = 0
          do i = 1, ns
             do j = 1, n
                mcon = max(mcon,c(i)%nstar(j)%ncon)
             end do
          end do
          allocate(iz1(n),ncon1(n),idcon1(mcon,n),iz2(n),ncon2(n),idcon2(mcon,n))
       end if

       ! run over structure pairs
       diff = 0d0
       do i = 1, ns
          ! save the atomic coordinates to x1
          do k = 1, n
             x1(:,k) = c(i)%atcel(k)%r + c(i)%molx0
          end do
          if (imethod == imethod_umeyama) then
             call c(i)%distmatrix(dref,conn=.true.)
          elseif (imethod == imethod_ullmann) then
             idcon1 = 0
             do k = 1, n
                iz1(k) = c(i)%spc(c(i)%atcel(k)%is)%z
                ncon1(k) = c(i)%nstar(k)%ncon
                idcon1(1:ncon1(k),k) = c(i)%nstar(k)%idcon(1:ncon1(k))
             end do
          end if

          do j = i+1, ns

             if (imethod == imethod_umeyama) then
                ddg = dref
                call c(j)%distmatrix(ddh,conn=.true.)
                call umeyama_graph_matching(n,ddg,ddh,isperm)
                do k = 1, n
                   x2(:,k) = c(j)%atcel(isperm(k))%r + c(j)%molx0
                end do
                diffmin = rmsd_walker(x1,x2)
             elseif (imethod == imethod_ullmann) then
                idcon2 = 0
                do k = 1, n
                   iz2(k) = c(j)%spc(c(j)%atcel(k)%is)%z
                   ncon2(k) = c(j)%nstar(k)%ncon
                   idcon2(1:ncon2(k),k) = c(j)%nstar(k)%idcon(1:ncon2(k))
                end do
                call ullmann_graph_matching(iz1,ncon1,idcon1,iz2,ncon2,idcon2,nlist,list)

                diffmin = huge(1d0)
                do k = 1, nlist
                   isperm = list(:,k)
                   do l = 1, n
                      x2(:,l) = c(j)%atcel(isperm(l))%r + c(j)%molx0
                   end do
                   diff2 = rmsd_walker(x1,x2)
                   if (diff2 < diffmin) diffmin = diff2
                end do
             else
                do k = 1, n
                   x2(:,k) = c(j)%atcel(k)%r + c(j)%molx0
                end do
                diffmin = rmsd_walker(x1,x2)
             end if

             ! calculate the RMS
             diff(i,j) = diffmin
             diff(j,i) = diff(i,j)
             if (epsreduce > 0d0) then
                if (diff(i,j) < epsreduce) irepeat(j) = i
             end if
          end do
       end do
       diff = diff * dunit0(iunit)
       deallocate(isperm,x1,x2)
       if (allocated(dref)) deallocate(dref)
       if (allocated(ddg)) deallocate(ddg)
       if (allocated(ddh)) deallocate(ddh)
       if (allocated(iz1)) deallocate(iz1)
       if (allocated(ncon1)) deallocate(ncon1)
       if (allocated(idcon1)) deallocate(idcon1)
       if (allocated(iz2)) deallocate(iz2)
       if (allocated(ncon2)) deallocate(ncon2)
       if (allocated(idcon2)) deallocate(idcon2)
    endif

    ! write output
    if (ns == 2) then
       write (uout,'("+ ",A," = ",A)') string(difstr), string(diff(1,2),'e',12,6)
    else
       do i = 0, (ns-1)/5
          write (uout,'(99(A," "))') string(tname,15,ioj_center), &
             (string(c(5*i+j)%file,15,ioj_center),j=1,min(5,ns-i*5))
          write (uout,'(99(A," "))') string(difstr,15,ioj_center), &
             (string(5*i+j,15,ioj_center),j=1,min(5,ns-i*5))
          do j = 1, ns
             do k = 1, min(5,ns-i*5)
                auxdif(k) = string(diff(j,5*i+k),'f',15,7,3)
                if (epsreduce > 0d0) then
                   if (diff(j,5*i+k) < epsreduce) auxdif(k) = string(" not-calculated")
                end if
             end do
             write (uout,'("  ",99(A," "))') string(c(j)%file,15,ioj_left), &
                (auxdif(k),k=1,min(5,ns-i*5))
          end do
          write (uout,*)
       end do
    endif
    write (uout,*)

    ! reduce
    if (epsreduce > 0d0) then
       write (uout,'("+ List of unique structures (",A,"): ")') string(count(irepeat == 0))
       do i = 1, ns
          if (irepeat(i) == 0) write (uout,'(A,": ",A," with multiplicity ",A)') &
             string(i), trim(c(i)%file), string(count(irepeat == i) + 1)
       end do
       write (uout,*)

       write (uout,'("+ List of repeated structures (",A,"): ")') string(ns-count(irepeat == 0))
       do i = 1, ns
          if (irepeat(i) > 0) write (uout,'(A,": ",A," same as ",A,": ",A)') &
             string(i), trim(c(i)%file), string(irepeat(i)), trim(c(irepeat(i))%file)
       end do
       write (uout,*)

       deallocate(irepeat)
    end if

    ! clean up
    deallocate(diff)
    doguess = doguess0

  end subroutine struct_compare

  !> Compare structures, allowing deformation of one crystal into the
  !> other such as may be cause by temperature or pressure effects
  !> without a phase change. Variable-cell version of the POWDIFF
  !> comparison method in struct_compare.
  module subroutine struct_comparevc(s,line)
    use tools_io, only: lgetword, equal
    type(system), intent(in) :: s
    character*(*), intent(in) :: line

    integer :: lp
    character(len=:), allocatable :: word

    ! header and initalization
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,'vcpwdf')) then
       call struct_comparevc_vcpwdf(s,line(lp:))
    elseif (equal(word,'vcgpwdf')) then
       call struct_comparevc_vcgpwdf(s,line(lp:))
    else
#ifdef HAVE_NLOPT
       call struct_comparevc_vcgpwdf(s,line)
#else
       call struct_comparevc_vcpwdf(s,line)
#endif
    end if

  end subroutine struct_comparevc

  !> COMPAREVC: vcpwdf version. Transform both cells to the primitive, then
  !> calculate the deformation that brings one into agreement with the other.
  !> Returns the minimal powder diffraction pattern (de Gelder) overlap. See:
  !>   Mayo et al., CrystEngComm 24 (2022) 8326-8338.
  subroutine struct_comparevc_vcpwdf(s,line)
    use spglib, only: spg_delaunay_reduce, spg_standardize_cell
    use global, only: fileroot
    use crystalmod, only: crystal, vcpwdf_compare
    use crystalseedmod, only: crystalseed
    use tools_math, only: matinv, m_c2x_from_cellpar, det3, crosscorr_triangle, &
       m_x2c_from_cellpar
    use tools_io, only: getword, faterr, ferror, uout, string, isreal, equal, lgetword
    type(system), intent(in) :: s
    character*(*), intent(in) :: line

    type(crystalseed) :: seed
    integer :: lp
    character(len=:), allocatable :: file1, file2, errmsg, word
    type(crystal) :: c1, c2, c1out, c2out
    real*8 :: diff
    logical :: ok, dowrite, noh
    real*8 :: powdiff_thr, max_elong, max_ang, max_vol

    real*8, parameter :: max_elong_def = 0.3d0 ! at most 30% elongation of cell lengths
    real*8, parameter :: max_ang_def = 20d0    ! at most 20 degrees change in angle
    real*8, parameter :: max_vol_def = 0.5d0 ! at most 50% change in volume

    ! header and initalization
    write (uout,'("* COMPAREVC: compare crystal structures allowing for deformed cells")')
    write (uout,'("+ VCPWDF method: exploring lattice deformations until a match is found")')
    write (uout,'("# Please cite:")')
    write (uout,'("#   Mayo et al., CrystEngComm 24 (2022) 8326-8338. (doi:10.1039/D2CE01080A)")')
    write (uout,*)

    ! read the input files
    lp = 1
    file1 = getword(line,lp)
    if (len_trim(file1) == 0) &
       call ferror('struct_comparevc_vcpwdf','Missing first structure file',faterr)
    file2 = getword(line,lp)
    if (len_trim(file2) == 0) &
       call ferror('struct_comparevc_vcpwdf','Missing second structure file',faterr)

    max_elong = max_elong_def
    max_ang = max_ang_def
    max_vol = max_vol_def
    powdiff_thr = -1d0
    dowrite = .false.
    noh = .false.
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'thr')) then
          ok = isreal(powdiff_thr,line,lp)
          if (.not.ok) call ferror('struct_comparevc_vcpwdf','Wrong THR',faterr)
       elseif (equal(word,'maxelong')) then
          ok = isreal(max_elong,line,lp)
          if (.not.ok) call ferror('struct_comparevc_vcpwdf','Wrong MAXELONG',faterr)
       elseif (equal(word,'maxang')) then
          ok = isreal(max_ang,line,lp)
          if (.not.ok) call ferror('struct_comparevc_vcpwdf','Wrong MAXANG',faterr)
       elseif (equal(word,'maxvol')) then
          ok = isreal(max_vol,line,lp)
          if (.not.ok) call ferror('struct_comparevc_vcpwdf','Wrong MAXVOL',faterr)
       elseif (equal(word,'write')) then
          dowrite = .true.
       elseif (equal(word,'noh')) then
          noh = .true.
       elseif (len_trim(word) > 0) then
          if (.not.ok) call ferror('struct_comparevc_vcpwdf','Unknown keyword',faterr)
       else
          exit
       end if
    end do

    ! read the structures, force symmetry recalculation
    if (equal(file1,".")) then
       write (uout,'("+ Structure from currently loaded file: ",A)') trim(s%c%file)
       c1 = s%c
       file1 = s%c%file
    else
       write (uout,'("+ Reading the structure from: ",A)') trim(file1)
       call seed%read_any_file(file1,-1,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_comparevc_vcpwdf','error reading geometry file: ' // file1,faterr)
       if (noh) call seed%strip_hydrogens()
       call c1%struct_new(seed,.true.)
       call c1%calcsym(.false.,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_comparevc_vcpwdf','error recalculating symmetry: ' // file1,faterr)
    end if

    ! read the second structure, force symmetry recalculation
    if (equal(file2,".")) then
       write (uout,'("+ Structure from currently loaded file: ",A)') trim(s%c%file)
       c2 = s%c
       file2 = s%c%file
    else
       write (uout,'("+ Reading the structure from: ",A)') trim(file2)
       call seed%read_any_file(file2,-1,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_comparevc_vcpwdf','error reading geometry file: ' // file2,faterr)
       if (noh) call seed%strip_hydrogens()
       call c2%struct_new(seed,.true.)
       call c2%calcsym(.false.,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_comparevc_vcpwdf','error recalculating symmetry: ' // file1,faterr)
    end if

    ! run the comparison
    if (dowrite) then
       call vcpwdf_compare(c1,c2,diff,errmsg,max_elong,max_ang,max_vol,powdiff_thr,c1out,c2out,verbose=.true.)
    else
       call vcpwdf_compare(c1,c2,diff,errmsg,max_elong,max_ang,max_vol,powdiff_thr,verbose=.true.)
    end if
    if (len_trim(errmsg) > 0) &
       call ferror('struct_comparevc_vcpwdf',errmsg,faterr)

    ! write the final structure  to output
    if (dowrite) then
       ! write both to a res file
       file1 = fileroot // "_structure_1.res"
       write (uout,'("+ Structure 1 written to file: ",A)') trim(file1)
       call c1out%write_res(file1,-1)

       file2 = fileroot // "_structure_2.res"
       write (uout,'("+ Structure 2 written to file: ",A)') trim(file2)
       call c2out%write_res(file2,-1)
    end if

  end subroutine struct_comparevc_vcpwdf

  !> COMPAREVC: vcgpwdf version. Do a global search over lattice
  !> deformations of the first structure to have a list of diffraction
  !> angles and intensities that best matches a second structure or an
  !> experimental pattern. See:
  !> Otero-de-la-Roza, J. Appl. Cryst. 57 (2024) 1401-1414 (doi:10.1107/S1600576724007489)
  subroutine struct_comparevc_vcgpwdf(s,line)
    use tools_io, only: uout, ferror, faterr
#ifdef HAVE_NLOPT
    use global, only: fileroot
    use crystalseedmod, only: crystalseed
    use crystalmod, only: crystal, xrpd_peaklist, xrpd_fpol_def,&
       gaussian_compare, xrpd_alpha_def, xrpd_lambda_def,&
       xrpd_th2ini_def, xrpd_th2end_def, xrpd_sigma_def
    use struct_drivers, only: struct_crystal_input
    use tools, only: qcksort
    use tools_io, only: getword, tictac, lgetword, equal, string,&
       fopen_read, fclose, getline, isreal, isinteger, fopen_write, ioj_center
    use types, only: realloc
#endif
    type(system), intent(in) :: s
    character*(*), intent(in) :: line

#ifdef HAVE_NLOPT
    integer :: lp
    integer :: i
    real*8 :: diff
    real*8, allocatable :: t(:), ih(:)
    character(len=:), allocatable :: word, file1, errmsg
    type(crystal) :: c1, c2
    integer :: imode, lu, maxfeval
    real*8 :: th2ini, th2end, alpha, lambda, besteps, max_elong, max_ang
    logical :: ok, readc2, dowrite
    type(crystalseed) :: seed
    type(xrpd_peaklist) :: p2

    ! global block
    integer :: neval
    real*8 :: lastval
    real*8 :: bestval
    integer :: nbesteval

    ! parameters
    integer, parameter :: imode_sp = 0
    integer, parameter :: imode_local = 1
    integer, parameter :: imode_global = 2
    integer, parameter :: final_npts = 10001

    include 'nlopt.f'
#endif

    ! header and initalization
    write (uout,'("* COMPAREVC: compare crystal structures allowing for deformed cells")')
    write (uout,'("+ VCGPWDF method: global lattice search for maximal powder diffraction overlap")')
    write (uout,'("# Please cite:")')
    write (uout,'("#   Otero-de-la-Roza, J. Appl. Cryst. 57 (2024) 1401-1414 (doi:10.1107/S1600576724007489)")')
    write (uout,*)

#ifndef HAVE_NLOPT
    call ferror("struct_comparevc_vcgpwdf","COMPARE VCGPWDF can only be used if the NLOPT library is available",faterr)
#else
    ! initialize
    lp = 1
    neval = 0
    lastval = -1d0
    bestval = 1.1d0
    nbesteval = 0
    lambda = xrpd_lambda_def

    ! read crystal structures
    lp = 1
    file1 = getword(line,lp)
    write (uout,'("  Crystal 1: ",A)') string(file1)
    if (file1 == ".") then
       c1 = s%c
    else
       call struct_crystal_input(file1,0,.false.,.false.,cr0=c1)
    end if
    word = getword(line,lp)
    write (uout,'("  Crystal 2: ",A)') string(word)
    if (word == ".") then
       th2ini = xrpd_th2ini_def
       th2end = xrpd_th2end_def
       c2 = s%c
       readc2 = .true.
    elseif (len(word) - index(word,".peaks") == 6) then
       call p2%from_peaks_file(word,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror("struct_comparevc_vcgpwdf",errmsg,faterr)
       th2ini = p2%th2(1) - 1d-2
       th2end = p2%th2(p2%npeak) + 1d-2
       readc2 = .false.
    else
       ! read crystal structure 2
       th2ini = xrpd_th2ini_def
       th2end = xrpd_th2end_def
       call struct_crystal_input(word,0,.false.,.false.,cr0=c2)
       readc2 = .true.
    end if
    write (uout,'("  Initial and final 2*theta: ",A," ",A)') &
       string(th2ini,'f',decimal=4), string(th2end,'f',decimal=4)

    ! read additional options
    call set_quick_params()
    alpha = xrpd_alpha_def
    dowrite = .false.
    imode = imode_global
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'sp')) then
          imode = imode_sp
       elseif (equal(word,'local')) then
          imode = imode_local
       elseif (equal(word,'global')) then
          imode = imode_global
       elseif (equal(word,'quick')) then
          call set_quick_params()
       elseif (equal(word,'safe')) then
          call set_safe_params()
       elseif (equal(word,'alpha')) then
          ok = isreal(alpha,line,lp)
          if (.not.ok) then
             call ferror("struct_comparevc_vcgpwdf","invalid alpha in gaucomp",faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'maxfeval')) then
          ok = isinteger(maxfeval,line,lp)
          if (.not.ok) then
             call ferror("struct_comparevc_vcgpwdf","invalid maxfeval in gaucomp",faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'besteps')) then
          ok = isreal(besteps,line,lp)
          if (.not.ok) then
             call ferror("struct_comparevc_vcgpwdf","invalid besteps in gaucomp",faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'lambda')) then
          ok = isreal(lambda,line,lp)
          if (.not.ok) then
             call ferror("struct_comparevc_vcgpwdf","invalid lambda in gaucomp",faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'maxelong')) then
          ok = isreal(max_elong,line,lp)
          if (.not.ok) then
             call ferror("struct_comparevc_vcgpwdf","invalid max_elong in gaucomp",faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'maxang')) then
          ok = isreal(max_ang,line,lp)
          if (.not.ok) then
             call ferror("struct_comparevc_vcgpwdf","invalid max_ang in gaucomp",faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'write')) then
          dowrite = .true.
       else
          exit
       end if
    end do

    ! pre-calculation
    if (readc2) then
       call p2%from_crystal(c2,th2ini,th2end,lambda,xrpd_fpol_def,.false.,.false.,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror("struct_comparevc_vcgpwdf",errmsg,faterr)
    end if

    ! run the comparison
    call gaussian_compare(c1,p2,imode,diff,errmsg,seedout=seed,verbose0=.true.,&
       alpha0=alpha,lambda0=lambda,maxfeval0=maxfeval,besteps0=besteps,&
       max_elong_def0=max_elong,max_ang_def0=max_ang)

    if (len_trim(errmsg) > 0) &
       call ferror("struct_comparevc_vcgpwdf",errmsg,faterr)

    write (uout,'("+ DIFF = ",A/)') string(max(diff,0d0),'f',decimal=10)

    ! write to output
    if (dowrite .and. imode /= imode_sp) then
       call c1%struct_new(seed,.true.)

       ! write structure to output
       word = fileroot // "-final.res"
       call c1%write_simple_driver(word)
       write (uout,'("+ Final structure written to ",A/)') trim(word)

       ! write diffraction patterns to output
       word = fileroot // "-final.xy"
       lu = fopen_write(word)
       call c1%powder(0,th2ini,th2end,lambda,xrpd_fpol_def,final_npts,&
          xrpd_sigma_def,.false.,t=t,ih=ih)
       write (lu,'("# ",A,A)') string("2*theta",13,ioj_center), string("Intensity",13,ioj_center)
       do i = 1, final_npts
          write (lu,'(A," ",A)') string(t(i),"f",15,7,ioj_center), &
             string(ih(i),"f",15,7,ioj_center)
       end do
       call fclose(lu)
       write (uout,'("+ Final diffraction pattern written to ",A/)') trim(word)
    end if
contains
  subroutine set_quick_params()
    use crystalmod, only: xrpd_besteps_def_quick, xrpd_maxfeval_def_quick,&
       xrpd_max_ang_def_quick, xrpd_max_elong_def_quick

    max_elong = xrpd_max_elong_def_quick
    max_ang = xrpd_max_ang_def_quick
    maxfeval = xrpd_maxfeval_def_quick
    besteps = xrpd_besteps_def_quick

  end subroutine set_quick_params
  subroutine set_safe_params()
    use crystalmod, only: xrpd_besteps_def_safe, xrpd_maxfeval_def_safe,&
       xrpd_max_elong_def_safe, xrpd_max_ang_def_safe

    max_elong = xrpd_max_elong_def_safe
    max_ang = xrpd_max_ang_def_safe
    maxfeval = xrpd_maxfeval_def_safe
    besteps = xrpd_besteps_def_safe

  end subroutine set_safe_params
#endif

  end subroutine struct_comparevc_vcgpwdf

  !> Calculate the atomic environment of a point or all the
  !> non-equivalent atoms in the unit cell.
  module subroutine struct_environ(s,line)
    use systemmod, only: system
    use global, only: eval_next, dunit0, iunit, iunitname0
    use tools_io, only: string, lgetword, equal, ferror, faterr, string, uout,&
       ioj_right, ioj_center, zatguess, isinteger
    use param, only: bohrtoa, icrd_crys
    type(system), intent(inout), target :: s
    character*(*), intent(in) :: line

    real*8 :: up2d
    integer :: nat
    integer, allocatable :: eid(:), ishell(:), lvec(:,:)
    real*8, allocatable :: dist(:)
    integer :: lp, lp2
    character(len=:), allocatable :: word
    real*8 :: x0(3), x0in(3)
    logical :: doatoms, ok, groupshell
    logical, allocatable :: isdone(:)
    integer :: i, j, k, l, idx

    integer :: iat, iby, iat_mode, iby_mode
    integer, parameter :: inone = 0
    integer, parameter :: iznuc = 1
    integer, parameter :: iid = 2
    integer, parameter :: iidcel = 3

    up2d = 5d0 / bohrtoa

    x0 = 0d0
    doatoms = .true.
    iat = 0
    iat_mode = inone
    iby = 0
    iby_mode = inone
    groupshell = .false.

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"dist")) then
          ok = eval_next(up2d,line,lp)
          if (.not.ok) then
             call ferror('struct_environ','Wrong DIST syntax',faterr,line,syntax=.true.)
             return
          end if
          up2d = up2d / dunit0(iunit)
       elseif (equal(word,"point")) then
          ok = eval_next(x0(1),line,lp)
          ok = ok .and. eval_next(x0(2),line,lp)
          ok = ok .and. eval_next(x0(3),line,lp)
          if (.not.ok) then
             call ferror('struct_environ','Wrong POINT syntax',faterr,line,syntax=.true.)
             return
          end if
          doatoms = .false.
          x0in = x0
          if (s%c%ismolecule) &
             x0 = s%c%c2x(x0 / dunit0(iunit) - s%c%molx0)
       elseif (equal(word,"atom")) then
          doatoms = .true.
          lp2 = lp
          word = lgetword(line,lp)
          iat = zatguess(word)
          if (iat < 0) then
             lp = lp2
             ok = isinteger(iat,line,lp)
             if (.not.ok) &
                call ferror('struct_environ','Syntax error in ENVIRON/ATOM',faterr,line,syntax=.true.)
             iat_mode = iid
          else
             iat_mode = iznuc
          end if
       elseif (equal(word,"celatom")) then
          ok = isinteger(iat,line,lp)
          if (.not.ok) &
             call ferror('struct_environ','Syntax error in ENVIRON/CELATOM',faterr,line,syntax=.true.)
          iat_mode = iidcel
          doatoms = .true.
       elseif (equal(word,"by")) then
          lp2 = lp
          word = lgetword(line,lp)
          iby = zatguess(word)
          if (iby < 0) then
             lp = lp2
             ok = isinteger(iby,line,lp)
             if (.not.ok) &
                call ferror('struct_environ','Syntax error in ENVIRON/BY',faterr,line,syntax=.true.)
             iby_mode = iid
          else
             iby_mode = iznuc
          end if
       elseif (equal(word,"shell").or.equal(word,"shells")) then
          groupshell = .true.
       elseif (len_trim(word) > 0) then
          call ferror('struct_environ','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    write (uout,'("* ENVIRON")')
    if (doatoms) then
       allocate(isdone(s%c%nneq))
       isdone = .false.
       do i = 1, s%c%ncel
          idx = s%c%atcel(i)%idx
          if (iat_mode == iid .and. idx /= iat .or. isdone(idx)) cycle
          if (iat_mode == iznuc .and. s%c%spc(s%c%at(idx)%is)%z /= iat .or. isdone(idx)) cycle
          if (iat_mode == iidcel .and. i /= iat) cycle
          isdone(idx) = .true.

          if (iby_mode == iid) then
             call s%c%list_near_atoms(s%c%atcel(i)%x,icrd_crys,.true.,nat,eid,dist,&
                lvec,ishell,up2d=up2d,nid0=iby,nozero=.true.)
          elseif (iby_mode == iznuc) then
             call s%c%list_near_atoms(s%c%atcel(i)%x,icrd_crys,.true.,nat,eid,dist,&
                lvec,ishell,up2d=up2d,iz0=iby,nozero=.true.)
          else
             call s%c%list_near_atoms(s%c%atcel(i)%x,icrd_crys,.true.,nat,eid,dist,&
                lvec,ishell,up2d=up2d,nozero=.true.)
          end if

          write (uout,'("+ Environment of atom ",A," (spc=",A,", nid=",A,") at ",3(A," "))') &
             string(i), string(s%c%spc(s%c%at(idx)%is)%name), string(idx), &
             (string(s%c%atcel(i)%x(j),'f',length=10,decimal=6),j=1,3)

          if (groupshell) then
             call output_by_shell()
          else
             call output_by_distance()
          end if
          if (iat_mode == iid) exit
       end do
       deallocate(isdone)
    else
       if (iby_mode == iid) then
          call s%c%list_near_atoms(x0,icrd_crys,.true.,nat,eid,dist,lvec,ishell,up2d=up2d,nid0=iby)
       elseif (iby_mode == iznuc) then
          call s%c%list_near_atoms(x0,icrd_crys,.true.,nat,eid,dist,lvec,ishell,up2d=up2d,iz0=iby)
       else
          call s%c%list_near_atoms(x0,icrd_crys,.true.,nat,eid,dist,lvec,ishell,up2d=up2d)
       end if

       write (uout,'("+ Environment of point ",3(A," "))') (string(x0in(j),'f',length=10,decimal=6),j=1,3)

       if (groupshell) then
          call output_by_shell()
       else
          call output_by_distance()
       end if
    end if

  contains
    subroutine output_by_distance()
      integer :: mshel, cidx, nidx
      real*8 :: xx(3)

      mshel = maxval(ishell)
      write (uout,'("# Up to distance (",A,"): ",A)') iunitname0(iunit), string(up2d*dunit0(iunit),'f',length=10,decimal=6)
      write (uout,'("# Number of atoms in the environment: ",A)') string(nat)
      write (uout,'("# Number of atomic shells in the environment: ",A)') string(mshel)
      write (uout,'("# nid = non-equivalent list atomic ID. id = complete list atomic ID plus lattice vector (lvec).")')
      write (uout,'("# name = atomic name (symbol). spc = atomic species. dist = distance. shel = shell ID.")')
      if (s%c%ismolecule) then
         write (uout,'("#nid   id      lvec     name   spc  dist(",A,") shel          Coordinates (",A,") ")') &
            iunitname0(iunit), iunitname0(iunit)
      else
         write (uout,'("#nid   id      lvec     name   spc  dist(",A,") shel      Coordinates (cryst. coord.) ")') iunitname0(iunit)
      end if
      do j = 1, nat
         cidx = eid(j)
         nidx = s%c%atcel(cidx)%idx
         if (s%c%ismolecule) then
            xx = (s%c%atcel(cidx)%r + s%c%molx0) * dunit0(iunit)
         else
            xx = s%c%atcel(cidx)%x + lvec(:,j)
         end if

         write (uout,'("  ",2(A," "),"(",A," ",A," ",A,")",99(" ",A))') string(nidx,4,ioj_center), string(cidx,4,ioj_center),&
            (string(lvec(k,j),2,ioj_right),k=1,3), string(s%c%spc(s%c%atcel(cidx)%is)%name,7,ioj_center),&
            string(s%c%atcel(cidx)%is,2,ioj_right), string(dist(j)*dunit0(iunit),'f',12,6,4), string(ishell(j),3,ioj_center),&
            (string(xx(k),'f',12,8,4),k=1,3)
      end do
      write (uout,*)
    end subroutine output_by_distance

    subroutine output_by_shell()
      integer :: mshel, cidx, nidx, nneig, ishl0
      real*8 :: xx(3)

      mshel = maxval(ishell)
      write (uout,'("# Up to distance (",A,"): ",A)') iunitname0(iunit), string(up2d*dunit0(iunit),'f',length=10,decimal=6)
      write (uout,'("# Number of atoms in the environment: ",A)') string(nat)
      write (uout,'("# Number of atomic shells in the environment: ",A)') string(mshel)
      write (uout,'("# ishl = shell ID. nneig = number of neighbors in the shell. nid = non-equivalent list atomic ID.")')
      write (uout,'("# name = atomic name (symbol). spc = atomic species. dist = distance. Rest of fields are for a single")')
      write (uout,'("# representative of the shell: id = complete list atomic ID plus lattice vector (lvec) and coordinates. ")')
      if (s%c%ismolecule) then
         write (uout,'("#ishl nneig nid   name   spc  dist(",A,")  id     lvec            Coordinates (",A,") ")') &
            iunitname0(iunit), iunitname0(iunit)
      else
         write (uout,'("#ishl nneig nid   name   spc  dist(",A,")  id     lvec        Coordinates (cryst. coord.) ")')&
            iunitname0(iunit)
      end if
      k = 1
      main: do j = 1, mshel
         ishl0 = k
         do while (ishell(k) == j)
            k = k + 1
            if (k > nat) exit
         end do
         if (k == 1) exit
         nneig = k - ishl0

         cidx = eid(k-1)
         nidx = s%c%atcel(cidx)%idx
         if (s%c%ismolecule) then
            xx = (s%c%atcel(cidx)%r + s%c%molx0) * dunit0(iunit)
         else
            xx = s%c%atcel(cidx)%x + lvec(:,k-1)
         end if

         write (uout,'(7(A," "),"(",3(A," "),")",99(A," "))') string(j,5,ioj_center),&
            string(nneig,5,ioj_center), string(nidx,4,ioj_center),&
            string(s%c%spc(s%c%atcel(cidx)%is)%name,7,ioj_center), string(s%c%atcel(cidx)%is,2,ioj_right),&
            string(dist(k-1)*dunit0(iunit),'f',12,6,4), string(cidx,4,ioj_center),&
            (string(lvec(l,k-1),2,ioj_right),l=1,3), (string(xx(l),'f',12,8,4),l=1,3)
      end do main
      write (uout,*)

    end subroutine output_by_shell

  end subroutine struct_environ

  !> Calculate the coordination numbers of all atoms.
  module subroutine struct_coord(s,line)
    use systemmod, only: system
    use global, only: bondfactor, eval_next
    use tools_io, only: ferror, faterr, zatguess, lgetword, equal, isinteger, ioj_center,&
       ioj_left, ioj_right, uout, string
    use param, only: atmcov, icrd_crys
    type(system), intent(inout), target :: s
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word
    integer :: lp, i, j, k
    real*8, allocatable :: rad(:)
    real*8 :: fac, dd
    logical :: ok
    integer :: nat
    integer, allocatable :: eid(:), coord2(:,:), coord2sp(:,:), coord3(:,:,:), coord3sp(:,:,:)
    real*8, allocatable :: up2dsp(:,:)

    ! allocate the default radii
    fac = bondfactor
    allocate(rad(s%c%nspc))
    rad = 0d0
    do i = 1, s%c%nspc
       if (s%c%spc(i)%z > 0) then
          rad(i) = atmcov(s%c%spc(i)%z)
       end if
    end do

    ! parse the input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"dist")) then
          ok = eval_next(dd,line,lp)
          if (.not.ok) then
             call ferror('struct_coord','Wrong DIST keyword',faterr,line,syntax=.true.)
             return
          end if
          fac = 1d0
          rad = 0.5d0 * dd
       elseif (equal(word,"fac")) then
          ok = eval_next(fac,line,lp)
          if (.not.ok) then
             call ferror('struct_coord','Wrong FAC keyword',faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror('struct_coord','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! some output
    write (uout,'("* COORD: calculation of coordination numbers")')
    write (uout,'("+ Atomic radii for all species")')
    write (uout,'("# Atoms i and j are coordinated if distance(i,j) <= fac*(radi+radj), fac = ",A)') &
       string(fac,'f',length=5,decimal=2,justify=ioj_left)
    write (uout,'("# ",99(A," "))') string("spc",3,ioj_center), &
       string("Z",3,ioj_center), string("name",7,ioj_center),&
       string("rad",length=5,justify=ioj_center)
    do i = 1, s%c%nspc
       write (uout,'("  ",99(A," "))') string(i,3,ioj_center), &
          string(s%c%spc(i)%z,3,ioj_center), string(s%c%spc(i)%name,7,ioj_center),&
          string(rad(i),'f',length=5,decimal=2,justify=ioj_right)
    end do
    write (uout,*)

    ! calculate the coordination numbers (per atom)
    allocate(up2dsp(s%c%nspc,2),coord2(s%c%nneq,s%c%nspc))
    up2dsp = 0d0
    coord2 = 0
    do i = 1, s%c%nneq
       up2dsp(:,2) = rad + rad(s%c%at(i)%is)
       call s%c%list_near_atoms(s%c%at(i)%x,icrd_crys,.false.,nat,eid,up2dsp=up2dsp,nozero=.true.)
       do j = 1, nat
          coord2(i,s%c%atcel(eid(j))%is) = coord2(i,s%c%atcel(eid(j))%is) + 1
       end do
    end do
    deallocate(rad,up2dsp)

    ! output coordination numbers (per non-equivalent atom)
    write (uout,'("+ Pair coordination numbers (per non-equivalent atom in the cell)")')
    write (uout,'("# id = non-equivalent atom ID. mult = multiplicity. rest = atomic species (by name).")')
    write (uout,'("# ",99(A," "))') string("nid",4,ioj_center), string("name",7,ioj_center), &
       string("mult",5,ioj_center), (string(s%c%spc(j)%name,5,ioj_left),j=1,s%c%nspc), "total"
    do i = 1, s%c%nneq
       write (uout,'("  ",99(A," "))') string(i,4,ioj_center), string(s%c%spc(s%c%at(i)%is)%name,7,ioj_center), &
          string(s%c%at(i)%mult,5,ioj_center), (string(coord2(i,j),5,ioj_center),j=1,s%c%nspc),&
          string(sum(coord2(i,1:s%c%nspc)),5,ioj_center)
    end do
    write (uout,*)

    ! calculate the coordination numbers (per species)
    allocate(coord2sp(s%c%nspc,s%c%nspc))
    coord2sp = 0
    do i = 1, s%c%nneq
       do j = 1, s%c%nspc
          coord2sp(s%c%at(i)%is,j) = coord2sp(s%c%at(i)%is,j) + coord2(i,j) * s%c%at(i)%mult
       end do
    end do

    write (uout,'("+ Pair coordination numbers (per species)")')
    write (uout,'("# ",99(A," "))') string("spc",8,ioj_center), (string(s%c%spc(j)%name,5,ioj_center),j=1,s%c%nspc)
    do i = 1, s%c%nspc
       write (uout,'("  ",99(A," "))') string(i,3,ioj_center), string(s%c%spc(i)%name,5,ioj_center), &
          (string((coord2sp(i,j)+coord2sp(j,i))/2,5,ioj_center),j=1,s%c%nspc)
    end do
    write (uout,*)
    deallocate(coord2sp)

    ! calculate the number of triplets
    allocate(coord3(s%c%nneq,s%c%nspc,s%c%nspc))
    coord3 = 0
    do i = 1, s%c%nneq
       do j = 1, s%c%nspc
          do k = 1, s%c%nspc
             if (j == k) then
                coord3(i,j,k) = coord3(i,j,k) + coord2(i,j) * (coord2(i,k)-1) / 2
             else
                coord3(i,j,k) = coord3(i,j,k) + coord2(i,j) * coord2(i,k)
             end if
          end do
       end do
    end do

    write (uout,'("+ Triplet coordination numbers (per non-equivalent atom in the cell)")')
    write (uout,'("# Y runs over atoms in the cell: nid = non-equivalent atom ID. mult = multiplicity.")')
    write (uout,'("#  ---    Y    ---     X     Z")')
    write (uout,'("# nid   name   mult   spc   spc  X-Y-Z")')
    do i = 1, s%c%nneq
       do j = 1, s%c%nspc
          do k = j, s%c%nspc
             write (uout,'("  ",99(A," "))') string(i,4,ioj_center), string(s%c%spc(s%c%at(i)%is)%name,7,ioj_center), &
                string(s%c%at(i)%mult,5,ioj_center), string(s%c%spc(j)%name,5,ioj_center), string(s%c%spc(k)%name,5,ioj_center),&
                string((coord3(i,j,k)+coord3(i,k,j))/2,5,ioj_center)
          end do
       end do
    end do
    write (uout,*)

    ! calculate the triplet coordination numbers (per species)
    allocate(coord3sp(s%c%nspc,s%c%nspc,s%c%nspc))
    coord3sp = 0
    do i = 1, s%c%nneq
       do j = 1, s%c%nspc
          do k = 1, s%c%nspc
             coord3sp(s%c%at(i)%is,j,k) = coord3sp(s%c%at(i)%is,j,k) + coord3(i,j,k) * s%c%at(i)%mult
          end do
       end do
    end do

    write (uout,'("+ Triplet coordination numbers (per species)")')
    write (uout,'("#   Y     X     Z   X-Y-Z")')
    do i = 1, s%c%nspc
       do j = 1, s%c%nspc
          do k = j, s%c%nspc
             write (uout,'("  ",99(A," "))') string(s%c%spc(j)%name,5,ioj_center), string(s%c%spc(i)%name,5,ioj_center), &
                string(s%c%spc(k)%name,5,ioj_center), string((coord3sp(i,j,k)+coord3sp(i,k,j))/2,5,ioj_center)
          end do
       end do
    end do
    write (uout,*)

    deallocate(coord3,coord3sp)

  end subroutine struct_coord

  !> Calculate the coordination polyedra.
  module subroutine struct_polyhedra(s,line)
    use iso_c_binding, only: c_int, c_double, c_ptr
    use systemmod, only: system
    use tools_io, only: equali, zatguess, ferror, faterr, getword, uout, string, ioj_left, ioj_center
    use tools_math, only: mixed
    use global, only: bondfactor, eval_next, dunit0, iunit, iunitname0
    ! use tools_io, only: ferror, faterr, zatguess, lgetword, equal, isinteger, ioj_center,&
    !    ioj_left, ioj_right, uout, string
    use param, only: atmcov, icrd_crys
    type(system), intent(inout), target :: s
    character*(*), intent(in) :: line

    character(len=:), allocatable :: at1, at2
    integer :: lp, is1, is2, iz1, iz2
    real*8 :: rdum, rmin, rmax, vol, xp1(3), xp2(3), xp3(3)
    logical :: ok
    integer, allocatable :: eid(:)
    integer :: i, j
    real*8, allocatable :: dist(:)
    integer(c_int) :: nat, nf
    real(c_double) :: x0(3)
    real(c_double), allocatable :: xstar(:,:)
    integer(c_int), allocatable :: iface(:,:), lvec(:,:)
    type(c_ptr), target :: fid

    interface
       ! The definitions and documentation for these functions are in doqhull.c
       subroutine runqhull_basintriangulate_step1(n,x0,xvert,nf,fid) bind(c)
         import c_int, c_double, c_ptr
         integer(c_int), value :: n
         real(c_double) :: x0(3)
         real(c_double) :: xvert(3,n)
         integer(c_int) :: nf
         type(c_ptr) :: fid
       end subroutine runqhull_basintriangulate_step1
       subroutine runqhull_basintriangulate_step2(nf,iface,fid) bind(c)
         import c_int, c_double, c_ptr
         integer(c_int), value :: nf
         integer(c_int) :: iface(3,nf)
         type(c_ptr), value :: fid
       end subroutine runqhull_basintriangulate_step2
    end interface

    ! Initialize
    lp = 1

    ! Literal interpretation of the species
    at1 = getword(line,lp)
    at2 = getword(line,lp)
    is1 = 0
    is2 = 0
    do i = 1, s%c%nspc
       if (is1 == 0) then
          if (equali(at1,s%c%spc(i)%name)) &
             is1 = i
       end if
       if (is2 == 0) then
          if (equali(at2,s%c%spc(i)%name)) &
             is2 = i
       end if
    end do

    ! If not found, convert to atomic number and try that instead
    if (is1 == 0) then
       iz1 = zatguess(at1)
    else
       iz1 = s%c%spc(is1)%z
    end if
    if (is2 == 0) then
       iz2 = zatguess(at2)
    else
       iz2 = s%c%spc(is2)%z
    end if
    if (iz1 <= 0) &
       call ferror("struct_polyhedra","unrecognized atomic species: " // at1,faterr)
    if (iz2 <= 0) &
       call ferror("struct_polyhedra","unrecognized atomic species: " // at2,faterr)

    ! default distance range
    rmin = 0d0
    rmax = (atmcov(iz1) + atmcov(iz2)) * bondfactor
    ok = eval_next(rdum,line,lp)
    if (ok) rmax = rdum / dunit0(iunit)
    ok = eval_next(rdum,line,lp)
    if (ok) then
       rmin = rmax
       rmax = rdum / dunit0(iunit)
    end if

    write (uout,'("* POLYHEDRA: calculation of coordination polyhedra")')
    if (is1 == 0) then
       write (uout,'("  Center of the polyhedra: Z = ",A)') string(iz1)
    else
       write (uout,'("  Center of the polyhedra: ",A)') trim(s%c%spc(is1)%name)
    end if
    if (is2 == 0) then
       write (uout,'("  Vertex of the polyhedra: Z = ",A)') string(iz2)
    else
       write (uout,'("  Vertex of the polyhedra: ",A)') trim(s%c%spc(is2)%name)
    end if
    write (uout,'("  Distance range(",A,"): ",A," to ",A)') string(iunitname0(iunit)),&
       string(rmin*dunit0(iunit),'f',decimal=4), string(rmax*dunit0(iunit),'f',decimal=4)
    write (uout, '("# nid = non-equivalent atom ID. at = atomic name (symbol)")')
    write (uout, '("# cryst coords = crystallographic coordinates. nv = number of vertices.")')
    write (uout, '("# rmin = minimum vertex distance (",A,"), rmax = maximum vertex distance (",A,").")') &
       string(iunitname0(iunit)), string(iunitname0(iunit))
    write (uout, '("# nf = number of faces. volume = polyhedron volume (",A,"^3).")') string(iunitname0(iunit))
    write (uout,'("#nid at mult   ----- Cryst. coords. -----     nv   rmin       rmax      nf   volume")')
    do i = 1, s%c%nneq
       if (is1 == 0) then
          if (s%c%spc(s%c%at(i)%is)%z /= iz1) cycle
       else
          if (s%c%at(i)%is /= is1) cycle
       end if

       ! find the environment
       if (is2 /= 0) then
          ! by species
          call s%c%list_near_atoms(s%c%at(i)%x,icrd_crys,.true.,nat,eid=eid,dist=dist,&
             lvec=lvec,up2d=rmax,ispc0=is2,nozero=.true.)
       else
          ! by z
          call s%c%list_near_atoms(s%c%at(i)%x,icrd_crys,.true.,nat,eid=eid,dist=dist,&
             lvec=lvec,up2d=rmax,iz0=iz2,nozero=.true.)
       end if
       if (nat <= 2) cycle

       ! project on a sphere and triangulate the convex polyhedron
       if (allocated(xstar)) deallocate(xstar)
       allocate(xstar(3,nat))
       x0 = s%c%at(i)%r
       do j = 1, nat
          xstar(:,j) = s%c%atcel(eid(j))%x + lvec(:,j)
          xstar(:,j) = s%c%x2c(xstar(:,j))
       end do
       call runqhull_basintriangulate_step1(nat,x0,xstar,nf,fid)
       if (allocated(iface)) deallocate(iface)
       allocate(iface(3,nf))
       call runqhull_basintriangulate_step2(nf,iface,fid)

       ! calculate the polyhedron volume
       vol = 0d0
       do j = 1, nf
          xp1 = xstar(:,iface(1,j)) - x0
          xp2 = xstar(:,iface(2,j)) - x0
          xp3 = xstar(:,iface(3,j)) - x0
          vol = vol + abs(mixed(xp1,xp2,xp3)) / 6d0
       end do

       write (uout,'(99(A," "))') string(i,3,ioj_left), &
          string(s%c%spc(s%c%at(i)%is)%name,4,ioj_center),string(s%c%at(i)%mult,3),&
          (string(s%c%at(i)%x(j),'f',length=10,decimal=6),j=1,3),&
          string(nat,3,ioj_center), string(minval(dist)*dunit0(iunit),'f',length=10,decimal=6),&
          string(maxval(dist)*dunit0(iunit),'f',length=10,decimal=6), string(nf,3,ioj_center),&
          string(vol*dunit0(iunit)**3,'f',length=10,decimal=6)
    end do
    write (uout,*)

  end subroutine struct_polyhedra

  !> Calculate the packing ratio of the crystal.
  module subroutine struct_packing(s,line)
    use systemmod, only: system
    use global, only: eval_next
    use tools_io, only: ferror, faterr, uout, lgetword, equal, string
    use param, only: atmvdw, atmcov
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line

    integer :: lp
    logical :: ok
    character(len=:), allocatable :: word
    real*8 :: prec
    real*8 :: vout
    integer :: whichuse = 0

    if (s%c%ismolecule) then
       call ferror("critic2","PACKING cannot be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! default values
    prec = 1d-2
    whichuse = 0 ! 0 = nnm, 1 = vdw, 2 = cov

    ! header
    write (uout,'("* PACKING")')

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"vdw")) then
          whichuse = 1
       elseif (equal(word,"cov")) then
          whichuse = 2
       elseif (equal(word,"prec")) then
          ok = eval_next(prec,line,lp)
          if (.not.ok) then
             call ferror('struct_packing','Wrong PREC syntax',faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror('struct_packing','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    if (whichuse == 0) then
       write (uout,'("+ Packing ratio (%): ",A)') string(s%c%get_pack_ratio(),'f',10,4)
    else
       if (whichuse == 1) then
          vout = s%c%vdw_volume(prec,atmvdw)
       else
          vout = s%c%vdw_volume(prec,atmcov)
       end if
       write (uout,'("+ Van der Waals volume: ",A," +- ",A)') &
          string(vout,'f',decimal=6), string(prec*vout,'f',decimal=6)
       write (uout,'("+ Interstitial volume (outside vdw spheres): ",A," +- ",A)') &
          string(s%c%omega-vout,'f',decimal=6), string(prec*vout,'f',decimal=6)
       write (uout,'("+ Cell volume: ",A)') string(s%c%omega,'f',decimal=6)
       write (uout,'("+ Packing ratio (%): ",A," +- ",A)') string(vout/s%c%omega*100d0,'f',decimal=6), &
          string(prec*vout/s%c%omega*100d0,'f',decimal=6)
    end if
    write (uout,*)

  end subroutine struct_packing

  !> Calculate the van der Waals volume of a crystal or molecule.
  module subroutine struct_vdw(s,line)
    use global, only: iunitname0, dunit0, iunit
    use systemmod, only: system
    use tools_io, only: ferror, faterr, uout, lgetword, equal, string
    use global, only: eval_next
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line

    integer :: lp
    logical :: ok
    character(len=:), allocatable :: word
    real*8 :: vvdw, fac3
    real*8 :: prec

    ! default values
    prec = 1d-2

    ! header
    write (uout,'("* VDW")')

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"prec")) then
          ok = eval_next(prec,line,lp)
          if (.not.ok) then
             call ferror('struct_packing','Wrong PREC syntax in VDW',faterr,line,syntax=.true.)
             return
          end if
       elseif (len_trim(word) > 0) then
          call ferror('struct_packing','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    ! calculate vdw volume
    fac3 = dunit0(iunit)**3
    vvdw = s%c%vdw_volume(prec) * fac3

    ! output
    write (uout,'("+ Requested relative std. deviation (sigma_{Vvdw}/Vvdw): ",A)') string(prec,'e',decimal=4)
    write (uout,'("+ Van der Waals volume (Vvdw): ",A," +- ",A," ",A,"^3")') &
       string(vvdw,'f',decimal=6), string(prec*vvdw,'f',decimal=6),&
       iunitname0(iunit)
    if (.not.s%c%ismolecule) then
       write (uout,'("+ Interstitial volume (outside vdw spheres): ",A," +- ",A," ",A,"^3")') &
          string(s%c%omega*fac3-vvdw,'f',decimal=6), string(prec*vvdw,'f',decimal=6),&
          iunitname0(iunit)
       write (uout,'("+ Cell volume: ",A," ",A,"^3")') string(s%c%omega*fac3,'f',decimal=6),&
          iunitname0(iunit)
       write (uout,'("+ Packing ratio (%): ",A," +- ",A)') string(vvdw/s%c%omega/fac3*100d0,'f',decimal=6), &
          string(prec*vvdw/s%c%omega/fac3*100d0,'f',decimal=6)
    end if
    write (uout,*)

  end subroutine struct_vdw

  !> Edit a molecular or crystal structure
  module subroutine struct_edit(s,verbose)
    use systemmod, only: system
    use tools_io, only: uout, ucopy, uin, getline, lgetword, equal, ferror, faterr,&
       getword, isinteger, string
    use global, only: eval_next_real, iunit, iunit_fractional, iunit_bohr, iunit_ang
    use types, only: realloc
    type(system), intent(inout) :: s
    logical, intent(in) :: verbose

    character(len=:), allocatable :: line, word, errmsg
    integer :: lp
    integer :: idx, nat, i, iunit_l, iaxis
    integer, allocatable :: iat(:)
    logical :: ok, dorelative, dofraction, changed
    real*8 :: x(3), rdum

    if (verbose) then
       write (uout,'("* Edit the crystal or molecular structure (EDIT)")')
       write (uout,*)
    end if

    ! transform to the primitive?
    errmsg = "Invalid syntax"
    changed = .false.
    do while (getline(uin,line,ucopy=ucopy))
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"delete")) then
          word = lgetword(line,lp)

          nat = 0
          if (allocated(iat)) deallocate(iat)
          if (equal(word,"atoms") .or. equal(word,"atom")) then
             ! read the list of atoms
             allocate(iat(10))
             do while (.true.)
                word = getword(line,lp)
                if (len_trim(word) == 0) exit
                ok = isinteger(idx,word)
                if (.not.ok) then
                   errmsg = "invalid atom number in delete atoms"
                   goto 999
                end if
                if (idx < 1 .or. idx > s%c%ncel) then
                   errmsg = "atomic ID " // string(idx) // " outside valid range"
                   goto 999
                end if
                nat = nat + 1
                if (nat > size(iat,1)) call realloc(iat,2*nat)
                iat(nat) = idx
             end do
          elseif (equal(word,"molecules") .or. equal(word,"molecule")) then
             ! read the list of molecules
             allocate(iat(10))
             do while (.true.)
                word = getword(line,lp)
                if (len_trim(word) == 0) exit
                ok = isinteger(idx,word)
                if (.not.ok) then
                   errmsg = "invalid molecule number in delete atoms"
                   goto 999
                end if
                if (idx < 1 .or. idx > s%c%nmol) then
                   errmsg = "molecule ID " // string(idx) // " outside valid range"
                   goto 999
                end if
                do i = 1, s%c%ncel
                   if (s%c%idatcelmol(1,i) == idx) then
                      nat = nat + 1
                      if (nat > size(iat,1)) call realloc(iat,2*nat)
                      iat(nat) = i
                   end if
                end do
             end do
          elseif (equal(word,"hydrogen") .or. equal(word,"hydrogens")) then
             allocate(iat(s%c%ncel))
             nat = 0
             do i = 1, s%c%ncel
                if (s%c%spc(s%c%atcel(i)%is)%z == 1) then
                   nat = nat + 1
                   iat(nat) = i
                end if
             end do
          else
             goto 999
          end if
          if (nat > 0) then
             ! transform
             call s%c%delete_atoms(nat,iat(1:nat))
             changed = .true.
          end if
       elseif (equal(word,"move")) then
          ! read atom ID and position
          ok = isinteger(idx,line,lp)
          if (.not.ok) then
             errmsg = "invalid atom number in move"
             goto 999
          end if
          if (idx < 1 .or. idx > s%c%ncel) then
             errmsg = "atomic ID " // string(idx) // " outside valid range"
             goto 999
          end if
          ok = eval_next_real(x(1),line,lp)
          ok = ok .and. eval_next_real(x(2),line,lp)
          ok = ok .and. eval_next_real(x(3),line,lp)

          ! parse the rest of the options
          if (s%c%ismolecule) then
             iunit_l = iunit
          else
             iunit_l = iunit_fractional
          end if
          dorelative = .false.
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,"relative")) then
                dorelative = .true.
             elseif (equal(word,"bohr")) then
                iunit_l = iunit_bohr
             elseif (equal(word,"ang").or.equal(word,"angstrom")) then
                iunit_l = iunit_ang
             elseif (len_trim(word) == 0) then
                exit
             else
                errmsg = "unknown extra keyword: " // word
                goto 999
             end if
          end do

          ! transform
          call s%c%move_atom(idx,x,iunit_l,dorelative)
          changed = .true.

       elseif (equal(word,"cellmove")) then
          ! read axis/volume
          word = lgetword(line,lp)
          if (equal(word,"a")) then
             iaxis = 1
          elseif (equal(word,"b")) then
             iaxis = 2
          elseif (equal(word,"c")) then
             iaxis = 3
          elseif (equal(word,"alpha")) then
             iaxis = -1
          elseif (equal(word,"beta")) then
             iaxis = -2
          elseif (equal(word,"gamma")) then
             iaxis = -3
          elseif (equal(word,"volume").or.equal(word,"vol").or.equal(word,"v")) then
             iaxis = 0
          else
             errmsg = "unknown keyword in cellmove: " // word
             goto 999
          end if

          ! read displacement
          ok = eval_next_real(rdum,line,lp)
          if (.not.ok) then
             errmsg = "invalid numerical argument in cellmove"
             goto 999
          end if

          ! parse the rest of the options
          dofraction = .false.
          dorelative = .false.
          iunit_l = iunit
          do while (.true.)
             word = lgetword(line,lp)
             if (equal(word,"relative")) then
                dorelative = .true.
             elseif (equal(word,"fraction")) then
                dofraction = .true.
             elseif (equal(word,"bohr")) then
                iunit_l = iunit_bohr
             elseif (equal(word,"ang").or.equal(word,"angstrom")) then
                iunit_l = iunit_ang
             elseif (len_trim(word) == 0) then
                exit
             else
                errmsg = "unknown extra keyword: " // word
                goto 999
             end if
          end do

          ! transform
          call s%c%move_cell(iaxis,rdum,iunit_l,dorelative,dofraction)
          changed = .true.

       elseif (equal(word,"end") .or. equal(word,"endedit")) then
          exit
       elseif (len_trim(word) > 0) then
          goto 999
       end if
    end do

    ! wrap up
    if (changed) then
       call s%reset_fields()
       if (verbose) &
          call s%report(.true.,.true.,.true.,.true.,.true.,.true.,.false.)
    end if

    return
999 continue

    call ferror('struct_edit',errmsg,faterr,line,syntax=.true.)

  end subroutine struct_edit

  !> Build a new crystal from the current crystal by cell transformation
  module subroutine struct_newcell(s,line,verbose)
    use systemmod, only: system
    use tools_math, only: matinv, cross, det3
    use global, only: iunitname0, dunit0, iunit, eval_next
    use tools_io, only: uout, ferror, faterr, lgetword, equal, string, isinteger, ioj_left,&
       ioj_right
    use tools, only: delaunay_reduction
    use param, only: icrd_crys
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line
    logical, intent(in) :: verbose

    character(len=:), allocatable :: word
    logical :: ok, doprim, doforce, donice, changed
    integer :: lp, lp2, dotyp, i, j, k, n, inice
    real*8 :: x0(3,3), t0(3), rdum(4), x2c(3,3), dd
    real*8 :: r, mm(3,3), dmax0
    logical :: doinv, dorefine
    integer :: nat
    integer, allocatable :: lvec(:,:)
    real*8, allocatable :: rmax(:), mmax(:,:,:)

    integer, parameter :: inice_def = 64

    if (s%c%ismolecule) then
       call ferror("struct_newcell","NEWCELL cannot be used with molecules",faterr,syntax=.true.)
       return
    end if

    write (uout,'("* Transformation to a new unit cell (NEWCELL)")')
    write (uout,*)

    ! transform to the primitive?
    lp = 1
    doprim = .false.
    doforce = .false.
    dorefine = .false.
    donice = .false.
    inice = inice_def
    dotyp = 0
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"standard")) then
          dotyp = 1
          doprim = .false.
       elseif (equal(word,"primitive")) then
          doforce = .false.
          dotyp = 1
          doprim = .true.
       elseif (equal(word,"primstd")) then
          doforce = .true.
          dotyp = 1
          doprim = .true.
       elseif (equal(word,"niggli")) then
          dotyp = 2
       elseif (equal(word,"delaunay")) then
          dotyp = 3
       elseif (equal(word,"refine")) then
          dorefine = .true.
       elseif (equal(word,"nice")) then
          donice = .true.
          lp2 = lp
          ok = isinteger(inice,line,lp)
          if (.not.ok) then
             lp = lp2
             inice = inice_def
          end if
          inice = max(inice,1)
       else
          lp = 1
          exit
       end if
    end do

    ! search for a nice cell
    if (donice) then
       allocate(rmax(inice),mmax(3,3,inice))

       dmax0 = 1.2d0 * real(inice,8)**(1d0/3d0) * max(norm2(s%c%m_xr2c(:,1)),norm2(s%c%m_xr2c(:,2)),&
          norm2(s%c%m_xr2c(:,3))) + 1d-1
       call s%c%list_near_lattice_points((/0d0,0d0,0d0/),icrd_crys,.true.,nat,&
          lvec=lvec,up2d=dmax0,nozero=.true.)

       rmax = 0d0
       mmax = 0d0
       do i = 1, nat
          mm(:,1) = lvec(:,i)
          do j = i+1, nat
             mm(:,2) = lvec(:,j)
             do k = j+1, nat
                mm(:,3) = lvec(:,k)

                dd = det3(mm)
                if (dd < 0d0) mm = -mm
                dd = abs(dd)
                if (dd < 1d-2 .or. dd > (inice+0.5d0)) cycle
                n = nint(dd)
                if (n > inice) cycle

                x2c = matmul(s%c%m_x2c,mm)

                r = 0.5d0 * n * s%c%omega / max(norm2(cross(x2c(:,1),x2c(:,2))),&
                   norm2(cross(x2c(:,1),x2c(:,3))),norm2(cross(x2c(:,2),x2c(:,3))))
                if (r > rmax(n)) then
                   rmax(n) = r
                   mmax(:,:,n) = mm
                end if
             end do
          end do
       end do

       write (uout,'("+ List of nice cells with increasing size")')
       write (uout,'("# n = size of the supercell is n times the input cell.")')
       write (uout,'("# rmax = radius of the largest sphere that fits in the supercell (",A,").")') &
          string(iunitname0(iunit))
       write (uout,'("# niceness = inverse of the cell skewness, higher is nicer, cubic cell = 1")')
       write (uout,'("# newcell transformation = use these parameters in a NEWCELL command to obtain this cell")')
       write (uout,'("#n   rmax   niceness  -- NEWCELL transformation --")')
       do i = 1, inice
          write (uout,'(3(A," ")," ",3(3(A," ")," "))') string(i,3,ioj_left),&
             string(rmax(i)*dunit0(iunit),'f',7,3), string(8d0 * rmax(i)**3 / (i*s%c%omega),'f',7,5),&
             ((string(nint(mmax(j,k,i)),length=2,justify=ioj_right),j=1,3),k=1,3)
       end do
       write (uout,*)
       return
    end if

    ! process the cell transformation
    t0 = 0d0
    x0 = 0d0
    if (dotyp == 1) then
       x0 = s%c%cell_standard(doprim,doforce,dorefine)
       changed = any(abs(x0) > 1d-5)
    elseif (dotyp == 2) then
       x0 = s%c%cell_niggli()
       changed = any(abs(x0) > 1d-5)
    elseif (dotyp == 3) then
       x0 = s%c%cell_delaunay()
       changed = any(abs(x0) > 1d-5)
    else
       ! read the vectors from the input
       ok = eval_next(rdum(1),line,lp)
       ok = ok .and. eval_next(rdum(2),line,lp)
       ok = ok .and. eval_next(rdum(3),line,lp)
       if (.not.ok) then
          call ferror("struct_newcell","Wrong syntax in NEWCELL",faterr,line,syntax=.true.)
          return
       end if

       lp2 = lp
       ok = eval_next(rdum(4),line,lp)
       if (ok) then
          x0(:,1) = rdum(1:3)
          x0(1,2) = rdum(4)
          ok = eval_next(x0(2,2),line,lp)
          ok = ok .and. eval_next(x0(3,2),line,lp)
          ok = ok .and. eval_next(x0(1,3),line,lp)
          ok = ok .and. eval_next(x0(2,3),line,lp)
          ok = ok .and. eval_next(x0(3,3),line,lp)
          if (.not.ok) then
             call ferror("struct_newcell","Wrong syntax in NEWCELL",faterr,line,syntax=.true.)
             return
          end if
       else
          lp = lp2
          do i = 1, 3
             x0(i,i) = rdum(i)
          end do
       end if

       doinv = .false.
       do while (.true.)
          word = lgetword(line,lp)
          if (equal(word,"origin")) then
             ok = eval_next(t0(1),line,lp)
             ok = ok .and. eval_next(t0(2),line,lp)
             ok = ok .and. eval_next(t0(3),line,lp)
             if (.not.ok) then
                call ferror("struct_newcell","Wrong ORIGIN syntax in NEWCELL",faterr,line,syntax=.true.)
                return
             end if
          elseif (equal(word,"inv").or.equal(word,"inverse")) then
             doinv = .true.
          elseif (len_trim(word) > 0) then
             call ferror('struct_newcell','Unknown extra keyword',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do
       if (doinv) call matinv(x0,3)

       ! transform to the new crystal
       call s%c%newcell(x0,t0)
       changed = .true.
    end if

    if (changed) then
       write (uout,'("  Lattice vectors of the new cell in the old setting (cryst. coord.):")')
       write (uout,'("    ",3(A," "))') (string(x0(i,1),'f',12,7,4),i=1,3)
       write (uout,'("    ",3(A," "))') (string(x0(i,2),'f',12,7,4),i=1,3)
       write (uout,'("    ",3(A," "))') (string(x0(i,3),'f',12,7,4),i=1,3)
       write (uout,'("  Origin translation: ",3(A," "))') (string(t0(i),'f',12,7,4),i=1,3)
       if (dotyp == 1 .and. dorefine) &
          write (uout,'("  Atomic positions have been refined.")')
       write (uout,*)

       ! reset all fields and properties to the default
       call s%reset_fields()

       ! report
       if (verbose) &
          call s%report(.true.,.true.,.true.,.true.,.true.,.true.,.false.)
    else
       write (uout,'("+ Cell transformation leads to the same cell: skipping."/)')
    end if

  end subroutine struct_newcell

  !> Load crystal or molecular vibrations from an external file
  module subroutine struct_vibrations(s,line,verbose)
    use tools_io, only: uout, getword, ferror, faterr
    use param, only: isformat_unknown
    type(system), intent(inout) :: s
    character*(*), intent(in) :: line
    logical, intent(in) :: verbose

    character(len=:), allocatable :: word, errmsg
    integer :: lp

    ! read the file name
    if (verbose) &
       write (uout,'("* VIBRATIONS: reading vibrational frequencies and modes")')
    lp=1
    word = getword(line,lp)
    if (verbose) &
       write (uout,'("  File: ",A)') trim(word)

    ! read the file
    call s%c%read_vibrations_file(word,isformat_unknown,errmsg)
    if (len_trim(errmsg) > 0) &
       call ferror("struct_vibrations",errmsg,faterr)
    write (uout,*)

  end subroutine struct_vibrations

  !> Try to determine the molecular cell from the crystal geometry
  module subroutine struct_molcell(s,line)
    use systemmod, only: system
    use global, only: rborder_def, eval_next, dunit0, iunit
    use tools_io, only: ferror, faterr, uout, string

    type(system), intent(inout) :: s
    character*(*), intent(in) :: line

    integer :: i, j, lp
    real*8 :: xmin(3), xmax(3), rborder, raux
    logical :: ok

    if (.not.s%c%ismolecule) then
       call ferror('struct_molcell','MOLCELL works with MOLECULE, not CRYSTAL.',faterr,syntax=.true.)
       return
    end if
    if (any(abs(s%c%bb-90d0) > 1d-5)) then
       call ferror('struct_molcell','MOLCELL only allowed for orthogonal cells.',faterr,syntax=.true.)
       return
    end if

    ! defaults
    rborder = rborder_def

    ! read optional border input
    lp = 1
    ok = eval_next(raux,line,lp)
    if (ok) rborder = raux / dunit0(iunit)

    ! find the encompassing cube
    xmin = 1d40
    xmax = -1d30
    do i = 1, s%c%ncel
       do j = 1, 3
          xmin(j) = min(xmin(j),s%c%at(i)%x(j))
          xmax(j) = max(xmax(j),s%c%at(i)%x(j))
       end do
    end do

    ! apply the border
    do j = 1, 3
       xmin(j) = max(xmin(j) - rborder / s%c%aa(j),0d0)
       xmax(j) = min(xmax(j) + rborder / s%c%aa(j),1d0)
       s%c%molborder(j) = min(xmin(j),1d0-xmax(j))
    end do

    ! some output
    write (uout,'("* MOLCELL: set up a molecular cell")')
    write (uout,'("+ Limit of the molecule within the cell (cryst coords):")')
    write (uout,'("  a axis: ",A," -> ",A)') trim(string(s%c%molborder(1),'f',10,4)), trim(string(1d0-s%c%molborder(1),'f',10,4))
    write (uout,'("  b axis: ",A," -> ",A)') trim(string(s%c%molborder(2),'f',10,4)), trim(string(1d0-s%c%molborder(2),'f',10,4))
    write (uout,'("  c axis: ",A," -> ",A)') trim(string(s%c%molborder(3),'f',10,4)), trim(string(1d0-s%c%molborder(3),'f',10,4))

  end subroutine struct_molcell

  !> Calculate the Effective Coordination Number (ECON).
  module subroutine struct_econ(s)
    use systemmod, only: system
    use crystalmod, only: crystal
    use global, only: iunitname0, dunit0, iunit
    use tools_io, only: uout, string, ioj_left, ioj_right, ferror, faterr
    use param, only: icrd_crys, bohrtoa
    type(system), intent(inout) :: s

    logical :: ok
    integer :: i, j, k, n
    integer :: nat
    integer, allocatable :: eid(:)
    real*8 :: dist0, econ, up2d
    real*8 :: wi, numer, econprev
    real*8, allocatable :: econij(:,:), ndij(:,:), dist(:)
    real*8, allocatable :: econij_noit(:,:), ndij_noit(:,:), mindist(:)
    character(len=:), allocatable :: str, namestr

    real*8, parameter :: up2dmax = 15d0 / bohrtoa
    real*8, parameter :: wthresh = 1d-10
    real*8, parameter :: epsfac = (1d0 + 14d0 * log(10d0))**(1d0/6d0)

    if (s%c%ismolecule) then
       call ferror('struct_econ','ECON only available for crystals',faterr,syntax=.true.)
       return
    end if

    allocate(econij(0:s%c%nspc,s%c%nneq),econij_noit(0:s%c%nspc,s%c%nneq))
    allocate(ndij(0:s%c%nspc,s%c%nneq),ndij_noit(0:s%c%nspc,s%c%nneq),mindist(0:s%c%nspc))
    econij = 0d0
    econij_noit = 0d0
    ndij = 0d0
    ndij_noit = 0d0

    ! initialize
    up2d = up2dmax

    main: do while(.true.)
       ! loop over non-equivalent atoms
       do i = 1, s%c%nneq
          ! skip if this is a single atom in a box
          if (s%c%ismolecule .AND. s%c%ncel == 1) cycle

          ! grab the environment for this atom
          call s%c%list_near_atoms(s%c%at(i)%x,icrd_crys,.true.,nat,eid,dist,up2d=up2d,nozero=.true.)

          ! get the minimum distance for all the species
          mindist = -1d0
          mindist(0) = dist(1)
          do k = 1, nat
             if (mindist(s%c%atcel(eid(k))%is) < 0d0) mindist(s%c%atcel(eid(k))%is) = dist(k)
             if (all(mindist >= 0d0)) exit
          end do

          ! loop over species
          do j = 0, s%c%nspc
             econprev = -1d0
             econ = 0d0
             dist0 = mindist(j)
             n = 0
             do while (abs(econ-econprev) > wthresh .or. n <= 2)
                n = n + 1
                econprev = econ

                ! loop over the environment atoms and consider only species j
                numer = 0d0
                econ = 0d0
                do k = 1, nat
                   if (j /= 0 .and. s%c%atcel(eid(k))%is /= j) cycle
                   wi = exp(1 - (dist(k)/dist0)**6)
                   numer = numer + dist(k) * wi
                   econ = econ + wi
                end do
                ok = (econ > 0d0)
                if (ok) dist0 = numer / econ

                ! if dist0 is dangerously close to the up2d, increase by a factor and try again
                if (.not.ok .or. dist0 > up2d / epsfac) then
                   up2d = 1.5d0 * up2d
                   cycle main
                end if

                ! save the non-iterative econ and average bond length values
                if (n == 1) then
                   ndij_noit(j,i) = dist0
                elseif (n == 2) then
                   econij_noit(j,i) = econ
                end if
             end do ! while (abs(econ-econprev) > wthresh)
             econij(j,i) = econ
             ndij(j,i) = dist0
          end do ! j = 0, s%c%nspc
       end do ! i = 1, s%c%nneq
       exit
    end do main

    ! write the output
    write (uout, '("* ECON")')
    write (uout, '("+ Effective Coordination Number (ECoN), per-species")')
    write (uout, '("# Please cite:")')
    write (uout, '("#   Nespolo, Acta Cryst. B, 72 (2016) 51. (doi:10.1107/S2052520615019472)")')
    write (uout, '("#   Hoppe, Z. Kristallogr. 150 (1979) 23. (doi:10.1524/zkri.1979.150.14.23)")')
    write (uout, '("# nid = non-equivalent atom ID. name = atomic name (symbol)")')
    write (uout, '("# spc = atomic species (* means all species)")')
    write (uout, '("# econ = effective coordination number species spc around non-equivalent")')
    write (uout, '("#        atom nid (iterative variant)")')
    write (uout, '("# 1econ = effective coordination number species spc around non-equivalent")')
    write (uout, '("#         atom nid (non-iterative variant using the shortest bond length))")')
    write (uout, '("# nd = mean weighted distance (",A,"), Eq. 3 in Nespolo, iterative variant. ")') iunitname0(iunit)
    write (uout, '("# 1nd = mean weighted distance (",A,"), non-iterative variant. ")') iunitname0(iunit)
    write (uout, '("# nid->spc   name(nid)->name(spc)     econ      1econ       nd         1nd")')
    do i = 1, s%c%nneq
       do j = 0, s%c%nspc
          if (j == 0) then
             str = "*   "
             namestr = "*      "
          else
             str = string(j,length=4,justify=ioj_left)
             namestr = string(s%c%spc(j)%name,length=7,justify=ioj_left)
          end if

          write (uout,'(A," -> ",A,A," -> ",99(A," "))') string(i,length=4,justify=ioj_right), str,&
             string(s%c%spc(s%c%at(i)%is)%name,length=9,justify=ioj_right), namestr,&
             string(econij(j,i),'f',length=10,decimal=4,justify=ioj_right),&
             string(econij_noit(j,i),'f',length=10,decimal=4,justify=ioj_right),&
             string(ndij(j,i)*dunit0(iunit),'f',length=10,decimal=4,justify=ioj_right),&
             string(ndij_noit(j,i)*dunit0(iunit),'f',length=10,decimal=4,justify=ioj_right)
       end do
    end do
    write (uout,*)

    deallocate(econij,ndij,econij_noit,ndij_noit,mindist)

  end subroutine struct_econ

  !> Identify atoms by their coordinates
  module subroutine struct_identify(s,line0,lp)
    use systemmod, only: system
    use global, only: iunit, iunit_bohr, iunit_ang, iunitname0, dunit0, &
       eval_next
    use tools_io, only: lgetword, getword, getline, equal, ferror,&
       faterr, uin, ucopy, uout, string, ioj_left, ioj_center, fopen_read,&
       lower
    use param, only: bohrtoa
    use types, only: realloc
    type(system), intent(in) :: s
    character*(*), intent(in) :: line0
    integer, intent(inout) :: lp

    logical :: ok
    character(len=:), allocatable :: line, word, lword
    real*8 :: x0(3), xmin(3), xmax(3), x0out(3)
    real*8, allocatable :: pointlist(:,:)
    integer :: i, j, n, idx
    logical :: found, doenv

    integer :: ldunit, unit, mm
    integer, parameter :: unit_au = 1
    integer, parameter :: unit_ang = 2
    integer, parameter :: unit_x = 3

    real*8, parameter :: eps = 1d-4

    ! default units
    if (s%c%ismolecule) then
       if (iunit == iunit_bohr) then
          ldunit = unit_au
       elseif (iunit == iunit_ang) then
          ldunit = unit_ang
       end if
    else
       ldunit = unit_x
    endif

    ! parse the first word
    doenv = .true.
    word = getword(line0,lp)
    lword = lower(word)
    if (equal(lword,'angstrom') .or.equal(lword,'ang')) then
       ldunit = unit_ang
    elseif (equal(lword,'bohr') .or.equal(lword,'au')) then
       ldunit = unit_au
    elseif (equal(lword,'cryst')) then
       ldunit = unit_x
    elseif (len_trim(lword) > 0) then
       doenv = .false.
    endif

    ! read the input coordinates
    allocate(pointlist(3,10))
    n = 0
    if (doenv) then
       word = getword(line0,lp)
       lword = lower(word)
       if (len_trim(lword) > 0) then
          call ferror('struct_identify','Unkwnon extra keyword',faterr,line0,syntax=.true.)
          return
       end if

       lp = 1
       ok = getline(uin,line,ucopy=ucopy)
       if (ok) then
          lword = lgetword(line,lp)
       else
          lword = ""
          lp = 1
       end if
       do while (ok.and..not.equal(lword,'endidentify').and..not.equal(lword,'end'))
          lp = 1
          ok = eval_next (x0(1), line, lp)
          ok = ok .and. eval_next (x0(2), line, lp)
          ok = ok .and. eval_next (x0(3), line, lp)
          if (ok) then
             ! this is a point, parse the last word
             lword = lgetword(line,lp)
             if (equal(lword,'angstrom') .or.equal(lword,'ang')) then
                unit = unit_ang
             elseif (equal(lword,'bohr') .or.equal(lword,'au')) then
                unit = unit_au
             elseif (equal(lword,'cryst')) then
                unit = unit_x
             else
                unit = ldunit
             endif

             if (unit == unit_ang) then
                x0 = s%c%c2x(x0 / bohrtoa - s%c%molx0)
             elseif (unit == unit_au) then
                x0 = s%c%c2x(x0 - s%c%molx0)
             endif
             n = n + 1
             if (n > size(pointlist,2)) then
                call realloc(pointlist,3,2*n)
             end if
             pointlist(:,n) = x0
          else
             ! this is an xyz file
             call readxyz()
          endif

          ! read next line
          lp = 1
          ok = getline(uin,line,ucopy=ucopy)
          if (ok) then
             lword = lgetword(line,lp)
          else
             line = ""
             lp = 1
          end if
       enddo
       lword = lgetword(line,lp)
       if (len_trim(lword) > 0) then
          call ferror('struct_identify','Unkwnon extra keyword',faterr,line,syntax=.true.)
          return
       end if
    else
       call readxyz()
    endif

    xmin = 1d40
    xmax = -1d40
    found = .false.
    ! identify the atoms
    write(uout,'("* IDENTIFY: match input coordinates to atoms or CPs in the structure")')
    if (.not.s%c%ismolecule) then
       write(uout,'("# (x,y,z) is the position in crystallographic coordinates ")')
    else
       write(uout,'("# (x,y,z) is the position in Cartesian coordinates (",A,")")') &
          iunitname0(iunit)
    end if
    write(uout,'("# id        x             y             z     mult name  ncp  cp")')

    do i = 1, n
       x0 = pointlist(:,i)
       x0out = x0
       if (s%c%ismolecule) x0out = (s%c%x2c(x0)+s%c%molx0) * dunit0(iunit)
       idx = s%f(s%iref)%identify_cp(x0,eps)
       mm = s%c%get_mult(x0)
       if (idx > 0) then
          write (uout,'(99(A," "))') string(i,length=4,justify=ioj_left), &
             (string(x0out(j),'f',length=13,decimal=8,justify=4),j=1,3), &
             string(mm,length=3,justify=ioj_center), &
             string(s%f(s%iref)%cpcel(idx)%name,length=5,justify=ioj_center), &
             string(s%f(s%iref)%cpcel(idx)%idx,length=4,justify=ioj_center), &
             string(idx,length=4,justify=ioj_center)
          do j = 1, 3
             xmin(j) = min(xmin(j),x0(j))
             xmax(j) = max(xmax(j),x0(j))
             found = .true.
          end do
       else
          write (uout,'(99(A," "))') string(i,length=4,justify=ioj_left), &
             (string(x0out(j),'f',length=13,decimal=8,justify=4),j=1,3), &
             string(mm,length=3,justify=ioj_center), &
             string(" --- not found --- ")
       endif
    end do
    deallocate(pointlist)

    if (found) then
       if (.not.s%c%ismolecule) then
          write(uout,'("#")')
          write(uout,'("+ Cube, x0 (cryst): ",3(A," "))') (string(xmin(j),'f',decimal=8),j=1,3)
          write(uout,'("+ Cube, x1 (cryst): ",3(A," "))') (string(xmax(j),'f',decimal=8),j=1,3)
          xmin = s%c%c2x(s%c%x2c(xmin) - 2)
          xmax = s%c%c2x(s%c%x2c(xmax) + 2)
          write(uout,'("+ Cube + 2 bohr, x0 (cryst): ",3(A," "))') (string(xmin(j),'f',decimal=8),j=1,3)
          write(uout,'("+ Cube + 2 bohr, x1 (cryst): ",3(A," "))') (string(xmax(j),'f',decimal=8),j=1,3)
          xmin = s%c%c2x(s%c%x2c(xmin) - 3)
          xmax = s%c%c2x(s%c%x2c(xmax) + 3)
          write(uout,'("+ Cube + 5 bohr, x0 (cryst): ",3(A," "))') (string(xmin(j),'f',decimal=8),j=1,3)
          write(uout,'("+ Cube + 5 bohr, x1 (cryst): ",3(A," "))') (string(xmax(j),'f',decimal=8),j=1,3)
       else
          xmin = s%c%x2c(xmin)+s%c%molx0
          xmax = s%c%x2c(xmax)+s%c%molx0
          write(uout,'("+ Cube, x0 (",A,"): ",3(A," "))') iunitname0(iunit), (string(xmin(j)*dunit0(iunit),'f',decimal=8),j=1,3)
          write(uout,'("+ Cube, x1 (",A,"): ",3(A," "))') iunitname0(iunit), (string(xmax(j)*dunit0(iunit),'f',decimal=8),j=1,3)
       end if
    end if
    write(uout,*)

  contains
    !> Read an xyz file and add the coordinates to pointlist
    subroutine readxyz()
      use types, only: realloc
      use tools_io, only: fclose

      integer :: lu, nat, i
      real*8 :: x0(3)
      character*20 :: atsym

      lu = fopen_read(word)
      read(lu,*) nat
      read(lu,*)
      do i = 1, nat
         read(lu,*) atsym, x0
         x0 = s%c%c2x(x0 / bohrtoa - s%c%molx0)
         n = n + 1
         if (n > size(pointlist,2)) then
            call realloc(pointlist,3,2*n)
         end if
         pointlist(:,n) = x0
      end do
      call fclose(lu)
    end subroutine readxyz

  end subroutine struct_identify

  !> Write a .mols file for neighcrys. Works even if Z' > 1 or if
  !> there are more than two molecules.
  module subroutine struct_makemols_neighcrys(line0,lp)
    use systemmod, only: system
    use tools_math, only: det3, cross
    use tools_io, only: getword, ferror, faterr, string, getline_raw, fopen_read, fclose, &
       fopen_write
    use types, only: realloc
    character*(*), intent(in) :: line0
    integer, intent(inout) :: lp

    character(len=:), allocatable :: errmsg
    character(len=:), allocatable :: file, line, line2
    logical :: laux, ok, isfs
    integer :: lu, i, j
    integer :: ncel, nmol, nmol2, idum, idax(3)
    integer, allocatable :: imap(:), idmol(:), iz(:)
    real*8, allocatable :: xx(:,:)
    real*8 :: rmat(3,3), xv(3), xnn
    character*10 :: cdum
    character*15, allocatable :: atlbl(:)

    real*8, parameter :: eps = 1d-6

    ! read the molecules in from the fort.21 file
    file = getword(line0,lp)
    inquire(file=file,exist=laux)
    if (.not.laux) &
       call ferror("makemols_neighcrys","file not found: " // string(file),faterr)

    ncel = 0
    nmol = 0
    lu = fopen_read(file)
    main: do while (getline_raw(lu,line))
       if (trim(line) == " Equivalent basis atoms") then
          ! default error message
          errmsg = "Unexpected error or end of line reading Equivalent basis atoms block"

          ! skip four lines
          ok = getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          ok = ok.and.getline_raw(lu,line)
          if (.not.ok) goto 999

          ! allocate space for the atomic information
          allocate(imap(10),idmol(10),iz(10),atlbl(10))

          ! read the atomic information and check against the loaded structure
          do while (.true.)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (len_trim(line) == 0) exit

             ncel = ncel + 1
             if (ncel > size(imap,1)) then
                call realloc(imap,2*ncel)
                call realloc(idmol,2*ncel)
                call realloc(iz,2*ncel)
                call realloc(atlbl,2*ncel)
             end if
             read(line,*,err=999,end=999) idum, iz(ncel), atlbl(ncel), cdum, idmol(ncel)
          end do
          if (ncel == 0) &
             call ferror("makemols_neighcrys","no atoms found in the fort.21 file",faterr)

          ! keep reading until we find the "crystallographic" keyword
          isfs = .false.
          do while (.true.)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (index(line,"foreshortened") > 0) isfs = .true.
             if (index(line,"crystallographic") > 0) exit
          end do

          ! reallocate and allocate coordinates
          call realloc(imap,ncel)
          call realloc(idmol,ncel)
          call realloc(iz,ncel)
          call realloc(atlbl,ncel)
          allocate(xx(3,ncel))

          ! read the coordinates
          errmsg = "Unexpected error or end of line reading the coordinates block"
          ncel = 0
          do while (.true.)
             ok = getline_raw(lu,line)
             if (.not.ok) goto 999
             if (len_trim(line) == 0) exit
             ok = ok .and. getline_raw(lu,line2)
             if (.not.ok) goto 999
             if (isfs) then
                ok = ok .and. getline_raw(lu,line2)
                if (.not.ok) goto 999
             end if

             ncel = ncel + 1
             if (ncel > size(imap,1)) then
                errmsg = "Too many atoms in coordinates block"
                goto 999
             end if
             read(line,*,err=999,end=999) idum, cdum, xx(:,ncel)
          end do

          ! we are done with this file
          nmol = maxval(idmol)
          exit main
       end if
    end do main
    call fclose(lu)
    errmsg = ""
    if (ncel == 0) &
       call ferror("makemols_neighcrys","no atoms found in the fort.21 file",faterr)
    if (nmol == 0) &
       call ferror("makemols_neighcrys","no molecules found in the fort.21 file",faterr)

    ! read the mols file
    file = getword(line0,lp)
    if (len_trim(file) == 0) then
       errmsg = "A name for the mols file is required"
       goto 999
    end if

    ! count the actual number of molecules
    nmol2 = 0
    do i = 1, nmol
       if (any(idmol(1:ncel) == i)) nmol2 = nmol2 + 1
    end do
    if (nmol2 == 0) &
       call ferror("makemols_neighcrys","no inequivalent molecules found in the fort.21 file",faterr)

    ! write the mols file
    lu = fopen_write(file)
    write (lu,'("MOLX ",A)') string(nmol2)
    do i = 1, nmol
       if (all(idmol(1:ncel) /= i)) cycle
       idax = 0
       do j = 1, ncel
          if (idmol(j) /= i) cycle

          if (idax(1) == 0) then
             idax(1) = j
             rmat(:,1) = xx(:,j)
          elseif (idax(2) == 0) then
             xnn = norm2(xx(:,j) - xx(:,idax(1)))
             if (xnn < eps) cycle
             idax(2) = j
             rmat(:,2) = (xx(:,j) - xx(:,idax(1))) / xnn
          else
             xnn = norm2(xx(:,j) - xx(:,idax(1)))
             if (xnn < eps) cycle
             rmat(:,3) = (xx(:,j) - xx(:,idax(1))) / xnn
             xv = cross(rmat(:,2),rmat(:,3))
             if (norm2(xv) >= eps) then
                idax(3) = j
                exit
             end if
          end if
       end do

       if (any(idax == 0)) then
          ! ! the molecule must be linear
          ! ! find any third atom in the system that is not collinear
          ! do j = 1, ncel
          !    if (idmol(j) == i) cycle
          !    xnn = norm2(xx(:,j) - xx(:,idax(1)))
          !    if (xnn < eps) cycle
          !    rmat(:,3) = (xx(:,j) - xx(:,idax(1))) / xnn
          !    xv = cross(rmat(:,2),rmat(:,3))
          !    if (norm2(xv) >= eps) then
          !       idax(3) = j
          !       exit
          !    end if
          ! end do
          call ferror("makemols_neighcrys","cannot be applied to a linear molecule",faterr)
       end if

       write (lu,'("X LINE  ",A," ",A," 1")') trim(atlbl(idax(1))), trim(atlbl(idax(2)))
       write (lu,'("Y PLANE ",A," ",A," 1"," ",A," 1")') &
          trim(atlbl(idax(1))), trim(atlbl(idax(2))), trim(atlbl(idax(3)))
    end do
    write (lu,'("ENDS")')
    call fclose(lu)

    return
999 continue

    call ferror('struct_makemols_neighcrys',errmsg,faterr,line,syntax=.true.)
    return

  end subroutine struct_makemols_neighcrys

  !> Put the atoms of one molecule in the same order as another.
  module subroutine struct_molreorder(line,lp)
    use systemmod, only: sy
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use tools, only: qcksort
    use tools_math, only: crosscorr_triangle, rmsd_walker, umeyama_graph_matching,&
       ullmann_graph_matching
    use tools_io, only: getword, string, ferror, faterr, lower, equal, uout, string
    use types, only: realloc
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    type(crystal) :: cref, cx
    type(crystal), allocatable :: c(:)
    type(crystalseed) :: seed
    character*1024 :: fname(2)
    character(len=:), allocatable :: word, lword, wfile, msg
    integer :: i, n, is, ns, nat
    real*8 :: rms1, rms2, rmsmin, q1(3,3), q2(3,3), xcm1(3), xcm2(3), xnew(3)
    integer, allocatable :: isuse(:), isperm(:,:)
    integer, allocatable :: cidxorig(:,:)
    real*8, allocatable :: x1(:,:), x2(:,:), rmsd(:)
    logical, allocatable :: isinv(:)
    real*8, allocatable :: dref(:,:), ddg(:,:), ddh(:,:), mrot(:,:,:)
    logical :: moveatoms, doinv, umeyama
    integer, allocatable :: iz1(:), iz2(:), ncon1(:), ncon2(:), idcon1(:,:), idcon2(:,:), list(:,:)
    integer :: mcon, nlist

    integer, parameter :: isuse_valid = 0
    integer, parameter :: isuse_different_nat = 1
    integer, parameter :: isuse_incompatible_z = 2
    integer, parameter :: isuse_permutation_not_found = 3

    write(uout,'("* MOLREORDER: reorder the atoms in a molecule or a molecular crystal")')

    ! read the input
    ns = 0
    wfile = ""
    moveatoms = .false.
    doinv = .false.
    umeyama = .false.
    do while(.true.)
       word = getword(line,lp)
       lword = lower(word)
       if (equal(lword,'write')) then
          wfile = getword(line,lp)
       elseif (equal(lword,'moveatoms')) then
          moveatoms = .true.
       elseif (equal(lword,'inv')) then
          doinv = .true.
       elseif (equal(lword,'umeyama')) then
          umeyama = .true.
       elseif (equal(lword,'ullmann')) then
          umeyama = .false.
       elseif (len_trim(word) == 0) then
          exit
       else
          ns = ns + 1
          if (ns > 2) &
             call ferror("struct_molreorder","cannot order more than two structures",faterr)
          fname(ns) = word
       end if
    end do
    if (ns /= 2) &
       call ferror("struct_molreorder","need two structures to order",faterr)

    ! get the reference structure and check it is a molecule
    if (trim(fname(1)) == ".") then
       cref = sy%c
    else
       call struct_crystal_input(fname(1),-1,.false.,.false.,cr0=cref)
    end if
    if (.not.cref%isinit) &
       call ferror("struct_molreorder","could not load structure" // string(fname(1)),faterr)
    if (.not.cref%ismolecule) &
       call ferror("struct_molreorder","the first structure is not a molecule",faterr)
    if (cref%nmol > 1) &
       call ferror("struct_molreorder","the first structure contains more than one molecule",faterr)

    ! get the other structure (do not guess the symmetry)
    if (trim(fname(2)) == ".") then
       cx = sy%c
    else
       call struct_crystal_input(fname(2),-1,.false.,.false.,cr0=cx)
    end if
    if (.not.cx%isinit) &
       call ferror("struct_molreorder","could not load structure" // string(fname(2)),faterr)
    if (.not.cx%ismolecule) then
       if (.not.cx%ismol3d) &
          call ferror("struct_molreorder","the target structure is not a molecular crystal",faterr)
       call cx%wholemols()
       if (any(cx%idxmol(1:cx%nmol) < 0)) &
          call ferror("struct_molreorder","the target structure must have whole molecules",faterr)
    end if

    ! read the fragments from cx into structures
    ns = cx%nmol
    allocate(c(ns))
    do i = 1, ns
       call seed%from_fragment(cx%mol(i),.true.)
       seed%border = 0d0
       call c(i)%struct_new(seed,.true.)
    end do

    ! allocate the use and permutation arrays
    allocate(isuse(ns),isperm(cref%ncel,ns))
    isuse = isuse_valid

    ! check that the number of atoms are equal
    do i = 1, ns
       if (cref%ncel /= c(i)%ncel) isuse(i) = isuse_different_nat
    end do
    nat = cref%ncel

    ! allocate some work space
    allocate(rmsd(ns),x1(3,nat),x2(3,nat),isinv(ns))

    ! make the mapping between structure+atom and the original cidx
    ! sort because the from_fragment routine also sorts
    allocate(cidxorig(nat,ns))
    do is = 1, ns
       if (isuse(is) /= isuse_valid) cycle
       do i = 1, nat
          cidxorig(i,is) = cx%mol(is)%at(i)%cidx
       end do
       call qcksort(cidxorig(:,is))
    end do

    ! calculate the permutations (isperm) with the appropriate method.
    ! if the permutation is not found, deactivate the structure.
    if (umeyama) then
       ! use the umeyama method
       allocate(dref(nat,nat),ddg(nat,nat),ddh(nat,nat))
       call cref%distmatrix(dref,conn=.true.)
       do is = 1, ns
          if (isuse(is) /= isuse_valid) cycle
          ddg = dref
          call c(is)%distmatrix(ddh,conn=.true.)
          call umeyama_graph_matching(nat,ddg,ddh,isperm(:,is))
          if (any(cref%spc(cref%at(1:nat)%is)%z /= c(is)%spc(c(is)%at(isperm(1:nat,is))%is)%z)) &
             isuse(is) = isuse_incompatible_z
       end do
       deallocate(dref,ddg,ddh)
    else
       ! use the ullmann method
       ! obtain the maximum number of bonds per atom
       mcon = 0
       do i = 1, cref%ncel
          mcon = max(mcon,cref%nstar(i)%ncon)
       end do
       do is = 1, ns
          if (isuse(is) /= isuse_valid) cycle
          do i = 1, c(is)%ncel
             mcon = max(mcon,c(is)%nstar(i)%ncon)
          end do
       end do

       ! allocate work space and prepare reference arrays
       allocate(iz1(nat),ncon1(nat),idcon1(mcon,nat),iz2(nat),ncon2(nat),idcon2(mcon,nat))
       idcon1 = 0
       do i = 1, nat
          iz1(i) = cref%spc(cref%atcel(i)%is)%z
          ncon1(i) = cref%nstar(i)%ncon
          idcon1(1:ncon1(i),i) = cref%nstar(i)%idcon(1:ncon1(i))
       end do

       ! run over all target structures
       do is = 1, ns
          if (isuse(is) /= isuse_valid) cycle

          ! arrays for the target structure
          idcon2 = 0
          do i = 1, nat
             iz2(i) = c(is)%spc(c(is)%atcel(i)%is)%z
             ncon2(i) = c(is)%nstar(i)%ncon
             idcon2(1:ncon2(i),i) = c(is)%nstar(i)%idcon(1:ncon2(i))
          end do

          ! run ullmann to get the list of candidates
          call ullmann_graph_matching(iz1,ncon1,idcon1,iz2,ncon2,idcon2,nlist,list)
          if (nlist == 0) then
             isuse(is) = isuse_permutation_not_found
          else
             ! select the candidate with the lowest rms
             rmsmin = huge(1d0)
             do i = 1, nlist
                call calculate_rms(is,list(:,i),.false.,rms1,q1)
                if (rms1 < rmsmin) then
                   isperm(:,is) = list(:,i)
                   rmsmin = rms1
                end if
                if (doinv) then
                   call calculate_rms(is,list(:,i),.true.,rms1,q1)
                   if (rms1 < rmsmin) then
                      isperm(:,is) = list(:,i)
                      rmsmin = rms1
                   end if
                end if
             end do
          end if
       end do

       ! wrap up
       deallocate(iz1,ncon1,idcon1,iz2,ncon2,idcon2)
    end if

    ! Get the rmsd from walker. If moveatoms, also save the rotation.
    if (moveatoms) allocate(mrot(3,3,ns))
    rmsd = huge(1d0)
    isinv = .false.
    do is = 1, ns
       if (isuse(is) /= isuse_valid) cycle

       call calculate_rms(is,isperm(:,is),.false.,rms1,q1)
       isinv(is) = .false.
       if (doinv) then
          call calculate_rms(is,isperm(:,is),.true.,rms2,q2)
          if (rms2 < rms1) then
             isinv(is) = .true.
             rms1 = rms2
             q1 = q2
          end if
       end if
       rmsd(is) = rms1
       if (moveatoms) mrot(:,:,is) = q1
    end do
    deallocate(x1,x2)

    ! report on the results
    write (uout,'("+ Reference structure: ",A)') string(cref%file)
    write (uout,'("+ Target structure: ",A)') string(fname(2))
    if (cx%ismolecule) then
       write (uout,'("  The target structure is a molecule with ",A," fragments.")') string(ns)
    else
       write (uout,'("  The target structure is a molecular crystal with ",A," fragments.")') string(ns)
    end if
    if (umeyama) then
       write (uout,'("+ Method: umeyama")')
    else
       write (uout,'("+ Method: ullmann")')
    end if
    write (uout,'("# List of reordered structures (",A,")")') string(ns)
    write (uout,'("#id -- Reorder result --")')
    do is = 1, ns
       if (isuse(is) == isuse_different_nat) then
          msg = "the target and reference structures have different number of atoms"
       elseif (isuse(is) == isuse_incompatible_z) then
          msg = "the target and reference structures have different atomic types"
       elseif (isuse(is) == isuse_permutation_not_found) then
          msg = "ullmann could not find an appropriate permutation for the atomic sequence"
       else
          msg = ""
          do i = 1, nat
             msg = msg // " " // string(isperm(i,is))
          end do
       end if
       if (isuse(is) /= isuse_valid) then
          write (uout,'(" ",A," skipped (",A,")")') string(is,2), string(msg)
       else
          write (uout,'(" ",A," success (rmse=",A,",inv=",A,",perm=",A,")")') &
             string(is,2), trim(string(rmsd(is),'e',decimal=4)), string(isinv(is)), string(msg)
       end if
    end do
    write (uout,*)
    deallocate(rmsd,isinv)

    ! write the final structure
    if (len_trim(wfile) > 0 .and. all(isuse == isuse_valid)) then
       write (uout,'("* WRITE file: ",A/)') string(wfile)

       xcm1 = cref%mol(1)%cmass(.false.)

       call cx%makeseed(seed,.false.)
       n = 0
       do is = 1, ns
          xcm2 = cx%mol(is)%cmass(.false.)
          do i = 1, nat
             n = n + 1
             if (moveatoms) then
                xnew = cref%mol(1)%at(i)%r - xcm1
                xnew = matmul(mrot(:,:,is),xnew)
                xnew = xnew + xcm2
                xnew = cx%c2x(xnew)
             else
                xnew = cx%atcel(cidxorig(isperm(i,is),is))%x
             end if
             seed%x(:,n) = xnew
             seed%is(n) = cx%atcel(cidxorig(isperm(i,is),is))%is
          end do
       end do

       call cx%struct_new(seed,.true.)
       call cx%wholemols()
       call cx%write_simple_driver(wfile)
    end if
    deallocate(isperm,cidxorig,isuse,c)

  contains

    ! Calculate the rms between the reference molecule and the target
    ! molecule is using permutation isperm. If inv, invert the
    ! coordinates of the target molecule. Returns the rmsd and,
    ! possibly, the rotation matrix. Requires nat and having x1 and x2
    ! allocated.
    subroutine calculate_rms(is,isperm,inv,rms,mrot)
      integer, intent(in) :: is
      integer, intent(in) :: isperm(nat)
      logical, intent(in) :: inv
      real*8, intent(out) :: rms
      real*8, optional, intent(out) :: mrot(3,3)

      integer :: i
      real*8 :: q(3,3)

      if (present(mrot)) mrot = 0d0
      rms = huge(1d0)
      if (isuse(is) /= isuse_valid) return

      do i = 1, nat
         x1(:,i) = cref%at(i)%r
         x2(:,i) = c(is)%at(isperm(i))%r
      end do
      if (inv) x2 = -x2
      rms = rmsd_walker(x1,x2,q)
      if (present(mrot)) mrot = q

    end subroutine calculate_rms

  end subroutine struct_molreorder

  !> Move the atoms in a molecular crystal to match the atomic
  !> positions of the given molecules.
  module subroutine struct_molmove(line,lp)
    use crystalmod, only: crystal
    use crystalseedmod, only: crystalseed
    use fragmentmod, only: fragment
    use tools, only: qcksort
    use tools_math, only: crosscorr_triangle, rmsd_walker, umeyama_graph_matching
    use tools_io, only: getword, string, ferror, faterr, lower, equal, uout, string
    use types, only: realloc
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    integer :: i, j, n, ns, nat, idmin
    type(crystal) :: cini, caux, cfin
    character*1024, allocatable :: fname(:)
    character(len=:), allocatable :: word
    type(fragment), allocatable :: mol(:)
    type(crystalseed) :: seed
    real*8 :: rms, q1(3,3), xcm1(3), xcm2(3), xnew(3), dmin, dd
    real*8, allocatable :: x1(:,:), x2(:,:)
    integer, allocatable :: idx(:)
    logical, allocatable :: used(:)

    write(uout,'("* MOLMOVE: move the atoms in a molecule or a molecular crystal")')

    ! read the input
    ns = 0
    allocate(fname(10))
    do while (.true.)
       word = getword(line,lp)
       if (len_trim(word) == 0) exit

       ns = ns + 1
       if (ns > size(fname,1)) call realloc(fname,2*ns)
       fname(ns) = word
    end do
    if (ns < 3) &
       call ferror("struct_molmove","need more than two structures",faterr)

    ! write some output
    write (uout,'("+ Target structure: ",A)') string(fname(ns-1))
    do i = 1, ns - 2
    end do

    ! read the target
    call struct_crystal_input(fname(ns-1),-1,.false.,.false.,cr0=cini)
    if (cini%ismolecule.or..not.cini%ismol3d) &
       call ferror("struct_molmove","MOLMOVE can only be applied to molecular crystals (missing files?)",faterr)
    if (cini%nmol /= ns-2) &
       call ferror("struct_molmove","wrong number of molecules",faterr)

    ! read the fragments
    allocate(mol(ns-2))
    do i = 1, ns-2
       call struct_crystal_input(fname(i),1,.false.,.false.,cr0=caux)
       if (.not.caux%ismolecule.or.caux%nmol /= 1) &
          call ferror("struct_molmove","all fragments must be single molecules",faterr)
       mol(i) = caux%mol(1)
       do j = 1, caux%mol(1)%nat
          mol(i)%at(j)%r = mol(i)%at(j)%r + caux%molx0
       end do
    end do

    ! identify the molecules
    allocate(idx(cini%nmol),used(cini%nmol))
    used = .false.
    do i = 1, cini%nmol
       xcm1 = cini%mol(i)%cmass(.false.)

       dmin = 1d40
       do j = 1, cini%nmol
          if (used(j)) cycle
          xcm2 = mol(j)%cmass(.false.)
          dd = dot_product(xcm1-xcm2,xcm1-xcm2)
          if (dd < dmin) then
             dmin = dd
             idmin = j
          end if
       end do
       idx(i) = idmin
       used(idmin) = .true.
    end do
    deallocate(used)

    ! build the final structure
    call cini%makeseed(seed,.false.)
    n = 0
    do i = 1, cini%nmol
       ! check
       if (cini%mol(i)%nat /= mol(idx(i))%nat) &
          call ferror("struct_molmove","for molecule "//string(i)// " the number of atoms do not match",faterr)

       ! calculate the rotation matrix
       nat = cini%mol(i)%nat
       allocate(x1(3,nat),x2(3,nat))
       do j = 1, nat
          if (cini%spc(cini%mol(i)%at(j)%is)%z /= mol(idx(i))%spc(mol(idx(i))%at(j)%is)%z) &
             call ferror("struct_molmove","for molecule "//string(i)//" atom number "//&
             string(j)//" does not match",faterr)

          x1(:,j) = cini%mol(i)%at(j)%r
          x2(:,j) = mol(idx(i))%at(j)%r
       end do
       rms = rmsd_walker(x1,x2,q1)
       deallocate(x1,x2)

       ! some output
       xcm1 = mol(idx(i))%cmass(.false.)
       xcm2 = cini%mol(i)%cmass(.false.)
       write (uout,'("  Molecule ",A,": ",A," with rms = ",A," bohr, d(cm) = ",A," bohr")') string(i), &
          string(fname(idx(i))), string(rms,'f',decimal=4), string(norm2(xcm1-xcm2),'f',decimal=4)

       ! move the atoms and put them in the seed
       do j = 1, nat

          xnew = mol(idx(i))%at(j)%r - xcm1
          xnew = matmul(q1,xnew)
          xnew = xnew + xcm2
          xnew = cini%c2x(xnew)

          n = n + 1
          seed%x(:,n) = xnew
          seed%is(n) = cini%mol(i)%at(j)%is
       end do
    end do
    deallocate(mol,idx)

    ! write the final structure
    write (uout,'("+ Final structure: ",A)') string(fname(ns))
    write (uout,*)
    call cfin%struct_new(seed,.true.)
    call cfin%write_simple_driver(fname(ns))
    deallocate(fname)

  end subroutine struct_molmove

  !> Calculate k-point grid to a certain Rk using VASP's recipe.
  module subroutine struct_kpoints(s,line)
    use tools_io, only: lgetword, equal, isreal, ferror, faterr, string, uout,&
       ioj_right
    type(system), intent(in) :: s
    character*(*), intent(in) :: line

    integer :: lp, nk(3), nkold(3), i
    real*8 :: rk, rkold, rkmax
    character(len=:), allocatable :: word
    logical :: ok

    real*8 :: rkmax_def = 100d0
    real*8 :: rkmax_step = 0.1d0

    rk = -1d0
    rkmax = rkmax_def
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (len_trim(word) == 0) then
          exit
       elseif (equal(word,"rkmax")) then
          ok = isreal(rkmax,line,lp)
          if (.not.ok) then
             call ferror('struct_kpoints','Syntax error in KPOINTS/RKMAX',faterr,syntax=.true.)
             return
          end if
       elseif (isreal(rk,word)) then
          exit
       else
          call ferror('struct_kpoints','Unknown keyword in KPOINTS: ' // word,faterr,syntax=.true.)
          return
       endif
    end do

    write (uout,'("* KPOINTS: calculate dimensions of uniform k-point grid.")')
    write (uout,*)
    if (rk > 0d0) then
       call s%c%get_kpoints(rk,nk)
       write (uout,'("+ Rk = ",A," | kpts = ",3(A," "))') string(rk,'f',decimal=1),&
          (string(nk(i)),i=1,3)
    elseif (rkmax > 0d0) then
       write (uout,'("# ---  Rk  ---   -- kpts --")')
       rkold = rkmax_step
       call s%c%get_kpoints(rkold,nkold)
       do while (rk < rkmax)
          rk = rk + rkmax_step
          call s%c%get_kpoints(rk,nk)
          if (any(nk /= nkold)) then
             write (uout,'(A," -> ",A,"    ",3(A," "))') string(rkold,'f',5,1,ioj_right), &
                string(rk-rkmax_step,'f',5,1,ioj_right), (string(nkold(i),2),i=1,3)
             nkold = nk
             rkold = rk
          end if
       end do
       write (uout,'(A," -> ",A,"    ",3(A," "))') string(rkold,'f',5,1,ioj_right), &
          string(rk-rkmax_step,'f',5,1,ioj_right), (string(nkold(i),2),i=1,3)
    end if
    write (uout,*)

  end subroutine struct_kpoints

  !> Calculate and print information about the Brilloun zone.
  module subroutine struct_bz(s)
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc
    use global, only: iunitname0, iunit, dunit0
    use tools_math, only: matinv
    use tools_io, only: uout, faterr, ferror, string, ioj_right, ioj_center
    use tools, only: delaunay_reduction
    use param, only: pi
    type(system), intent(in) :: s

    interface
       ! The definitions and documentation for these functions are in doqhull.c
       subroutine runqhull_voronoi_step1(n,xstar,nf,nv,mnfv,fid) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
         integer(c_int), value :: n
         real(c_double) :: xstar(3,n)
         integer(c_int) :: nf, nv, mnfv
         type(c_ptr) :: fid
       end subroutine runqhull_voronoi_step1
       subroutine runqhull_voronoi_step2(nf,nv,mnfv,ivws,xvws,nfvws,fvws,fid) bind(c)
         use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
         integer(c_int), value :: nf, nv, mnfv
         integer(c_int) :: ivws(nf)
         real(c_double) :: xvws(3,nv)
         integer(c_int) :: nfvws(mnfv)
         integer(c_int) :: fvws(mnfv)
         type(c_ptr), value :: fid
       end subroutine runqhull_voronoi_step2
    end interface

    real*8 :: m_x2c_r(3,3), m_c2x_r(3,3), g(3,3), aa(3), bb(3), rmat(3,4), rdel(3,3)
    real*8 :: xstar(3,14), bz_ineighc(3,14), x0(3)
    integer :: i, j, n
    integer(c_int), allocatable :: ivws(:)
    real(c_double), allocatable :: xvws(:,:)
    integer :: bz_nf, bz_nv, bz_mnfv, bz_nside(14), bz_ineighx(3,14)
    integer, allocatable :: bz_iside(:,:)
    real*8, allocatable :: bz_x(:,:)
    type(c_ptr), target :: fid

    real*8, parameter :: eps_dnorm = 1d-5 !< minimum lattice vector length
    character*2, parameter :: lvecname(3) = (/"a*","b*","c*"/)

    ! header
    write (uout,'("* BZ: calculate Brillouin-Zone geometry.")')

    ! find the reduced cell in reciprocal space
    m_x2c_r = transpose(s%c%m_x2c)
    call matinv(m_x2c_r,3)
    m_c2x_r = transpose(s%c%m_x2c)
    call delaunay_reduction(m_x2c_r,rmat,rbas=rdel)

    ! cell lengths and angles
    g = matmul(transpose(m_x2c_r),m_x2c_r)
    do i = 1, 3
       aa(i) = sqrt(g(i,i))
    end do
    bb(1) = acos(g(2,3) / aa(2) / aa(3)) * 180d0 / pi
    bb(2) = acos(g(1,3) / aa(1) / aa(3)) * 180d0 / pi
    bb(3) = acos(g(1,2) / aa(1) / aa(2)) * 180d0 / pi

    ! construct star of lattice vectors -> use Delaunay reduction
    ! see 9.1.8 in ITC.
    n = 14
    xstar(:,1)  = rmat(:,1)
    xstar(:,2)  = rmat(:,2)
    xstar(:,3)  = rmat(:,3)
    xstar(:,4)  = rmat(:,4)
    xstar(:,5)  = rmat(:,1)+rmat(:,2)
    xstar(:,6)  = rmat(:,1)+rmat(:,3)
    xstar(:,7)  = rmat(:,2)+rmat(:,3)
    xstar(:,8)  = -(rmat(:,1))
    xstar(:,9)  = -(rmat(:,2))
    xstar(:,10) = -(rmat(:,3))
    xstar(:,11) = -(rmat(:,4))
    xstar(:,12) = -(rmat(:,1)+rmat(:,2))
    xstar(:,13) = -(rmat(:,1)+rmat(:,3))
    xstar(:,14) = -(rmat(:,2)+rmat(:,3))
    do i = 1, 14
       xstar(:,i) = matmul(m_x2c_r,xstar(:,i))
       if (norm2(xstar(:,i)) < eps_dnorm) &
          call ferror("wigner","Lattice vector too short. Please, check the unit cell definition.",faterr)
    end do

    call runqhull_voronoi_step1(n,xstar,bz_nf,bz_nv,bz_mnfv,fid)
    allocate(ivws(bz_nf),bz_iside(bz_mnfv,bz_nf),xvws(3,bz_nv))
    bz_nside = 0
    call runqhull_voronoi_step2(bz_nf,bz_nv,bz_mnfv,ivws,xvws,bz_nside(1:bz_nf),bz_iside,fid)

    if (allocated(bz_x)) deallocate(bz_x)
    allocate(bz_x(3,bz_nv))
    do i = 1, bz_nv
       bz_x(:,i) = matmul(m_c2x_r,xvws(:,i))
    end do

    bz_ineighc = 0
    bz_ineighx = 0
    do i = 1, bz_nf
       bz_ineighc(:,i) = xstar(:,ivws(i))
       bz_ineighx(:,i) = nint(matmul(m_c2x_r,xstar(:,ivws(i))))
    end do
    deallocate(ivws,xvws)

    write (uout,'("+ Reciprocal lattice vectors (",A,"^-1)")') iunitname0(iunit)
    do i = 1, 3
       write (uout,'("    ",A,": ",3(A," "))') lvecname(i),&
          (string(m_x2c_r(j,i)/dunit0(iunit),'f',length=16,decimal=10,justify=5),j=1,3)
    end do
    write (uout,'("  Reciprocal lattice parameters (bohr^-1): ",3(A,"  "))') &
       string(aa(1),'f',decimal=6), string(aa(2),'f',decimal=6), string(aa(3),'f',decimal=6)
    write (uout,'("  Reciprocal lattice angles (degrees): ",3(A,"  "))') &
       string(bb(1),'f',decimal=3), string(bb(2),'f',decimal=3), string(bb(3),'f',decimal=3)
    write (uout,*)

    write (uout,'("+ Vertex of the BZ in reciprocal cryst. coords. (",A,")")') string(bz_nv)
    write (uout,'("# id = vertex ID. xyz = vertex cryst. coords. d = vertex distance to origin (",A,"^-1).")') iunitname0(iunit)
    write (uout,'(5("  ",A))') string("id",length=3,justify=ioj_right),&
       string("x",length=11,justify=ioj_center),&
       string("y",length=11,justify=ioj_center),&
       string("z",length=11,justify=ioj_center),&
       string("d ("//iunitname0(iunit)//"^-1)",length=14,justify=ioj_center)
    do i = 1, bz_nv
       x0 = matmul(m_x2c_r,bz_x(:,i))
       write (uout,'(5("  ",A))') string(i,length=3,justify=ioj_right), &
          (string(bz_x(j,i),'f',length=11,decimal=6,justify=4),j=1,3), &
          string(norm2(x0)/dunit0(iunit),'f',length=14,decimal=8,justify=4)
    enddo
    write (uout,*)

    write (uout,'("+ Faces of the BZ cell (",A,")")') string(bz_nf)
    write (uout,'("# Face ID: vertexID1 vertexID2 ...")')
    do i = 1, bz_nf
       write (uout,'("  ",A,": ",999(A," "))') string(i,length=2,justify=ioj_right), &
          (string(bz_iside(j,i),length=2),j=1,bz_nside(i))
    end do
    write (uout,*)

    write (uout,'("+ Lattice vectors for the BZ neighbors")')
    write (uout,'("# FaceID: Voronoi lattice vector (cryst. coords.)")')
    do i = 1, bz_nf
       write (uout,'("  ",A,": ",99(A," "))') string(i,length=2,justify=ioj_right), &
          (string(bz_ineighx(j,i),length=2,justify=ioj_right),j=1,3)
    end do
    write (uout,*)

  end subroutine struct_bz

  !> Tools for editing and processing 1D x-ray powder diffraction
  !> patterns.
  module subroutine struct_xrpd(line0)
    use crystalmod, only: david_sivia_calculate_background, xrpd_peaklist
    use global, only: fileroot
    use tools_io, only: uout, getword, lgetword, equal, ferror, faterr, isinteger, string,&
       file_read_xy, fopen_write, fclose, isreal
    use tools, only: mergesort, qcksort
    use types, only: realloc
    character*(*), intent(in) :: line0

    character(len=:), allocatable :: word, xyfile, file, errmsg
    integer :: lp, nknot, n, lu, i, nadj
    logical :: ok
    real*8, allocatable :: x(:), y(:), yout(:), yfit(:)
    real*8 :: ymax_peakdetect, rms
    type(xrpd_peaklist) :: p

    integer, parameter :: nknot_def = 20

    ! header
    write (uout,'("* XRPD: processing of X-ray powder diffraction data")')

    lp = 1
    word = lgetword(line0,lp)
    if (equal(word,"background")) then
       ! background fit using David-Sivia's method
       !   BACKGROUND file-xy.s file-newxy.s [nknot.i]
       write (uout,'("+ BACKGROUND: Fitting background to experimental XRPD data")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   David and Sivia, J. Appl. Crystallogr. 34 (2001) 318-324 (doi:10.1107/S0021889801004332)")')

       ! read file names and header
       xyfile = getword(line0,lp)
       write (uout,'("  Reading the pattern from file: ",A)') trim(xyfile)
       file = getword(line0,lp)
       write (uout,'("  Writting background to file: ",A)') trim(file)
       ok = isinteger(nknot,line0,lp)
       if (.not.ok) nknot = nknot_def
       write (uout,'("  Number of knots: ",A)') string(nknot)

       ! read the pattern
       call file_read_xy(xyfile,n,x,y,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_xrpd',errmsg,faterr)

       ! calculate the background
       yout = david_sivia_calculate_background(n,x,y,errmsg,nknot=nknot)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_xrpd',errmsg,faterr)

       ! write the background to the file
       lu = fopen_write(file)
       write (lu,'("## x  yobs-ybackground  ybackground  yobs")')
       do i = 1, n
          write (lu,'(4(A," "))') string(x(i),'f',decimal=10),&
             string(y(i)-yout(i),'e',decimal=10), string(yout(i),'e',decimal=10),&
             string(y(i),'e',decimal=10)
       end do
       call fclose(lu)

    elseif (equal(word,"fit")) then
       ! fit an experimental XRPD profile
       !   FIT file-xy.s [ymax_peakdetect.r] [nadj.i]

       ! header
       write (uout,'("+ FIT: Fitting a peak profile to experimental XRPD data")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   Otero-de-la-Roza, J. Appl. Cryst. 57 (2024) 1401-1414 (doi:10.1107/S1600576724007489)")')

       ! read the input options
       file = getword(line0,lp)
       ymax_peakdetect = -huge(1d0)
       ok = isreal(ymax_peakdetect,line0,lp)
       if (ok) then
          nadj = 2
          ok = isinteger(nadj,line0,lp)
          if (nadj /= 1 .and. nadj /= 2) &
             call ferror('struct_xrpd','Number of adjacent points (NADJ) must be 1 or 2',faterr)
       end if

       ! run the peak fit
       write (uout,'("  Reading the pattern from file: ",A)') trim(file)
       call p%from_profile_file(file,rms,errmsg,.true.,ymax_peakdetect,nadj,xorig=x,yorig=y,ycalc=yfit)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_xrpd',errmsg,faterr)

       ! write the profile to disk
       file = fileroot // "_fit.dat"
       write (uout,'("  Writing original and fitted profile to file: ",A)') trim(file)
       lu = fopen_write(file)
       write (lu,'("## x yorig ycalc")')
       do i = 1, size(x,1)
          write (lu,'(3(A," "))') string(x(i),'f',decimal=10), string(y(i),'f',decimal=10),&
             string(yfit(i),'f',decimal=10)
       end do
       call fclose(lu)

       ! final list of peaks to disk
       file = fileroot // ".peaks"
       write (uout,'("+ List of peaks written to file: ",A)') file
       call p%write(file)
    elseif (equal(word,"refit")) then
       ! refit an experimental XRPD profile with a user-provided peaks file
       !   REFIT file-xy.s file.peaks

       ! read file names and header
       write (uout,'("+ REFIT: refit a powder pattern with a user-provided peaks file")')
       write (uout,'("# Please cite:")')
       write (uout,'("#   Otero-de-la-Roza, J. Appl. Cryst. 57 (2024) 1401-1414 (doi:10.1107/S1600576724007489)")')
       xyfile = getword(line0,lp)
       write (uout,'("  Reading the pattern from file: ",A)') trim(xyfile)
       file = getword(line0,lp)
       write (uout,'("  Reading the peaks from file: ",A)') trim(file)

       ! read the peaks file
       call p%from_peaks_file(file,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_xrpd',errmsg,faterr)

       ! read the pattern
       call file_read_xy(xyfile,n,x,y,errmsg)
       if (len_trim(errmsg) > 0) &
          call ferror('struct_xrpd',errmsg,faterr)

       ! fit the pattern
       call p%from_profile_file(xyfile,rms,errmsg,pkinput=p,xorig=x,&
          yorig=y,ycalc=yfit)

       ! write the profile to disk
       file = fileroot // "_fit.dat"
       write (uout,'("  Writing original and fitted profile to file: ",A)') trim(file)
       lu = fopen_write(file)
       write (lu,'("## x yorig ycalc")')
       do i = 1, size(x,1)
          write (lu,'(3(A," "))') string(x(i),'f',decimal=10), string(y(i),'f',decimal=10),&
             string(yfit(i),'f',decimal=10)
       end do
       call fclose(lu)

       ! final list of peaks to disk
       file = fileroot // ".peaks"
       write (uout,'("+ List of peaks written to file: fit.peaks")')
       call p%write(file)
    else
       call ferror("struct_xrpd","unknown keyword in XRPD: " // word,faterr)
    end if
    write (uout,*)

  end subroutine struct_xrpd

end submodule proc
