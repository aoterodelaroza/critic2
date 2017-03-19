! Copyright (c) 2015 Alberto Otero de la Roza
! <aoterodelaroza@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>.
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.

! Geometry of the crystal: variables and tools.
module struct
  implicit none

  private
  public :: struct_crystal_input
  public :: struct_clearsym
  public :: struct_charges
  public :: struct_write
  public :: struct_powder
  public :: struct_rdf
  public :: struct_compare
  public :: struct_environ
  public :: struct_packing
  public :: struct_newcell
  public :: struct_molcell

contains

  !xx! top-level routines
  !> Parse the input of the crystal keyword
  subroutine struct_crystal_input(c,line,mol,allownofile,verbose) 
    use struct_basic, only: crystal
    use struct_readers, only: struct_read_cif, struct_read_res, struct_read_cube,&
       struct_read_wien, struct_read_wien, struct_read_potcar, struct_read_vasp,&
       struct_read_abinit, struct_read_elk, struct_read_qeout, struct_read_qein,&
       struct_read_crystalout,&
       struct_read_library, struct_read_mol, struct_read_siesta, struct_read_dftbp,&
       struct_read_xsf, parse_crystal_env, parse_molecule_env, is_espresso
    use global, only: doguess, iunit_isdef, iunit, iunit_ang, iunit_bohr,&
       iunitname0, iunitname, dunit0, dunit, rborder_def, eval_next
    use tools_io, only: getword, equal, ferror, faterr, zatguess, lgetword,&
       string, uin, isinteger
    use param, only: dirsep, maxzat0

    character*(*), intent(in) :: line
    logical, intent(in) :: mol
    type(crystal), intent(inout) :: c
    logical, intent(in) :: allownofile
    logical, intent(in) :: verbose

    integer :: lp, lp2, istruct
    character(len=:), allocatable :: word, word2, wext1, wext2, subline, aux
    integer :: ntyp, nn
    character*5 :: ztyp(maxzat0)
    real*8 :: rborder, raux
    logical :: docube, ok

    ! Initialize the structure
    call c%end()
    call c%init()

    ! save whether this is a crystal or a molecule
    c%ismolecule = mol
    c%molx0 = 0d0

    ! enforce MOLECULE keyword particular settings
    if (c%ismolecule) then
       ! deactivate symmetry
       doguess = 0
       ! default unit is ang
       if (iunit_isdef) iunit = iunit_ang
    else
       if (iunit_isdef) iunit = iunit_bohr
    end if
    iunitname = trim(iunitname0(iunit))
    dunit = dunit0(iunit)

    ! read and parse
    lp=1
    word = getword(line,lp)
    aux = word(index(word,dirsep,.true.)+1:)
    wext1 = aux(index(aux,'.',.true.)+1:)
    aux = word(index(word,dirsep,.true.)+1:)
    wext2 = aux(index(aux,'_',.true.)+1:)
    c%file = ""

    if (equal(wext1,'cif')) then
       aux = getword(line,lp)
       if (equal(aux,"")) then
          aux = " "
       end if
       call struct_read_cif(c,word,aux,.false.,mol)
       call c%set_cryscar()
       c%file = word

    else if (equal(wext1,'res')) then
       call struct_read_res(c,word,.false.,mol)
       call c%set_cryscar()
       c%file = word

    else if (equal(wext1,'cube')) then
       call struct_read_cube(c,word,.false.,mol)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'struct')) then
       call struct_read_wien(c,word,.false.,mol)
       call c%set_cryscar()
       if (c%neqv == 0) then ! some structs may come without symmetry
          call struct_read_wien(c,word,.true.,mol)
          call c%set_cryscar()
       endif
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'POSCAR').or.equal(wext1,'CONTCAR').or.equal(wext1,'CHGCAR').or.&
             equal(wext1,'CHG').or.equal(wext1,'ELFCAR').or.equal(wext1,'AECCAR0').or.&
             equal(wext1,'AECCAR2')) then
       wext1 = getword(line,lp)
       wext2 = wext1(index(wext1,dirsep,.true.)+1:)
       wext2 = wext2(index(wext2,'.',.true.)+1:)
       if (equal(wext2,'POTCAR')) then
          call struct_read_potcar(wext1,ntyp,ztyp)
          aux = getword(line,lp)
          if (len_trim(aux) > 0) then
             call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
             return
          end if
       else
          ntyp = 0
          ztyp = ""
          do while(.true.)
             nn = zatguess(wext1)
             if (nn >= 0) then
                ntyp = ntyp + 1
                ztyp(ntyp) = string(wext1)
             else
                if (len_trim(wext1) > 0) then
                   call ferror('struct_crystal_input','Unknown atom type in CRYSTAL',faterr,line,syntax=.true.)
                   return
                end if
                exit
             end if
             wext1 = getword(line,lp)
          end do
       end if
       call struct_read_vasp(c,word,ntyp,ztyp,mol)
       c%file = word

    else if (equal(wext1,'DEN').or.equal(wext2,'DEN').or.equal(wext1,'ELF').or.equal(wext2,'ELF').or.&
       equal(wext1,'POT').or.equal(wext2,'POT').or.equal(wext1,'VHA').or.equal(wext2,'VHA').or.&
       equal(wext1,'VHXC').or.equal(wext2,'VHXC').or.equal(wext1,'VXC').or.equal(wext2,'VXC').or.&
       equal(wext1,'GDEN1').or.equal(wext2,'GDEN1').or.equal(wext1,'GDEN2').or.equal(wext2,'GDEN2').or.&
       equal(wext1,'GDEN3').or.equal(wext2,'GDEN3').or.equal(wext1,'LDEN').or.equal(wext2,'LDEN').or.&
       equal(wext1,'KDEN').or.equal(wext2,'KDEN').or.equal(wext1,'PAWDEN').or.equal(wext2,'PAWDEN')) then
       call struct_read_abinit(c,word,mol)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'OUT')) then
       call struct_read_elk(c,word,mol)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'out')) then
       if (is_espresso(word)) then
          ok = isinteger(istruct,line,lp)
          if (.not.ok) istruct = 0
          call struct_read_qeout(c,word,mol,istruct)
       else
          call struct_read_crystalout(c,word,mol)
       end if
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word
    else if (equal(wext1,'in')) then
       call struct_read_qein(c,word,mol)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(word,'library')) then
       subline = line(lp:)
       call struct_read_library(c,subline,mol,ok)
       if (.not.ok) return
       if (mol) then
          c%file = "molecular library (" // trim(line(lp:)) // ")"
       else
          c%file = "crystal library (" // trim(line(lp:)) // ")"
       end if

    else if (equal(wext1,'xyz').or.equal(wext1,'wfn').or.equal(wext1,'wfx').or.&
       equal(wext1,'fchk').or.equal(wext1,'molden')) then
       docube = .false.
       rborder = rborder_def 
       do while(.true.)
          lp2 = 1
          word2 = lgetword(line,lp)
          if (equal(word2,'cubic').or.equal(word2,'cube')) then
             docube = .true.
          elseif (eval_next(raux,word2,lp2)) then
             rborder = raux / dunit
          elseif (len_trim(word2) > 1) then
             call ferror('struct_crystal_input','Unknown extra keyword',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do

       call struct_read_mol(c,word,wext1,rborder,docube)
       call c%set_cryscar()
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'STRUCT_OUT') .or. equal(wext1,'STRUCT_IN')) then
       call struct_read_siesta(c,word,mol)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'xsf')) then
       call struct_read_xsf(c,word,mol)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (equal(wext1,'gen')) then
       docube = .false.
       rborder = rborder_def
       do while(.true.)
          lp2 = 1
          word2 = lgetword(line,lp)
          if (equal(word2,'cubic').or.equal(word2,'cube')) then
             docube = .true.
          elseif (eval_next(raux,word2,lp2)) then
             rborder = raux / dunit
          elseif (len_trim(word2) > 1) then
             call ferror('struct_crystal_input','Unknown extra keyword',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do

       call struct_read_dftbp(c,word,mol,rborder,docube)
       aux = getword(line,lp)
       if (len_trim(aux) > 0) then
          call ferror('struct_crystal_input','Unknown extra keyword in CRYSTAL',faterr,line,syntax=.true.)
          return
       end if
       c%file = word

    else if (len_trim(word) < 1) then
       if (.not.mol) then
          if (.not.allownofile) then
             call ferror('struct_crystal_input','Attempted to parse CRYSTAL environment but allownofile=.true.',faterr,syntax=.true.)
             return
          end if
          call parse_crystal_env(c,uin,ok)
          if (.not.ok) return
       else
          if (.not.allownofile) then
             call ferror('struct_crystal_input','Attempted to parse MOLECULE environment but allownofile=.true.',faterr,syntax=.true.)
             return
          end if
          call parse_molecule_env(c,uin,ok)
          if (.not.ok) return
       endif
       c%file = "<input>"

    else
       call ferror('struct_crystal_input','unrecognized file format',faterr,line,syntax=.true.)
       return
    end if

    call c%struct_fill(.true.,.true.,doguess,-1,.false.,.true.,.false.)
    if (verbose) call c%struct_report()

  end subroutine struct_crystal_input

  ! use the P1 space group
  subroutine struct_clearsym()
    use struct_basic, only: cr
    use types, only: atom, realloc
    use tools_io, only: uout
    use param, only: eyet
    type(atom) :: aux(cr%nneq)
    integer :: i

    ! nullify the space group
    cr%havesym = 0
    cr%neqv = 1
    cr%rotm(:,:,1) = eyet
    cr%ncv = 1
    if (allocated(cr%cen)) deallocate(cr%cen)
    allocate(cr%cen(3,4))
    cr%cen = 0d0

    ! convert ncel to nneq
    aux = cr%at(1:cr%nneq)
    cr%nneq = cr%ncel
    call realloc(cr%at,cr%ncel)
    do i = 1, cr%ncel
       cr%at(i) = aux(cr%atcel(i)%idx)
       cr%at(i)%x = cr%atcel(i)%x
    end do

    write (uout,'("* CLEARSYM: clear all symmetry operations and rebuild the atom list.")')
    write (uout,'("            The new crystal description follows"/)')
    call cr%struct_fill(.true.,.true.,0,-1,.false.,.true.,.false.)
    call cr%struct_report()

  end subroutine struct_clearsym

  subroutine struct_charges(line,oksyn)
    use struct_basic, only: cr
    use global, only: eval_next
    use tools_io, only: lgetword, equal, ferror, faterr, getword, zatguess

    character*(*), intent(in) :: line
    logical, intent(out) :: oksyn

    character(len=:), allocatable :: word
    integer :: lp, lp2, nn, i
    logical :: ok, do1
    real*8 :: xx

    oksyn = .false.
    lp = 1
    word = lgetword(line,lp)
    if (equal(word,'q') .or. equal(word,'zpsp') .or. equal(word,'qat')) then
       do1 = equal(word,'zpsp')
       do while (.true.)
          lp2 = lp
          ok = eval_next(nn,line,lp)
          if (ok) then
             ok = eval_next(xx,line,lp)
             if (.not.ok) then
                call ferror('struct_charges','Incorrect Q/QAT/ZPSP syntax',faterr,line,syntax=.true.)
                return
             end if
             if (do1) then
                cr%at(nn)%zpsp = nint(xx)
             else
                cr%at(nn)%qat = xx
             end if
          else 
             lp = lp2
             word = getword(line,lp)
             if (len_trim(word) < 1) exit
             nn = zatguess(word)
             if (nn == -1) then
                call ferror('struct_charges','Unknown atomic symbol in Q/QAT/ZPSP',faterr,line,syntax=.true.)
                return
             end if
             ok = eval_next(xx,line,lp)
             if (.not.ok) then
                call ferror('struct_charges','Incorrect Q/QAT/ZPSP syntax',faterr,line,syntax=.true.)
                return
             end if
             do i = 1, cr%nneq
                if (cr%at(i)%z == nn) then
                   if (do1) then
                      cr%at(i)%zpsp = nint(xx)
                   else
                      cr%at(i)%qat = xx
                   end if
                endif
             end do
          endif
       end do
    elseif (equal(word,'nocore')) then
       cr%at(1:cr%nneq)%zpsp = -1
       word = getword(line,lp)
       if (len_trim(word) > 0) then
          call ferror('critic','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       end if
    endif
    oksyn = .true.

  end subroutine struct_charges

  ! Write the crystal structure to a file
  subroutine struct_write(c,line)
    use struct_writers, only: struct_write_mol, struct_write_3dmodel, struct_write_gaussian,&
       struct_write_espresso, struct_write_vasp, struct_write_abinit, struct_write_elk,&
       struct_write_tessel, struct_write_critic, struct_write_cif, struct_write_escher,&
       struct_write_gulp, struct_write_lammps, struct_write_siesta_fdf, struct_write_siesta_in,&
       struct_write_dftbp_hsd, struct_write_dftbp_gen, struct_write_d12
    use struct_basic, only: crystal
    use global, only: eval_next, dunit
    use tools_io, only: getword, equal, lower, lgetword, ferror, faterr, uout, &
       string
    type(crystal), intent(inout) :: c
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word, wext, file, wroot
    integer :: lp, ix(3), lp2, iaux, nmer
    logical :: doborder, molmotif, dodreiding, docell, domolcell, ok
    logical :: onemotif, environ, lnmer
    real*8 :: rsph, xsph(3), rcub, xcub(3), renv

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
             renv = renv / dunit
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
             rsph = rsph / dunit
             if (c%ismolecule) then
                xsph = xsph / dunit - c%molx0
                xsph = c%c2x(xsph)
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
             rcub = rcub / dunit
             if (c%ismolecule) then
                xcub = xcub / dunit - c%molx0
                xcub = c%c2x(xcub)
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
          call struct_write_mol(c,file,wext,ix,doborder,onemotif,molmotif,&
             environ,renv,lnmer,nmer,rsph,xsph,rcub,xcub)
       else
          call struct_write_3dmodel(c,file,wext,ix,doborder,onemotif,molmotif,&
             docell,domolcell,rsph,xsph,rcub,xcub)
       end if
    elseif (equal(wext,'gau')) then
       ! gaussian periodic boundary conditions
       write (uout,'("* WRITE Gaussian file: ",A)') string(file)
       call struct_write_gaussian(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'in')) then
       ! espresso
       write (uout,'("* WRITE espresso file: ",A)') string(file)
       call struct_write_espresso(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'poscar') .or. equal(wext,'contcar')) then
       ! vasp
       write (uout,'("* WRITE VASP file: ",A)') string(file)
       call struct_write_vasp(file,c,.true.)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'abin')) then
       ! abinit
       write (uout,'("* WRITE abinit file: ",A)') string(file)
       call struct_write_abinit(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'elk')) then
       ! elk
       write (uout,'("* WRITE elk file: ",A)') string(file)
       call struct_write_elk(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'tess')) then
       ! tessel
       write (uout,'("* WRITE tess file: ",A)') string(file)
       call struct_write_tessel(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'incritic').or.equal(wext,'cri')) then
       ! critic2
       write (uout,'("* WRITE critic2 input file: ",A)') string(file)
       call struct_write_critic(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'cif')) then
       ! cif
       write (uout,'("* WRITE cif file: ",A)') string(file)
       call struct_write_cif(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'d12')) then
       ! d12
       write (uout,'("* WRITE crystal file: ",A)') string(file)
       call struct_write_d12(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'m')) then
       ! escher
       write (uout,'("* WRITE escher file: ",A)') string(file)
       call struct_write_escher(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'gin')) then
       ! gulp
       dodreiding = .false.
       do while(.true.)
          word = lgetword(line,lp)
          if (equal(word,'dreiding')) then
             dodreiding = .true.
          elseif (len_trim(word) > 1) then
             call ferror('struct_write','Unknown extra keyword',faterr,line,syntax=.true.)
             return
          else
             exit
          end if
       end do
       write (uout,'("* WRITE gulp file: ",A)') string(file)
       call struct_write_gulp(file,c,dodreiding)
    elseif (equal(wext,'lammps')) then
       write (uout,'("* WRITE lammps file: ",A)') string(file)
       call struct_write_lammps(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'fdf')) then
       write (uout,'("* WRITE fdf file: ",A)') string(file)
       call struct_write_siesta_fdf(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'struct_in')) then
       write (uout,'("* WRITE STRUCT_IN file: ",A)') string(file)
       call struct_write_siesta_in(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'hsd')) then
       write (uout,'("* WRITE hsd file: ",A)') string(file)
       call struct_write_dftbp_hsd(file,c)
       ok = check_no_extra_word()
       if (.not.ok) return
    elseif (equal(wext,'gen')) then
       write (uout,'("* WRITE gen file: ",A)') string(file)
       call struct_write_dftbp_gen(file,c)
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

  !> Calculate the powder diffraction pattern for the current
  !structure.
  subroutine struct_powder(line,c)
    use struct_basic, only: crystal
    use global, only: fileroot, eval_next
    use tools_io, only: ferror, faterr, uout, lgetword, equal, getword, &
       fopen_write, string, ioj_center, string, fclose
    use param, only: pi
    character*(*), intent(in) :: line
    type(crystal), intent(in) :: c

    integer :: i, lp, lu, np, npts
    real*8 :: th2ini, th2end, lambda, fpol, sigma
    real*8, allocatable :: t(:), ih(:), th2p(:), ip(:)
    integer, allocatable :: hvecp(:,:)
    character(len=:), allocatable :: root, word
    logical :: ok

    if (c%ismolecule) then
       call ferror("struct_powder","POWDER can not be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! default values
    th2ini = 5d0 
    th2end = 90d0 
    lambda = 1.5406d0
    fpol = 0d0 ! polarization, fpol = 0.95 for synchrotron
    npts = 10001
    sigma = 0.05d0
    root = trim(fileroot) // "_xrd"

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

    call c%powder(th2ini,th2end,npts,lambda,fpol,sigma,t,ih,th2p,ip,hvecp)
    np = size(th2p)

    ! write the data file 
    lu = fopen_write(trim(root) // ".dat")
    write (lu,'("# ",A,A)') string("2*theta",13,ioj_center), string("Intensity",13,ioj_center)
    do i = 1, npts
       write (lu,'(A,X,A)') string(t(i),"f",15,7,ioj_center), &
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
       write (uout,'(99(A,X))') string(i,3), string(hvecp(1,i),2), string(hvecp(2,i),2),&
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
  subroutine struct_rdf(line,c)
    use struct_basic, only: crystal
    use global, only: fileroot, eval_next
    use tools_io, only: faterr, ferror, uout, lgetword, equal, fopen_write,&
       ioj_center, getword, string, fclose
    character*(*), intent(in) :: line
    type(crystal), intent(in) :: c

    real*8 :: rend
    character(len=:), allocatable :: root, word
    logical :: ok
    integer :: npts
    real*8, allocatable :: t(:), ih(:)

    integer :: lp, lu, i
    ! integer :: i, lp, lu, np
    ! real*8 :: th2ini, th2end, lambda, fpol, sigma
    ! integer, allocatable :: hvecp(:,:)

    if (c%ismolecule) then
       call ferror("struct_rdf","RDF can not be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! default values
    rend = 25d0
    root = trim(fileroot) // "_rdf"
    npts = 10001

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
             call ferror('struct_powder','Incorrect TH2END',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,"npts")) then
          ok = eval_next(npts,line,lp)
          if (.not.ok) then
             call ferror('struct_powder','Incorrect NPTS',faterr,line,syntax=.true.)
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

    call c%rdf(rend,npts,t,ih)

    ! write the data file 
    lu = fopen_write(trim(root) // ".dat")
    write (lu,'("# ",A,A)') string("r (bohr)",13,ioj_center), string("RDF(r)",13,ioj_center)
    do i = 1, npts
       write (lu,'(A,X,A)') string(t(i),"f",15,7,ioj_center), &
          string(ih(i),"f",15,7,ioj_center)
    end do
    call fclose(lu)
    deallocate(t,ih)

    ! write the gnuplot input file
    lu = fopen_write(trim(root) // ".gnu")
    write (lu,'("set terminal postscript eps color enhanced ""Helvetica"" 25")')
    write (lu,'("set output """,A,".eps""")') trim(root)
    write (lu,*)
    write (lu,'("set xlabel ""r (bohr)""")')
    write (lu,'("set ylabel ""RDF(r)""")')
    write (lu,'("set xrange [0:",A,"]")') string(rend,"f")
    write (lu,'("set style data lines")')
    write (lu,'("set grid")')
    write (lu,'("unset key")')
    write (lu,'("plot """,A,".dat"" w lines")') string(root)
    write (lu,*)
    call fclose(lu)

    if (allocated(t)) deallocate(t)
    if (allocated(ih)) deallocate(ih)
    
  end subroutine struct_rdf

  !> Compare two crystal structures using the powder diffraction
  !> patterns or the radial distribution functions. Uses the
  !> similarity based on cross-correlation functions proposed in
  !>   de Gelder et al., J. Comput. Chem., 22 (2001) 273.
  subroutine struct_compare(line)
    use struct_basic, only: crystal, isformat_unknown, cr
    use struct_readers, only: struct_detect_format
    use global, only: doguess, eval_next, dunit, iunit, iunitname0
    use tools_math, only: crosscorr_triangle, rmsd_walker
    use tools_io, only: getword, equal, faterr, ferror, uout, string, ioj_center,&
       ioj_left, string
    use types, only: realloc
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word, tname, difstr, difout, diftyp
    integer :: doguess0
    integer :: lp, i, j, k, n
    integer :: ns, imol, isformat
    type(crystal), allocatable :: c(:)
    real*8 :: tini, tend, nor, h, xend
    real*8, allocatable :: t(:), ih(:), th2p(:), ip(:), iha(:,:)
    integer, allocatable :: hvecp(:,:)
    real*8, allocatable :: diff(:,:), xnorm(:), x1(:,:), x2(:,:)
    logical :: ok, usedot
    logical :: dopowder, ismol, laux
    character*1024, allocatable :: fname(:)

    real*8, parameter :: sigma0 = 0.2d0
    real*8, parameter :: lambda0 = 1.5406d0
    real*8, parameter :: fpol0 = 0d0
    integer, parameter :: npts = 10001
    real*8, parameter :: th2ini = 5d0
    real*8, parameter :: th2end0 = 50d0
    real*8, parameter :: rend0 = 25d0

    ! initialized
    doguess0 = doguess
    lp = 1
    ns = 0
    dopowder = .true.
    xend = -1d0
    allocate(fname(1))

    ! read the input optons
    doguess = 0
    imol = -1
    do while(.true.)
       word = getword(line,lp)
       if (equal(word,'xend')) then
          ok = eval_next(xend,line,lp)
          if (.not.ok) then
             call ferror('struct_compare','incorrect TH2END',faterr,syntax=.true.)
             return
          end if
       elseif (equal(word,'powder')) then
          dopowder = .true.
       elseif (equal(word,'rdf')) then
          dopowder = .false.
       elseif (equal(word,'molecule')) then
          imol = 1
       elseif (equal(word,'crystal')) then
          imol = 0
       elseif (len_trim(word) > 0) then
          ns = ns + 1
          if (ns > size(fname)) &
             call realloc(fname,2*ns)
          fname(ns) = word
       else
          exit
       end if
    end do
    if (ns < 2) &
       call ferror('struct_compare','At least 2 structures are needed for the comparison',faterr)

    ! determine whether to use crystal or molecule comparison
    ismol = .true.
    usedot = .false.
    do i = 1, ns
       if (.not.equal(fname(i),".")) then
          call struct_detect_format(fname(i),isformat,laux)
          ismol = ismol .and. laux
          if (isformat == isformat_unknown) &
             call ferror("struct_compare","unknown file format: " // string(fname(i)),faterr)
          inquire(file=fname(i),exist=laux)
          if (.not.laux) &
             call ferror("struct_compare","file not found: " // string(fname(i)),faterr)
       else
          usedot = .true.
       end if
    end do
    if (imol == 0) then
       ismol = .false.
    elseif (imol == 1) then
       ismol = .true.
    end if
    if (usedot .and. (cr%ismolecule .neqv. ismol)) &
       call ferror("struct_compare","current structure (.) incompatible with molecule/crystal in compare",faterr)
    if (usedot.and..not.cr%isinit) &
       call ferror('struct_compare','Current structure is not initialized.',faterr)
       
    ! Read the structures and header
    write (uout,'("* COMPARE: compare structures")')
    if (ismol) then
       tname = "Molecule"
    else
       tname = "Crystal"
    end if
    allocate(c(ns))
    do i = 1, ns
       if (equal(fname(i),".")) then
          write (uout,'("  ",A," ",A,": <current>")') string(tname), string(i,2)
          c(i) = cr
       else
          write (uout,'("  ",A," ",A,": ",A)') string(tname), string(i,2), string(fname(i)) 
          call struct_crystal_input(c(i),fname(i),ismol,.false.,.false.)
          if (.not.c(i)%isinit) &
             call ferror("struct_compare","could not load crystal structure" // string(fname(i)),faterr)
       end if
    end do

    ! rest of the header and default variables
    if (ismol) then
       diftyp = "Molecule"
       difstr = "RMS"
       write (uout,'("# RMS of the atomic positions in ",A)') iunitname0(iunit)
    else
       diftyp = "Crystal"
       difstr = "DIFF"
       if (dopowder) then
          write (uout,'("# Using cross-correlated POWDER diffraction patterns.")')
          if (xend < 0d0) xend = th2end0
       else
          write (uout,'("# Using cross-correlated radial distribution functions (RDF).")')
          if (xend < 0d0) xend = rend0
       end if
       write (uout,'("# Two structures are exactly equal if DIFF = 0.")')
    end if

    ! allocate space for difference/rms values
    allocate(diff(ns,ns))
    diff = 1d0

    if (.not.ismol) then
       ! crystals
       allocate(iha(10001,ns))
       do i = 1, ns
          ! calculate the powder diffraction pattern
          if (dopowder) then
             call c(i)%powder(th2ini,xend,npts,lambda0,fpol0,sigma0,t,ih,th2p,ip,hvecp)

             ! normalize the integral of abs(ih)
             tini = ih(1)**2
             tend = ih(npts)**2
             nor = (2d0 * sum(ih(2:npts-1)**2) + tini + tend) * (xend - th2ini) / 2d0 / real(npts-1,8)
             iha(:,i) = ih / sqrt(nor)
          else
             call c(i)%rdf(xend,npts,t,ih)
             iha(:,i) = ih
          end if
       end do
       if (allocated(t)) deallocate(t)
       if (allocated(ih)) deallocate(ih)
       if (allocated(th2p)) deallocate(th2p)
       if (allocated(ip)) deallocate(ip)
       if (allocated(hvecp)) deallocate(hvecp)
       
       ! self-correlation
       allocate(xnorm(ns))
       h =  (xend-th2ini) / real(npts-1,8)
       do i = 1, ns
          xnorm(i) = crosscorr_triangle(h,iha(:,i),iha(:,i),1d0)
       end do
       xnorm = sqrt(abs(xnorm))

       ! calculate the overlap between diffraction patterns
       diff = 0d0
       do i = 1, ns
          do j = i+1, ns
             diff(i,j) = max(1d0 - crosscorr_triangle(h,iha(:,i),iha(:,j),1d0) / xnorm(i) / xnorm(j),0d0)
             diff(j,i) = diff(i,j)
          end do
       end do
       deallocate(xnorm)
    else
       ! molecules
       diff = 0d0
       do i = 1, ns
          do j = i+1, ns
             if (c(i)%ncel == c(j)%ncel) then
                n = c(i)%ncel
                allocate(x1(3,n),x2(3,n))
                do k = 1, n
                   x1(:,k) = c(i)%atcel(k)%r + c(i)%molx0
                   x2(:,k) = c(j)%atcel(k)%r + c(j)%molx0
                end do
                diff(i,j) = rmsd_walker(x1,x2)
                deallocate(x1,x2)
             else
                diff(i,j) = -1d0
             end if
             diff(j,i) = diff(i,j)
          end do
       end do
       diff = diff * dunit
    endif
   
    ! write output 
    if (ns == 2) then
       write (uout,'("+ ",A," = ",A)') string(difstr), string(diff(1,2),'e',12,6)
    else
       ! do i = 1, ns
       !    do j = i+1, ns
       !       if (diff(i,j) < 0d0) then
       !          difout = "<not comparable>"
       !       else
       !          difout = string(diff(i,j),'f',10,6)
       !       end if
       !       write (uout,'("+ ",A,"(",A,": ",A," | ",A,": ",A,") = ",A)') &
       !          string(difstr), string(i), string(c(i)%file), string(j), &
       !          string(c(j)%file), difout
       !    end do
       ! end do
       do i = 0, (ns-1)/5
          write (uout,'(99(A,X))') string(diftyp,15,ioj_center), &
             (string(c(5*i+j)%file,15,ioj_center),j=1,min(5,ns-i*5))
          write (uout,'(99(A,X))') string(difstr,15,ioj_center), &
             (string(5*i+j,15,ioj_center),j=1,min(5,ns-i*5))
          do j = 1, ns
             write (uout,'(2X,99(A,X))') string(c(j)%file,15,ioj_left), &
                (string(diff(j,5*i+k),'f',15,7,3),k=1,min(5,ns-i*5))
          end do
          write (uout,*)
       end do
    endif
    write (uout,*)

    ! clean up
    deallocate(diff)
    doguess = doguess0

  end subroutine struct_compare

  !> Calculate the atomic environment of a point or all the 
  !> non-equivalent atoms in the unit cell.
  subroutine struct_environ(line)
    use struct_basic, only: cr
    use global, only: eval_next, dunit, iunit, iunitname0
    use tools_io, only: string, lgetword, equal, ferror, faterr, string, uout,&
       ioj_right, ioj_center, zatguess, isinteger
    character*(*), intent(in) :: line

    integer :: lp, lp2
    integer :: nn, i, j, k, l
    real*8 :: x0(3), xout(3), x0in(3)
    logical :: doatoms, ok
    character(len=:), allocatable :: word
    integer, allocatable :: nneig(:), wat(:)
    real*8, allocatable :: dist(:), xenv(:,:,:)

    integer :: iat, iby, iat_mode, iby_mode
    integer, parameter :: inone = 0
    integer, parameter :: iznuc = 1
    integer, parameter :: iid = 2

    nn = 10
    x0 = 0d0
    doatoms = .true.
    iat = 0
    iat_mode = inone
    iby = 0
    iby_mode = inone

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"shells")) then
          ok = eval_next(nn,line,lp)
          if (.not.ok) then
             call ferror('struct_environ','Wrong SHELLS syntax',faterr,line,syntax=.true.)
             return
          end if
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
          if (cr%ismolecule) &
             x0 = cr%c2x(x0 / dunit - cr%molx0)
       elseif (equal(word,"at")) then
          lp2 = lp
          word = lgetword(line,lp)
          iat = zatguess(word)
          if (iat < 0) then
             lp = lp2
             ok = isinteger(iat,line,lp)
             if (.not.ok) &
                call ferror('struct_environ','Syntax error in ENVIRON/AT',faterr,line,syntax=.true.)
             iat_mode = iid
          else
             iat_mode = iznuc
          end if
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
       elseif (len_trim(word) > 0) then
          call ferror('struct_environ','Unknown extra keyword',faterr,line,syntax=.true.)
          return
       else
          exit
       end if
    end do

    write (uout,'("* ENVIRON")')
    allocate(nneig(nn),wat(nn),dist(nn))
    if (doatoms) then
       write (uout,'("+ Atomic environments")')
       if (.not.cr%ismolecule) then
          write (uout,'("     Atom     neig       d(",A,")     nneq  type           position (cryst)")') &
             iunitname0(iunit)
       else
          write (uout,'("     Atom     neig       d(",A,")     nneq  type           position (",A,")")') &
             iunitname0(iunit), iunitname0(iunit)
       end if
       do i = 1, cr%nneq
          if (iat_mode == iid) then
             if (iat /= i) cycle
          elseif (iat_mode == iznuc) then
             if (iat /= cr%at(i)%z) cycle
          end if
          call cr%pointshell(cr%at(i)%x,nn,nneig,wat,dist,xenv)
          do j = 1, nn
             if (iby_mode == iid) then
                if (iby /= wat(j)) cycle
             elseif (iby_mode == iznuc) then
                if (iby /= cr%at(wat(j))%z) cycle
             end if
             xout = xenv(:,1,j)
             if (cr%ismolecule) xout = (cr%x2c(xout)+cr%molx0) * dunit
             if (j == 1) then
                write (uout,'(I3,1X,"(",A4,")",4X,I4,3X,F12.7,3X,I4,3X,A4,3A)') &
                   i, cr%at(i)%name, nneig(j), dist(j)*dunit, wat(j), &
                   cr%at(wat(j))%name, (string(xout(k),'f',12,7,ioj_right),k=1,3)
             else
                if (wat(j) /= 0) then
                   write (uout,'(5X,"...",6X,I4,3X,F12.7,3X,I4,3X,A4,3A)') &
                      nneig(j), dist(j)*dunit, wat(j), cr%at(wat(j))%name,&
                      (string(xout(k),'f',12,7,ioj_right),k=1,3)
                end if
             end if
          end do
       end do
       write (uout,*)
    else
       call cr%pointshell(x0,nn,nneig,wat,dist,xenv)
       ! List of atomic environments
       write (uout,'("+ Atomic environments of (",A,",",A,",",A,")")') &
          string(x0in(1),'f'), string(x0in(2),'f'), string(x0in(3),'f')
       if (.not.cr%ismolecule) then
          write (uout,'("     Atom     neig       d(",A,")     nneq  type           position (cryst)")') &
             iunitname0(iunit)
       else
          write (uout,'("     Atom     neig       d(",A,")     nneq  type           position (",A,")")') &
             iunitname0(iunit), iunitname0(iunit)
       end if
       do j = 1, nn
          if (iby_mode == iid) then
             if (iby /= wat(j)) cycle
          elseif (iby_mode == iznuc) then
             if (iby /= cr%at(wat(j))%z) cycle
          end if
          xout = xenv(:,1,j)
          if (cr%ismolecule) xout = (cr%x2c(xout)+cr%molx0) * dunit
          if (wat(j) /= 0) then
             write (uout,'(5X,"...",6X,I4,3X,F12.7,3X,I4,3X,A4,3A)') &
                nneig(j), dist(j)*dunit, wat(j), cr%at(wat(j))%name, &
                (string(xout(k),'f',12,7,ioj_right),k=1,3)
          end if
       end do
       write (uout,*)

       ! Detailed list of neighbors
       write (uout,'("+ Neighbors of (",A,",",A,",",A,")")') &
          string(x0in(1),'f'), string(x0in(2),'f'), string(x0in(3),'f')
       write (uout,'(" Atom Id           position (cryst. coords)      Distance (",A,")")') &
          iunitname0(iunit)
       do j = 1, nn
          if (wat(j) == 0) cycle
          if (iby_mode == iid) then
             if (iby /= wat(j)) cycle
          elseif (iby_mode == iznuc) then
             if (iby /= cr%at(wat(j))%z) cycle
          end if
          do k = 1, nneig(j)
             xout = xenv(:,k,j)
             if (cr%ismolecule) xout = (cr%x2c(xout)+cr%molx0) * dunit
             write (uout,'(6(A,X))') &
                string(cr%at(wat(j))%name,5,ioj_center), string(wat(j),2),&
                (string(xout(l),'f',12,7,ioj_right),l=1,3),&
                string(dist(j)*dunit,'f',12,5,ioj_right)
          end do
       end do
       write (uout,*)
    end if
    deallocate(nneig,wat,dist)
    if (allocated(xenv)) deallocate(xenv)
    
  end subroutine struct_environ

  !> Calculate the packing ratio of the crystal.
  subroutine struct_packing(line)
    use struct_basic, only: cr
    use global, only: eval_next
    use tools_io, only: ferror, faterr, uout, lgetword, equal, string
    use param, only: atmvdw

    character*(*), intent(in) :: line

    integer :: lp
    logical :: dovdw, found, ok
    character(len=:), allocatable :: word
    integer :: i, j, n(3), ii(3), ntot, iaux, idx
    integer :: lvec(3)
    real*8 :: prec, alpha, x(3), dist
    real*8 :: vout, dv

    if (cr%ismolecule) then
       call ferror("critic2","PACKING can not be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! default values    
    dovdw = .false.
    prec = 1d-1

    ! header
    write (uout,'("* PACKING")')

    ! parse input
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"vdw")) then
          dovdw = .true.
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
    
    if (.not.dovdw) then
       write (uout,'("+ Packing ratio (%): ",A)') string(cr%get_pack_ratio(),'f',10,4)
    else
       ! prepare the grid
       write (uout,'("+ Est. precision in the % packing ratio: ",A)') string(prec,'e',10,3)
       prec = prec * cr%omega / 100d0
       write (uout,'("+ Est. precision in the interstitial volume: ",A)') trim(string(prec,'f',100,4))
       alpha = (prec / cr%omega)**(1d0/3d0)
       n = ceiling(cr%aa / alpha)
       write (uout,'("+ Using a volume grid with # of nodes: ",3(A," "))') &
          (string(n(j)),j=1,3)
       
       ! run over all the points in the grid
       ntot = n(1)*n(2)*n(3)
       vout = 0d0
       dv = cr%omega / ntot

       !$omp parallel do reduction(+:vout) private(ii,iaux,x,found,idx,dist,lvec) schedule(dynamic)
       do i = 0, ntot-1
          ! unpack the index
          ii(1) = modulo(i,n(1))
          iaux = (i - ii(1)) / (n(1))
          ii(2) = modulo(iaux,n(2))
          iaux = (iaux - ii(2)) / (n(2))
          ii(3) = modulo(iaux,n(3))
          ! calculate the point in the cell
          x = real(ii,8) / real(n)

          found = .false.
          do j = 1, cr%nneq
             idx = j
             call cr%nearest_atom(x,idx,dist,lvec)
             found = (dist < atmvdw(cr%at(j)%z))
             if (found) exit
          end do
          if (.not.found) then
             vout = vout + dv
          end if
       end do
       !$omp end parallel do
       write (uout,'("+ Interstitial volume (outside vdw spheres): ",A)') &
          trim(string(vout,'f',100,4))
       write (uout,'("+ Cell volume: ",A)') trim(string(cr%omega,'f',100,4))
       write (uout,'("+ Packing ratio (%): ",A)') string((cr%omega-vout)/cr%omega*100,'f',10,4)
    end if
    write (uout,*)

  end subroutine struct_packing

  !> Build a new crystal from the current crystal by cell transformation
  subroutine struct_newcell(line)
    use struct_basic, only: cr
    use global, only: eval_next
    use tools_math, only: matinv
    use tools_io, only: ferror, faterr, lgetword, equal
    character*(*), intent(in) :: line

    character(len=:), allocatable :: word
    logical :: ok, doprim
    integer :: lp, lp2, dotyp, i
    real*8 :: x0(3,3), t0(3), rdum(4)
    logical :: doinv

    if (cr%ismolecule) then
       call ferror("struct_newcell","NEWCELL can not be used with molecules",faterr,syntax=.true.)
       return
    end if

    ! transform to the primitive?
    lp = 1
    doprim = .false.
    dotyp = 0
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,"standard")) then
          dotyp = 1
          doprim = .false.
       elseif (equal(word,"primitive")) then
          dotyp = 1
          doprim = .true.
       elseif (equal(word,"niggli")) then
          dotyp = 2
       elseif (equal(word,"delaunay")) then
          dotyp = 3
       else
          lp = 1
          exit
       end if
    end do

    if (dotyp == 1) then
       call cr%cell_standard(doprim,.true.)
    elseif (dotyp == 2) then
       call cr%cell_niggli(.true.)
    elseif (dotyp == 3) then
       call cr%cell_delaunay(.true.)
    end if
    if (dotyp > 0) return

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
       x0 = 0d0
       do i = 1, 3
          x0(i,i) = rdum(i)
       end do
    end if

    t0 = 0d0
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
    if (doinv) x0 = matinv(x0)

    call cr%newcell(x0,t0,.true.)

  end subroutine struct_newcell

  !> Try to determine the molecular cell from the crystal geometry
  subroutine struct_molcell(line)
    use struct_basic, only: cr
    use global, only: rborder_def, eval_next, dunit
    use tools_io, only: ferror, faterr, uout, string

    character*(*), intent(in) :: line
    integer :: i, j, lp
    real*8 :: xmin(3), xmax(3), rborder, raux
    logical :: ok

    if (.not.cr%ismolecule) then
       call ferror('struct_molcell','Use MOLCELL with MOLECULE, not CRYSTAL.',faterr,syntax=.true.)
       return
    end if
    if (any(abs(cr%bb-90d0) > 1d-5)) then
       call ferror('struct_molcell','MOLCELL only allowed for orthogonal cells.',faterr,syntax=.true.)
       return
    end if

    ! defaults
    rborder = rborder_def 

    ! read optional border input
    lp = 1
    ok = eval_next(raux,line,lp)
    if (ok) rborder = raux / dunit

    ! find the encompassing cube
    xmin = 1d40
    xmax = -1d30
    do i = 1, cr%ncel
       do j = 1, 3
          xmin(j) = min(xmin(j),cr%at(i)%x(j))
          xmax(j) = max(xmax(j),cr%at(i)%x(j))
       end do
    end do

    ! apply the border
    do j = 1, 3
       xmin(j) = max(xmin(j) - rborder / cr%aa(j),0d0)
       xmax(j) = min(xmax(j) + rborder / cr%aa(j),1d0)
       cr%molborder(j) = min(xmin(j),1d0-xmax(j))
    end do

    ! some output
    write (uout,'("* MOLCELL: set up a molecular cell")')
    write (uout,'("+ Limit of the molecule within the cell (cryst coords):")')
    write (uout,'("  a axis: ",A," -> ",A)') trim(string(cr%molborder(1),'f',10,4)), trim(string(1d0-cr%molborder(1),'f',10,4))
    write (uout,'("  b axis: ",A," -> ",A)') trim(string(cr%molborder(2),'f',10,4)), trim(string(1d0-cr%molborder(2),'f',10,4))
    write (uout,'("  c axis: ",A," -> ",A)') trim(string(cr%molborder(3),'f',10,4)), trim(string(1d0-cr%molborder(3),'f',10,4))
    
  end subroutine struct_molcell

end module struct
