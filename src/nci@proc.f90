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

! non-covalent interaction index plots
submodule (nci) proc
  implicit none

  !xx! private procedures
  ! function field_has_core(s,ifield)
  ! subroutine grid_frame(xmat,isortho,rchol,delta,ok)
  ! subroutine header_atoms(opts,s,nat,iz,xat)
  ! subroutine write_cube_header(lu,l1,l2,res,isortho,nat,iz,xat)
  ! subroutine write_cube_body(lu,n,c,scale)
  ! function check_no_extra_word(line,lp,routine)
  ! function read_fragment(c,lu) result(fr)
  ! subroutine promol_checks(opts,s,errmsg)
  ! subroutine promol_point(opts,s,x,dopromol,nfrag,res,i,j,k)

  ! NCI checkpoint file format identifier and version
  character(len=8), parameter :: chk_magic = "CRI2NCCK"
  integer, parameter :: chk_version = 1

  ! tolerance for the real data in the checkpoint fingerprint
  real*8, parameter :: chk_eps = 1d-10

  ! the cube bodies and the dat file are written in (i,j,k) order with
  ! k fastest, the transpose of the grid order; they are transposed
  ! through a buffer this many i-planes at a time
  integer, parameter :: nblk = 16

contains

  !> The NCIPLOT keyword: parse the options, calculate the NCI grids
  !> (or read them from the checkpoint file), and write the output
  !> files.
  module subroutine nciplot()
    use systemmod, only: sy
    use fieldmod, only: type_grid
    use global, only: fileroot, eval_next, dunit0, quiet, iunit, iunitname0, abort_requested
    use fragmentmod, only: realloc_fragment
    use tools_io, only: getline, lgetword, equal, uin, faterr, ferror, ucopy, &
       string, getword, uout, tictac
    type(nci_options) :: opts
    type(nci_result) :: res
    character(len=:), allocatable :: line, word, oname, errmsg, chkfile
    logical :: ok, ok2, nstep_from_grid
    integer :: lp
    integer :: i, j, k
    real*8, allocatable :: splot(:,:,:)
    logical(1), allocatable :: inter(:,:,:)
    real*8 :: sumq, sumqp, sumv, rdum1, rdum2
    logical :: lchk, usecore, usechk, domolmotif, periodic
    integer :: findlimits, ithres, istep
    type(fragment) :: fr0
    ! cube limits
    real*8 :: x0(3), x1(3), x0x(3), x1x(3), xinc(3), xmat(3,3), xv(3), xaux(3)
    integer :: nstep(3)
    ! rthres
    real*8 :: rthres

    ! default values
    call nci_options_init(opts,sy,sy%iref)
    if (sy%c%ismolecule) then
       x0 = sy%c%molborder
       x1 = 1d0 - sy%c%molborder
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
       periodic = .false.
    else
       x0 = 0d0
       x1 = 1d0
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
       periodic = .true.
    end if
    xinc = 0.1d0
    oname = fileroot
    nstep = -1
    findlimits = -1
    ithres = -1
    istep = -1
    usechk = .true.
    rthres = 2.d0
    domolmotif = .false.
    allocate(opts%fr(1))

    ! the interpolation keyword cannot use the Fourier transform if the
    ! field has a core contribution
    usecore = field_has_core(sy,sy%iref)

    do while(.true.)
       ok = getline(uin,line,.true.,ucopy)
       lp = 1
       word = lgetword(line,lp)

       if (equal(word,'oname')) then
          oname = line(lp:)
       elseif (equal(word,'molmotif')) then
          domolmotif = .true.
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'cutoffs')) then
          ok = eval_next(opts%rhocut,line,lp)
          ok = ok .and. eval_next(opts%dimcut,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong cutoffs keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'rhoparam')) then
          ok = eval_next(opts%rhoparam,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong rhoparam keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'rhoparam2')) then
          ok = eval_next(opts%rhoparam2,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong rhoparam2 keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'cutplot')) then
          ok = eval_next(opts%rhoplot,line,lp)
          ok = ok .and. eval_next(opts%dimplot,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong cutplot keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'srhorange')) then
          ok = eval_next(rdum1,line,lp)
          ok = ok .and. eval_next(rdum2,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong srhorange keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return

          opts%srhorange(1) = min(rdum1,rdum2)
          opts%srhorange(2) = max(rdum1,rdum2)
       elseif (equal(word,'void')) then
          ok = eval_next(opts%rho_void,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong void keyword',faterr,line,syntax=.true.)
             return
          end if
          opts%dovoid = (opts%rho_void > 0d0)
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'rthres')) then
          ok = eval_next(rthres,line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong rthres keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          rthres = rthres / dunit0(iunit)
       elseif (equal(word,'increments')) then
          istep = 0
          ok = eval_next(xinc(1),line,lp)
          ok = ok .and. eval_next(xinc(2),line,lp)
          ok = ok .and. eval_next(xinc(3),line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong increments keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          xinc = xinc / dunit0(iunit)
       elseif (equal(word,'nstep')) then
          istep = 1
          ok = eval_next(nstep(1),line,lp)
          ok = ok .and. eval_next(nstep(2),line,lp)
          ok = ok .and. eval_next(nstep(3),line,lp)
          if (.not.ok) then
             call ferror('nciplot','wrong nstep keyword',faterr,line,syntax=.true.)
             return
          end if
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          do i = 1, 3
             nstep(i) = max(nstep(i),2)
          end do
       elseif (equal(word,'onlyneg')) then
          opts%onlyneg = .true.
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'nochk')) then
          usechk = .false.
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
       elseif (equal(word,'cube')) then
          findlimits = 0
          periodic = .false.
          ok = eval_next(x0(1),line,lp)
          if (ok) then
             ok = ok .and. eval_next(x0(2),line,lp)
             ok = ok .and. eval_next(x0(3),line,lp)
             ok = ok .and. eval_next(x1(1),line,lp)
             ok = ok .and. eval_next(x1(2),line,lp)
             ok = ok .and. eval_next(x1(3),line,lp)
             if (.not.ok) then
                call ferror('nciplot','wrong cube syntax',faterr,line,syntax=.true.)
                return
             end if
             ok = check_no_extra_word(line,lp,'nciplot')
             if (.not.ok) return
             if (.not.sy%c%ismolecule) then
                x0 = sy%c%x2c(x0)
                x1 = sy%c%x2c(x1)
             else
                x0 = x0 / dunit0(iunit) - sy%c%molx0
                x1 = x1 / dunit0(iunit) - sy%c%molx0
             endif
             ithres = 0
          else
             x0 = 1d30
             x1 = -1d30
             ok = .false.
             do while (.true.)
                word = getword(line,lp)
                inquire(file=word,exist=ok2)
                if (.not.ok2) exit
                ok = .true.
                fr0 = sy%c%identify_fragment_from_xyz(word,errmsg)
                if (len_trim(errmsg) > 0) call ferror('nciplot',errmsg,faterr)
                do j = 1, fr0%nat
                   do k = 1, 3
                      x0(k) = min(x0(k),fr0%at(j)%r(k))
                      x1(k) = max(x1(k),fr0%at(j)%r(k))
                   end do
                end do
             end do
             ithres = 1
          endif
          if (.not.ok) then
             call ferror('nciplot','wrong cube keyword',faterr,line,syntax=.true.)
             return
          end if
       elseif (equal(word,'fragment')) then
          if (findlimits == -1) findlimits = 1
          if (ithres == -1) ithres = 1
          opts%nfrag = opts%nfrag + 1
          if (opts%nfrag > size(opts%fr)) call realloc_fragment(opts%fr,2*opts%nfrag)
          word = getword(line,lp)
          if (equal(word,'')) then
             opts%fr(opts%nfrag) = read_fragment(sy%c,uin)
          else
             opts%fr(opts%nfrag) = sy%c%identify_fragment_from_xyz(word,errmsg)
             if (len_trim(errmsg) > 0) call ferror('nciplot',errmsg,faterr)
          endif

       elseif (equal(word,'interpolation')) then
          if (sy%f(sy%iref)%type /= type_grid) then
             call ferror('nciplot','interpolation keyword can only be used if the reference field is a grid',&
                faterr,line,syntax=.true.)
             return
          end if

          word = lgetword(line,lp)
          if (equal(word,'fourier')) then
             if (usecore) then
                call ferror('nciplot','interpolation keyword cannot be used with reference fields&
                   & that have a core contribution',faterr,line,syntax=.true.)
                return
             end if
             opts%fourierint = .true.
          elseif (equal(word,'grid')) then
             opts%fourierint = .false.
          else
             call ferror('nciplot','wrong interpolation keyword',faterr,line,syntax=.true.)
             return
          end if

       elseif (equal(word,'endnciplot').or.equal(word,'end')) then
          ok = check_no_extra_word(line,lp,'nciplot')
          if (.not.ok) return
          exit
       else
          call ferror ('nciplot','unknown keyword',faterr,line,syntax=.true.)
          return
       end if
    end do
    if (opts%nfrag > 0) then
       call realloc_fragment(opts%fr,opts%nfrag)
    else
       deallocate(opts%fr)
    end if

    ! header
    write (uout,'("* NCIPLOT: non-covalent interactions")')
    write (uout,'("  Please cite:")')
    write (uout,'("    E. R. Johnson et al., J. Am. Chem. Soc., 132 (2010) 6498. (10.1021/ja100936w)")')
    write (uout,'("    A. Otero-de-la-Roza et al., Phys. Chem. Chem. Phys., 14 (2012) 12165. (10.1039/C2CP41395G)")')
    if (.not.quiet) call tictac("start nciplot")

    ! redefine the cube using the fragments
    if (findlimits == 1) then
       if (opts%nfrag < 1) call ferror('nciplot','no fragments but findlimits=1',faterr)
       x0 = 1d30
       x1 = -1d30
       do i = 1, opts%nfrag
          do j = 1, opts%fr(i)%nat
             do k = 1, 3
                x0(k) = min(x0(k),opts%fr(i)%at(j)%r(k))
                x1(k) = max(x1(k),opts%fr(i)%at(j)%r(k))
             end do
          end do
       end do
       periodic = .false.
    endif

    ! apply the threshold
    if (ithres == 1) then
       x0 = x0 - rthres
       x1 = x1 + rthres
    end if

    ! convert to crystallographic coordinates
    x0 = sy%c%c2x(x0)
    x1 = sy%c%c2x(x1)
    x0x = x0
    x1x = x1

    ! determine the grid step
    nstep_from_grid = .false.
    if (periodic) then
       ! whole cells: the grid of the reference field, or xinc
       if (istep == -1 .and. sy%f(sy%iref)%type == type_grid) then
          nstep_from_grid = .true.
          nstep = sy%f(sy%iref)%grid%n
       end if
       do i = 1, 3
          xv = 0d0
          xv(i) = 1d0
          xv = sy%c%x2c(xv)
          if (.not.nstep_from_grid .and. istep /= 1) nstep(i) = ceiling(norm2(xv) / xinc(i))
          xmat(:,i) = xv / real(nstep(i),8)
          xinc(i) = norm2(xmat(:,i))
       end do
       x0 = sy%c%x2c(x0)
    else
       ! the box, with the corners sorted
       x0 = sy%c%x2c(x0)
       x1 = sy%c%x2c(x1)
       xaux = x0
       x0 = min(xaux,x1)
       x1 = max(xaux,x1)

       ! use the grid step
       if (istep == 0 .or. istep == -1) then
          nstep = ceiling(abs(x1 - x0) / xinc)
       else
          xinc = (x1-x0) / real(nstep,8)
       end if
       xmat = 0d0
       do i = 1, 3
          xmat(i,i) = xinc(i)
       end do
    end if
    opts%periodic = periodic
    opts%x0 = x0
    opts%xmat = xmat
    opts%nstep = nstep

    ! do stuff
    write(uout,'("+ File root: ",A)') trim(oname)
    write(uout,'("  Density cube file: ",A)') trim(oname) // "-dens.cube"
    write(uout,'("  RDG cube file: ",A)') trim(oname) // "-grad.cube"
    if (nci_dopromol(opts)) then
       write(uout,'("  Promolecular density cube file: ",A)') trim(oname) // "-pdens.cube"
    end if
    write(uout,'("  RDG vs. density file: ",A)') trim(oname) // ".dat"
    write(uout,'("  VMD visualization script: ",A)') trim(oname) // ".vmd"
    write(uout,'("+ Density cutoff (for dat file): ",A)') string(opts%rhocut,'f',decimal=6)
    write(uout,'("  RDG cutoff (for dat file): ",A)') string(opts%dimcut,'f',decimal=6)
    write(uout,'("  Density plot cutoff: ",A)') string(opts%rhoplot,'f',decimal=6)
    write(uout,'("  RDG plot cutoff: ",A)') string(opts%dimplot,'f',decimal=6)
    write(uout,'("+ Cube origin (crystal coord.): ",3(A," "))') (string(x0x(j),'f',decimal=6),j=1,3)
    write(uout,'("  Cube end (crystal coord.): ",3(A," "))') (string(x1x(j),'f',decimal=6),j=1,3)
    write(uout,'("  Cube nodes in each direction: ",3(A," "))') (string(nstep(j)),j=1,3)
    write(uout,'("  Cube steps in each direction: ",3(A," "))') (string(xinc(j),'f',decimal=6),j=1,3)
    write(uout,'("  Cube side lengths in each direction (",A,"): ",3(A," "))') &
       iunitname0(iunit), (string(xinc(j)*nstep(j)*dunit0(iunit),'f',decimal=4),j=1,3)
    if (nstep_from_grid) &
       write(uout,'("  Values of NSTEP derived from reference field grid dimensions.")')
    if (opts%fourierint) then
       write(uout,'("  Using fast Fourier transform for derivative calculation.")')
    else
       write(uout,'("  Using scalar field differentiation for derivative calculation.")')
    end if
    write(uout,*)

    ! calculate density, rdg,... read from chkpoint if available
    chkfile = trim(oname) // ".chk_nci"
    lchk = .false.
    if (usechk) then
       inquire(file=chkfile,exist=lchk)
       if (lchk) then
          write(uout,'("  Reading the checkpoint file: ",A/)') trim(chkfile)
          lchk = nci_chk_read(chkfile,opts,sy,res)
          if (.not.lchk) write(uout,'("  Checkpoint file is not compatible -- restarting."/)')
       end if
    end if
    if (.not.lchk) then
       call nci_compute(opts,sy,res,errmsg)
       if (len_trim(errmsg) > 0) call ferror('nciplot',errmsg,faterr)

       ! a cancelled calculation (the GUI) writes nothing
       if (abort_requested) then
          write (uout,'("+ The calculation was cancelled: no output files were written."/)')
          return
       end if

       ! save the ncichk file
       if (usechk) then
          write(uout,'("  Writing the checkpoint file: ",A/)') trim(chkfile)
          call nci_chk_write(chkfile,opts,sy,res,errmsg)
          if (len_trim(errmsg) > 0) call ferror('nciplot',errmsg,faterr)
       end if
    end if

    ! apply the cutoffs
    call nci_masks(opts,res,inter,splot,sumq,sumqp,sumv,errmsg)
    if (len_trim(errmsg) > 0) call ferror('nciplot',errmsg,faterr)
    if (opts%dovoid) then
       write (uout,'("* Void charge (a.u.): ",A)') string(sumq,'f',decimal=6)
       write (uout,'("* Void promolecular charge (a.u.): ",A)') string(sumqp,'f',decimal=6)
       write (uout,'("* Void volume (bohr^3): ",A)') string(sumv,'f',decimal=6)
    endif

    ! write the output files
    call nci_write_files(oname,opts,sy,res,inter,splot,domolmotif,errmsg)
    if (len_trim(errmsg) > 0) call ferror('nciplot',errmsg,faterr)

    if (.not.quiet) call tictac("End NCIPLOT")

  end subroutine nciplot

  !> Set the default options of an NCI calculation of the density of
  !> field ifield in system s: the cutoffs depend on whether the
  !> density is promolecular, and a grid field without core
  !> contribution has its derivatives calculated by Fourier transform.
  module subroutine nci_options_init(opts,s,ifield)
    use systemmod, only: system
    type(nci_options), intent(out) :: opts
    type(system), intent(in) :: s
    integer, intent(in) :: ifield

    call nci_options_set_field(opts,s,ifield)

  end subroutine nci_options_init

  !> Use field ifield of system s as the density of options opts: the
  !> field and the defaults that depend on it (the cutoffs of a
  !> promolecular density differ, and a grid without core densities
  !> is differentiated by Fourier transform). The other options are
  !> left as they are.
  module subroutine nci_options_set_field(opts,s,ifield)
    use systemmod, only: system
    use fieldmod, only: type_grid
    type(nci_options), intent(inout) :: opts
    type(system), intent(in) :: s
    integer, intent(in) :: ifield

    opts%ifield = ifield
    if (ifield /= 0) then
       opts%dimcut = 2.0d0
       opts%dimplot = 0.5d0
       opts%rhoplot = 0.1d0
    else
       opts%dimcut = 1.0d0
       opts%dimplot = 0.3d0
       opts%rhoplot = 0.12d0
    end if
    opts%fourierint = (s%f(ifield)%type == type_grid .and. .not.field_has_core(s,ifield))

  end subroutine nci_options_set_field

  !> Whether an NCI calculation with options opts calculates the
  !> promolecular density grid: only when the fragment or void
  !> criteria need it.
  pure module function nci_dopromol(opts) result(dopromol)
    type(nci_options), intent(in) :: opts
    logical :: dopromol

    dopromol = opts%nfrag > 0 .or. opts%dovoid

  end function nci_dopromol

  !> Calculate the NCI grids of system s with options opts: the
  !> density times the sign of the second Hessian eigenvalue, the
  !> reduced density gradient, and, if nci_dopromol, the promolecular
  !> density of the system and of each fragment. If pre is given and
  !> holds the promolecular grids for this grid and these fragments
  !> (nci_compute_promol; the caller checks the fragments), they are
  !> moved into res instead of being calculated again. Reports
  !> progress and can be cancelled (abort_poll); a cancelled
  !> calculation leaves the grids incomplete. errmsg is non-empty on
  !> error.
  module subroutine nci_compute(opts,s,res,errmsg,pre)
    use systemmod, only: system
    use fieldmod, only: field, type_grid
    use global, only: abort_poll, progress_start, progress_step
    use tools_math, only: eigsym
    use types, only: scalar_value, field_evaluation_avail
    use param, only: pi, vsmall, ifformat_as_ft_grad, ifformat_as_ft_xx,&
       ifformat_as_ft_yy, ifformat_as_ft_zz
    type(nci_options), intent(in) :: opts
    type(system), intent(inout), target :: s
    type(nci_result), intent(out) :: res
    character(len=:), allocatable, intent(out) :: errmsg
    type(nci_result), intent(inout), optional :: pre

    type(field) :: fgrho, fxx(3)
    type(scalar_value) :: sv, svg
    type(field_evaluation_avail) :: request0, request2
    logical :: dopromol
    integer :: n(3), i, j, k, l, istat, ier, nfrag, ierdiag
    real*8 :: x(3), ehess(3), dimgrad

    ! named constants
    real*8, parameter :: fthirds = 4d0/3d0
    real*8, parameter :: const = 2.d0*(3.d0*pi**2)**(1.d0/3.d0)

    ! checks
    n = opts%nstep
    nfrag = opts%nfrag
    dopromol = nci_dopromol(opts)
    call promol_checks(opts,s,errmsg)
    if (len_trim(errmsg) > 0) return
    if (.not.s%goodfield(opts%ifield)) then
       errmsg = "the density field is not available"
       return
    end if
    if (opts%fourierint .and. s%f(opts%ifield)%type /= type_grid) then
       errmsg = "Fourier transform derivatives need a grid field"
       return
    end if

    ! allocate memory for the grids; the promolecular grids of pre
    ! are taken if they are for this grid
    res%ifield = opts%ifield
    res%fourierint = opts%fourierint
    res%x0 = opts%x0
    res%xmat = opts%xmat
    res%n = n
    allocate(res%srho(n(1),n(2),n(3)),res%rdg(n(1),n(2),n(3)),stat=istat)
    if (istat == 0 .and. dopromol .and. present(pre)) then
       if (promol_matches(pre)) then
          call move_alloc(pre%rhoat,res%rhoat)
          if (nfrag > 0) call move_alloc(pre%rhofrag,res%rhofrag)
          dopromol = .false.
          nfrag = 0
       end if
    end if
    if (istat == 0 .and. dopromol) allocate(res%rhoat(n(1),n(2),n(3)),stat=istat)
    if (istat == 0 .and. nfrag > 0) allocate(res%rhofrag(n(1),n(2),n(3),nfrag),stat=istat)
    if (istat /= 0) then
       errmsg = "could not allocate memory for the NCI grids"
       return
    end if
    call request0%field_only()
    call request2%field_nder2()

    ! Initialize gradrho out of the omp loop
    if (opts%fourierint) then
       call fgrho%load_as_fftgrid(s%c,-1,"",s%f(opts%ifield)%grid,ifformat_as_ft_grad)
       call fxx(1)%load_as_fftgrid(s%c,-1,"",s%f(opts%ifield)%grid,ifformat_as_ft_xx)
       call fxx(2)%load_as_fftgrid(s%c,-1,"",s%f(opts%ifield)%grid,ifformat_as_ft_yy)
       call fxx(3)%load_as_fftgrid(s%c,-1,"",s%f(opts%ifield)%grid,ifformat_as_ft_zz)
       ! Set trilinear interpolation to prevent wonky negative values in the derivatives
       ! of the field in the low-density regions (which may have been noisy to begin with).
       call fgrho%grid%setmode("trilinear")
       do l = 1, 3
          call fxx(l)%grid%setmode("trilinear")
       end do
    endif

    ! calculate density, rdg, rhoat
    ierdiag = 0
    call progress_start(n(3)*n(2),"lines")
    !$omp parallel do private (x,sv,svg,ehess,dimgrad,ier) schedule(dynamic)
    do k = 1, n(3)
       do j = 1, n(2)
          if (abort_poll()) cycle
          do i = 1, n(1)
             x = opts%x0 + (i-1)*opts%xmat(:,1) + (j-1)*opts%xmat(:,2) + (k-1)*opts%xmat(:,3)

             ! calculate properties at x: rho and rdg
             if (opts%fourierint) then
                call s%f(opts%ifield)%grd(x,request0,sv)
                call fgrho%grd(x,request0,svg)
                dimgrad = svg%f / (const*max(sv%f,vsmall)**fthirds)
                do l = 1, 3
                   call fxx(l)%grd(x,request0,svg)
                   ehess(l) = svg%f
                end do
                if (count(ehess > 0d0) >= 2) then
                   ehess(2) = 1d0
                else
                   ehess(2) = -1d0
                endif
             else
                ! strangely enough, eig takes about the same time as counting the signs
                ! and using sylvester's law of inertia.
                call s%f(opts%ifield)%grd(x,request2,sv)
                call eigsym(sv%hf,3,ehess,ier)
                if (ier /= 0) then
                   !$omp atomic write
                   ierdiag = 1
                end if
                dimgrad = sv%gfmod / (const*max(sv%f,vsmall)**fthirds)
             end if
             res%rdg(i,j,k) = dimgrad
             res%srho(i,j,k) = sign(sv%f,ehess(2))

             ! promolecular density of the system and the fragments
             call promol_point(opts,s,x,dopromol,nfrag,res,i,j,k)
          end do
          call progress_step()
       end do
    end do
    !$omp end parallel do
    if (ierdiag /= 0) errmsg = "error in the diagonalization of the Hessian"

  contains

    !> Whether p holds the promolecular grids of this calculation.
    function promol_matches(p) result(ok)
      type(nci_result), intent(in) :: p
      logical :: ok

      ok = all(p%n == n) .and. all(p%x0 == opts%x0) .and. all(p%xmat == opts%xmat)
      if (ok) ok = allocated(p%rhoat)
      if (ok .and. nfrag > 0) ok = allocated(p%rhofrag)
      if (ok .and. nfrag > 0) ok = (size(p%rhofrag,4) == nfrag)

    end function promol_matches

  end subroutine nci_compute

  !> Calculate only the promolecular grids of an NCI calculation of
  !> system s with options opts on its grid: the promolecular density
  !> of the system and of each fragment (res%rhoat, res%rhofrag), for
  !> the fragment criterion, without the density field. A later
  !> nci_compute with the same grid and fragments can take them (pre).
  !> Reports progress and can be cancelled (abort_poll). errmsg is
  !> non-empty on error.
  module subroutine nci_compute_promol(opts,s,res,errmsg)
    use systemmod, only: system
    use global, only: abort_poll, progress_start, progress_step
    type(nci_options), intent(in) :: opts
    type(system), intent(inout), target :: s
    type(nci_result), intent(out) :: res
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: n(3), i, j, k, istat, nfrag
    real*8 :: x(3)

    n = opts%nstep
    nfrag = opts%nfrag
    call promol_checks(opts,s,errmsg)
    if (len_trim(errmsg) > 0) return

    res%x0 = opts%x0
    res%xmat = opts%xmat
    res%n = n
    allocate(res%rhoat(n(1),n(2),n(3)),stat=istat)
    if (istat == 0 .and. nfrag > 0) allocate(res%rhofrag(n(1),n(2),n(3),nfrag),stat=istat)
    if (istat /= 0) then
       errmsg = "could not allocate memory for the promolecular grids"
       return
    end if

    call progress_start(n(3)*n(2),"lines")
    !$omp parallel do private (x) schedule(dynamic)
    do k = 1, n(3)
       do j = 1, n(2)
          if (abort_poll()) cycle
          do i = 1, n(1)
             x = opts%x0 + (i-1)*opts%xmat(:,1) + (j-1)*opts%xmat(:,2) + (k-1)*opts%xmat(:,3)
             call promol_point(opts,s,x,.true.,nfrag,res,i,j,k)
          end do
          call progress_step()
       end do
    end do
    !$omp end parallel do

  end subroutine nci_compute_promol

  !> The fragment criterion of options opts as a smooth field on the
  !> promolecular grids of res: g = min(rhoparam - max_i rho_i/sum,
  !> sum/rhoat - rhoparam2), with sum the promolecular density of all
  !> the fragments, so the points that pass the criterion (nci_masks)
  !> are those with g >= 0. Points with no fragment density (beyond the
  !> cutoffs of all the atoms) get -1: the masks keep them if there is
  !> no promolecular density either, but there is no NCI plot in the
  !> vacuum, and the region would end in sheets at the cutoffs.
  module subroutine nci_frag_region(opts,res,g)
    use param, only: vsmall
    type(nci_options), intent(in) :: opts
    type(nci_result), intent(in) :: res
    real*8, allocatable, intent(out) :: g(:,:,:)

    integer :: n(3), i, j, k, nfrag
    real*8 :: sumf

    n = res%n
    allocate(g(n(1),n(2),n(3)))
    g = -1d0
    if (.not.allocated(res%rhofrag) .or. .not.allocated(res%rhoat)) return
    nfrag = size(res%rhofrag,4)
    if (nfrag < 1) return

    !$omp parallel do private (sumf) schedule(static)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             sumf = sum(res%rhofrag(i,j,k,1:nfrag))
             if (sumf <= 0d0) cycle
             g(i,j,k) = min(opts%rhoparam - maxval(res%rhofrag(i,j,k,1:nfrag)) / sumf,&
                sumf / max(res%rhoat(i,j,k),vsmall) - opts%rhoparam2)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine nci_frag_region

  !> Apply the fragment and void criteria and the plot cutoffs of
  !> options opts to the NCI grids res. On output, inter (if present)
  !> is true at the grid points that pass the fragment and void criteria, and
  !> splot is the reduced density gradient with the value
  !> nci_rdg_masked at the points removed from the isosurface plot (by
  !> those criteria, the density cutoff, the signed density range, or
  !> the only-negative option). With a void cutoff, voidq, voidqp, and
  !> voidv are the density and promolecular density integrated over
  !> the void and the void volume (a.u.); zero otherwise. errmsg is
  !> non-empty if the criteria need grids res does not have.
  module subroutine nci_masks(opts,res,inter,splot,voidq,voidqp,voidv,errmsg)
    use tools_math, only: det3
    type(nci_options), intent(in) :: opts
    type(nci_result), intent(in) :: res
    logical(1), allocatable, intent(out), optional :: inter(:,:,:)
    real*8, allocatable, intent(out) :: splot(:,:,:)
    real*8, intent(out) :: voidq, voidqp, voidv
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: n(3), i, j, k, nfrag
    real*8 :: sumf, dv
    logical :: lint, dointer

    ! checks
    errmsg = ""
    voidq = 0d0
    voidqp = 0d0
    voidv = 0d0
    n = res%n
    nfrag = opts%nfrag
    if (opts%dovoid .and. .not.allocated(res%rhoat)) then
       errmsg = "the void criterion needs the promolecular density grid: recalculate"
       return
    end if
    if (nfrag > 0) then
       if (.not.allocated(res%rhofrag) .or. .not.allocated(res%rhoat)) then
          errmsg = "the fragment criterion needs the fragment density grids: recalculate"
          return
       elseif (size(res%rhofrag,4) /= nfrag) then
          errmsg = "the fragments changed since the grids were calculated: recalculate"
          return
       end if
    end if
    dointer = present(inter)
    if (dointer) allocate(inter(n(1),n(2),n(3)))
    allocate(splot(n(1),n(2),n(3)))

    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             ! the fragment and void criteria
             lint = .true.
             if (nfrag > 0) then
                sumf = sum(res%rhofrag(i,j,k,1:nfrag))
                lint = (sumf >= opts%rhoparam2 * res%rhoat(i,j,k)) .and. &
                   all(res%rhofrag(i,j,k,1:nfrag) <= sumf*opts%rhoparam)
             end if
             if (opts%dovoid) then
                if (res%rhoat(i,j,k) < opts%rho_void) then
                   voidq = voidq + abs(res%srho(i,j,k))
                   voidqp = voidqp + res%rhoat(i,j,k)
                   voidv = voidv + 1
                else
                   lint = .false.
                end if
             end if
             if (dointer) inter(i,j,k) = lint

             ! the plot cutoffs
             if (.not.lint .or. abs(res%srho(i,j,k)) > opts%rhoplot .or.&
                res%srho(i,j,k) < opts%srhorange(1) .or. res%srho(i,j,k) > opts%srhorange(2) .or.&
                (opts%onlyneg .and. res%srho(i,j,k) > 0)) then
                splot(i,j,k) = nci_rdg_masked
             else
                splot(i,j,k) = res%rdg(i,j,k)
             end if
          end do
       end do
    end do

    ! the volume element of the grid
    dv = abs(det3(res%xmat))
    voidq = voidq * dv
    voidqp = voidqp * dv
    voidv = voidv * dv

  end subroutine nci_masks

  !> Whether grid point (i,j,k) of the NCI grids res is in the s vs.
  !> density plot (the dat file): it passes the fragment and void
  !> criteria (inter, from nci_masks) and the density and RDG cutoffs
  !> of options opts.
  pure module function nci_in_dat(opts,res,inter,i,j,k) result(ok)
    type(nci_options), intent(in) :: opts
    type(nci_result), intent(in) :: res
    logical(1), intent(in) :: inter(:,:,:)
    integer, intent(in) :: i, j, k
    logical :: ok

    ok = inter(i,j,k) .and. abs(res%srho(i,j,k)) < opts%rhocut .and. res%rdg(i,j,k) < opts%dimcut

  end function nci_in_dat

  !> The domains of the isosurface plot of the NCI grids res with
  !> options opts: the connected regions of the grid where the masked
  !> reduced density gradient splot (nci_masks) is below the isovalue
  !> (grid_domains; wrapping around on a periodic grid). On output,
  !> ndom domains numbered by decreasing size, lbl the domain of each
  !> grid point (0 = none), and for each domain its number of points
  !> npts, its volume vol (bohr^3) and the range rmin..rmax of the
  !> signed density in it. errmsg is non-empty in case of error.
  module subroutine nci_domains(opts,res,splot,ndom,lbl,npts,vol,rmin,rmax,errmsg)
    use grid3mod, only: grid_domains
    use tools_math, only: det3
    type(nci_options), intent(in) :: opts
    type(nci_result), intent(in) :: res
    real*8, intent(in) :: splot(:,:,:)
    integer, intent(out) :: ndom
    integer, allocatable, intent(out) :: lbl(:,:,:)
    integer, allocatable, intent(out) :: npts(:)
    real*8, allocatable, intent(out) :: vol(:), rmin(:), rmax(:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, j, k, id, n(3)

    n = res%n
    errmsg = ""
    ndom = 0
    if (any(shape(splot) /= n)) then
       errmsg = "the masked grid does not match the result"
       return
    end if
    call grid_domains(splot,opts%dimplot,opts%periodic,ndom,lbl,npts,errmsg)
    if (len_trim(errmsg) > 0) return

    allocate(vol(ndom),rmin(ndom),rmax(ndom))
    vol = real(npts,8) * abs(det3(res%xmat))
    rmin = huge(1d0)
    rmax = -huge(1d0)
    do k = 1, n(3)
       do j = 1, n(2)
          do i = 1, n(1)
             id = lbl(i,j,k)
             if (id == 0) cycle
             rmin(id) = min(rmin(id),res%srho(i,j,k))
             rmax(id) = max(rmax(id),res%srho(i,j,k))
          end do
       end do
    end do

  end subroutine nci_domains

  !> Write the NCI cube files for root oname: the signed density times
  !> 100 (oname-dens.cube), the reduced density gradient of the plot,
  !> splot (oname-grad.cube), and, if it was calculated, the
  !> promolecular density (oname-pdens.cube).
  module subroutine nci_write_cubes(oname,opts,s,res,splot,errmsg)
    use systemmod, only: system
    use tools_io, only: fopen_write, fclose
    character*(*), intent(in) :: oname
    type(nci_options), intent(in) :: opts
    type(system), intent(in) :: s
    type(nci_result), intent(in) :: res
    real*8, intent(in) :: splot(:,:,:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: lugc, ludc, lupc, nat
    logical :: isortho, ok
    real*8 :: rchol(3,3), delta(3)
    integer, allocatable :: iz(:)
    real*8, allocatable :: xat(:,:)

    errmsg = ""
    call grid_frame(res%xmat,isortho,rchol,delta,ok)
    if (.not.ok) then
       errmsg = "region too small"
       return
    end if
    call header_atoms(opts,s,nat,iz,xat)

    ! open the files
    lugc = fopen_write(trim(oname)//"-grad.cube",errstop=.false.)
    ludc = fopen_write(trim(oname)//"-dens.cube",errstop=.false.)
    lupc = 0
    if (allocated(res%rhoat)) lupc = fopen_write(trim(oname)//"-pdens.cube",errstop=.false.)
    if (lugc < 0 .or. ludc < 0 .or. lupc < 0) then
       errmsg = "error opening the cube files"
       if (lugc >= 0) call fclose(lugc)
       if (ludc >= 0) call fclose(ludc)
       if (lupc > 0) call fclose(lupc)
       return
    end if

    ! write the headers and the bodies (the body closes the file)
    call write_cube_header(lugc,'grad_cube','3d plot, reduced density gradient',res,isortho,nat,iz,xat)
    call write_cube_header(ludc,'dens_cube','3d plot, density',res,isortho,nat,iz,xat)
    if (lupc > 0) &
       call write_cube_header(lupc,'pdens_cube','3d plot, promolecular density',res,isortho,nat,iz,xat)
    call write_cube_body(ludc,res%n,res%srho,100d0)
    call write_cube_body(lugc,res%n,splot,1d0)
    if (lupc > 0) call write_cube_body(lupc,res%n,res%rhoat,1d0)

  end subroutine nci_write_cubes

  !> Write the s vs. rho points of the NCI grids res to file: the
  !> signed density and the reduced density gradient at the grid
  !> points that pass the fragment and void criteria (inter) and are
  !> under the dat cutoffs, in (i,j,k) order with k fastest.
  module subroutine nci_write_dat(file,opts,res,inter,errmsg)
    use tools_io, only: fopen_write, fclose
    character*(*), intent(in) :: file
    type(nci_options), intent(in) :: opts
    type(nci_result), intent(in) :: res
    logical(1), intent(in) :: inter(:,:,:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: lu, i, j, k, i0, ii, nb, n(3)
    real*8, allocatable :: bs(:,:,:), bg(:,:,:)
    logical(1), allocatable :: bl(:,:,:)

    errmsg = ""
    lu = fopen_write(file,errstop=.false.)
    if (lu < 0) then
       errmsg = "error opening the dat file"
       return
    end if

    ! transpose nblk i-planes at a time, so the grids are read in memory order
    n = res%n
    allocate(bs(n(3),n(2),nblk),bg(n(3),n(2),nblk),bl(n(3),n(2),nblk))
    do i0 = 1, n(1), nblk
       nb = min(nblk,n(1)-i0+1)
       do k = 1, n(3)
          do j = 1, n(2)
             do ii = 1, nb
                i = i0 + ii - 1
                bl(k,j,ii) = nci_in_dat(opts,res,inter,i,j,k)
                if (bl(k,j,ii)) then
                   bs(k,j,ii) = res%srho(i,j,k)
                   bg(k,j,ii) = res%rdg(i,j,k)
                end if
             end do
          end do
       end do
       do ii = 1, nb
          do j = 1, n(2)
             do k = 1, n(3)
                if (bl(k,j,ii)) write(lu,'(1p,E15.7,E15.7)') bs(k,j,ii), bg(k,j,ii)
             end do
          end do
       end do
    end do
    call fclose(lu)

  end subroutine nci_write_dat

  !> Write the VMD script (oname.vmd) and the xyz file with the atoms
  !> in the plot (oname_cell.xyz) for an NCI calculation with options
  !> opts in system s. The atoms are those of the fragments, of the
  !> unit cell (periodic grids; with domolmotif, the molecules are
  !> completed), or inside the region box. Both are written in the
  !> frame VMD uses for the cubes. The script names the xyz and cube
  !> files from oname, or from refroot if given (the files next to
  !> the script by their base name, say).
  module subroutine nci_write_vmd(oname,opts,s,domolmotif,errmsg,refroot)
    use systemmod, only: system
    use tools_io, only: fopen_write, fclose
    use types, only: realloc
    use param, only: bohrtoa, icrd_crys
    character*(*), intent(in) :: oname
    type(nci_options), intent(in) :: opts
    type(system), intent(inout) :: s
    logical, intent(in) :: domolmotif
    character(len=:), allocatable, intent(out) :: errmsg
    character*(*), intent(in), optional :: refroot

    character(len=:), allocatable :: rref

    type(fragment) :: fr0
    type(fragment), allocatable :: fr(:)
    logical, allocatable :: isdiscrete(:)
    integer :: luvmd, nmol, nat0, i, j, k, l, imin(3), imax(3)
    real*8 :: x(3), xx0(3), xx1(3), rchol(3,3), delta(3), xlo(3), xhi(3)
    real*8, allocatable :: fx(:,:)
    integer, allocatable :: fis(:), fidx(:), fcidx(:), flvec(:,:)
    logical :: isortho, ok

    real*8, parameter :: rthres_xyz = 1d0
    real*8, parameter :: xlist0(3,12) = reshape((/&
       0d0, 0d0, 0d0,&
       0d0, 0d0, 0d0,&
       0d0, 0d0, 0d0,&
       0d0, 0d0, 1d0,&
       0d0, 1d0, 1d0,&
       0d0, 1d0, 0d0,&
       1d0, 1d0, 0d0,&
       1d0, 0d0, 0d0,&
       1d0, 0d0, 1d0,&
       1d0, 1d0, 1d0,&
       1d0, 1d0, 1d0,&
       1d0, 1d0, 1d0/),shape(xlist0))
    real*8, parameter :: xlist1(3,12) = reshape((/&
       1d0, 0d0, 0d0,&
       0d0, 1d0, 0d0,&
       0d0, 0d0, 1d0,&
       0d0, 1d0, 1d0,&
       0d0, 1d0, 0d0,&
       1d0, 1d0, 0d0,&
       1d0, 0d0, 0d0,&
       1d0, 0d0, 1d0,&
       0d0, 0d0, 1d0,&
       0d0, 1d0, 1d0,&
       1d0, 0d0, 1d0,&
       1d0, 1d0, 0d0/),shape(xlist1))

    errmsg = ""

    ! compute the atomic positions for the xyz file
    if (opts%nfrag == 0) then
       if (opts%periodic) then
          fr0 = s%c%listatoms_cells((/1,1,1/),.true.)
          if (domolmotif) then
             call s%c%listmolecules(fr0,nmol,fr,isdiscrete)
             call fr0%merge_array(fr,.false.)
          end if
       else
          call grid_box(opts,xlo,xhi)
          xx0 = s%c%c2x(xlo - rthres_xyz)
          xx1 = s%c%c2x(xhi + rthres_xyz)
          imin = floor(min(xx0,xx1)) - 1
          imax = ceiling(max(xx0,xx1)) + 1
          ! the cell of a molecule is an artifact: no lattice images
          if (s%c%ismolecule) then
             imin = 0
             imax = 0
          end if

          ! gather the atoms inside the box
          nat0 = 0
          allocate(fx(3,10),fis(10),fidx(10),fcidx(10),flvec(3,10))
          do i = imin(1), imax(1)
             do j = imin(2), imax(2)
                do k = imin(3), imax(3)
                   do l = 1, s%c%ncel
                      x = s%c%x2c(s%c%atcel(l)%x + real((/i,j,k/),8))
                      if (all(x > xlo) .and. all(x < xhi)) then
                         nat0 = nat0 + 1
                         if (nat0 > size(fis)) then
                            call realloc(fx,3,2*nat0)
                            call realloc(fis,2*nat0)
                            call realloc(fidx,2*nat0)
                            call realloc(fcidx,2*nat0)
                            call realloc(flvec,3,2*nat0)
                         end if
                         fx(:,nat0) = s%c%atcel(l)%x + real((/i,j,k/),8)
                         fcidx(nat0) = l
                         fidx(nat0) = s%c%atcel(l)%idx
                         flvec(:,nat0) = (/i,j,k/)
                         fis(nat0) = s%c%atcel(l)%is
                      end if
                   end do
                end do
             end do
          end do
          call fr0%build(s%c%nspc,s%c%spc,nat0,fx(:,1:nat0),icrd_crys,fis(1:nat0),&
             fidx(1:nat0),fcidx(1:nat0),s%c%m_x2c,flvec(:,1:nat0))
          deallocate(fx,fis,fidx,fcidx,flvec)
       end if
    else
       call fr0%merge_array(opts%fr(1:opts%nfrag),.false.)
    end if
    if (fr0%nat == 0) then
       errmsg = "no atoms in region"
       return
    end if

    ! the frame VMD uses for the cubes
    call grid_frame(opts%xmat,isortho,rchol,delta,ok)
    if (.not.ok) then
       errmsg = "region too small"
       return
    end if

    ! Transform the atoms to agree with vmd
    do i = 1, fr0%nat
       fr0%at(i)%r = matmul(fr0%at(i)%r + delta,rchol)
    end do

    ! open the vmd script
    ! the files as the script names them: the root as given, or refroot
    rref = trim(oname)
    if (present(refroot)) rref = trim(refroot)
    luvmd = fopen_write(trim(oname)//".vmd",errstop=.false.)
    if (luvmd < 0) then
       errmsg = "error opening the vmd file"
       return
    end if

    ! write the fragment to the file
    call fr0%writexyz(trim(oname) // "_cell.xyz",.false.)

    ! write vmd script
    write (luvmd,'("#!/usr/local/bin/vmd")')
    write (luvmd,'("")')
    write (luvmd,'("# Display settings            ")')
    write (luvmd,'("display projection   Orthographic            ")')
    write (luvmd,'("display nearclip set 0.000000            ")')
    write (luvmd,'("display depthcue   off")')
    write (luvmd,'("light 2 on")')
    write (luvmd,'("light 3 on")')
    write (luvmd,'("axes location off")')
    write (luvmd,'("stage location off")')
    write (luvmd,'("color Display {Background} white")')
    write (luvmd,'("")')
    write (luvmd,'("# load the xyz with the atoms")')
    write (luvmd,'("mol new ",A," type xyz first 0 last -1 step 1 filebonds 0 autobonds 1 waitfor all")') rref//"_cell.xyz"
    write (luvmd,'("")')
    write (luvmd,'("# Render cell atoms as CPK")')
    write (luvmd,'("mol top 0")')
    write (luvmd,'("mol modstyle 0 0 CPK 1.000000 0.300000 50.000000 50.000000")')
    write (luvmd,'("mol modcolor 0 0 Element")')
    write (luvmd,'("mol modmaterial 0 0 Opaque")')
    write (luvmd,'("")')
    if (opts%nfrag == 0 .and. opts%periodic) then
       write (luvmd,'("# Cell")')
       write (luvmd,'("draw color red")')
       do i = 1, 12
          xx0 = matmul(s%c%x2c(xlist0(:,i)) + delta,rchol)
          xx1 = matmul(s%c%x2c(xlist1(:,i)) + delta,rchol)
          write (luvmd,'("draw cylinder {",3(F14.4," "),"} {",3(F14.4," "),"} radius 0.05 resolution 20")') &
             xx0 * bohrtoa, xx1 * bohrtoa
       end do
       write (luvmd,'("")')
    end if
    write (luvmd,'("# Load the cubes")')
    write (luvmd,'("mol new ",A," type cube first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all")')&
       rref//"-dens.cube"
    write (luvmd,'("mol addfile ",A," type cube first 0 last -1 step 1 filebonds 0 autobonds 0 waitfor all")')&
       rref//"-grad.cube"
    write (luvmd,'("")')
    write (luvmd,'("# representation of the atoms")')
    write (luvmd,'("mol top 1")')
    write (luvmd,'("mol modstyle 0 top Isosurface ",F8.2," 1 0 0 1 1")') opts%dimplot
    write (luvmd,'("mol modcolor 0 top Volume 0")')
    write (luvmd,'("mol modmaterial 0 top Opaque")')
    write (luvmd,'("mol scaleminmax top 0 ",F8.2," ",F8.2)') -opts%rhoplot/2d0*100d0, opts%rhoplot/2d0*100d0
    write (luvmd,'("mol selupdate 0 top 0")')
    write (luvmd,'("mol colupdate 0 top 0")')
    write (luvmd,'("color scale method BGR")')
    write (luvmd,'("")')
    write (luvmd,'("# center")')
    write (luvmd,'("mol top 0")')
    write (luvmd,'("display resetview")')
    write (luvmd,'("rotate x by 15")')
    write (luvmd,'("rotate y by 20")')
    write (luvmd,'("rotate z by 3")')
    call fclose(luvmd)

  end subroutine nci_write_vmd


  !> Write the files of an NCI calculation with options opts in
  !> system s, named from root: the s vs. density data (root.dat, the
  !> points inter from nci_masks), the cube files (root-dens.cube,
  !> root-grad.cube from the masked RDG splot, root-pdens.cube if the
  !> promolecular density was calculated), and the VMD script with
  !> its xyz file (root.vmd, root_cell.xyz; the script names the
  !> files from refroot if given, see nci_write_vmd). errmsg is
  !> non-empty if a file could not be written.
  module subroutine nci_write_files(root,opts,s,res,inter,splot,domolmotif,errmsg,refroot)
    use systemmod, only: system
    character*(*), intent(in) :: root
    type(nci_options), intent(in) :: opts
    type(system), intent(inout) :: s
    type(nci_result), intent(in) :: res
    logical(1), intent(in) :: inter(:,:,:)
    real*8, intent(in) :: splot(:,:,:)
    logical, intent(in) :: domolmotif
    character(len=:), allocatable, intent(out) :: errmsg
    character*(*), intent(in), optional :: refroot

    call nci_write_dat(trim(root) // ".dat",opts,res,inter,errmsg)
    if (len_trim(errmsg) > 0) return
    call nci_write_cubes(root,opts,s,res,splot,errmsg)
    if (len_trim(errmsg) > 0) return
    call nci_write_vmd(root,opts,s,domolmotif,errmsg,refroot)

  end subroutine nci_write_files

  !> Read the NCI grids from checkpoint file file into res; returns
  !> whether they were read. They are not if the file does not exist
  !> or cannot be read, or if it was written for another density
  !> (structure and field fingerprint, derivative method), grid, or
  !> set of fragments than options opts in system s, or lacks the
  !> promolecular grids opts needs. res is empty on failure.
  module function nci_chk_read(file,opts,s,res) result(ok)
    use systemmod, only: system
    use tools_io, only: fopen_read, fclose
    character*(*), intent(in) :: file
    type(nci_options), intent(in) :: opts
    type(system), intent(in) :: s
    type(nci_result), intent(out) :: res
    logical :: ok

    integer :: lu, ios, i, n(3), ni, nr, ns, nattot, nfrag, version
    integer, allocatable :: ihdr(:), ihdr0(:), iz(:), iz0(:), natfr(:), natfr0(:)
    real*8, allocatable :: rhdr(:), rhdr0(:), xat(:,:), xat0(:,:)
    character(len=:), allocatable :: str, str0
    character(len=8) :: magic
    logical :: llat, llfrag, lfourier, dopromol
    real*8 :: x0(3), xmat(3,3)

    ok = .false.
    inquire(file=file,exist=ok)
    if (.not.ok) return
    ok = .false.
    if (.not.s%goodfield(opts%ifield)) return
    n = opts%nstep
    nfrag = opts%nfrag
    dopromol = nci_dopromol(opts)

    ! open the file
    lu = fopen_read(file,"unformatted",errstop=.false.)
    if (lu < 0) return

    chk: block
      ! the format and version
      read (lu,iostat=ios) magic, version
      if (ios /= 0 .or. magic /= chk_magic .or. version /= chk_version) exit chk

      ! the structure and field fingerprint
      call s%f(opts%ifield)%chk_fingerprint(ihdr0,rhdr0,str0)
      read (lu,iostat=ios) ni, nr, ns
      if (ios /= 0 .or. ni /= size(ihdr0) .or. nr /= size(rhdr0) .or. ns /= len(str0)) exit chk
      allocate(ihdr(ni),rhdr(nr))
      allocate(character(len=ns) :: str)
      read (lu,iostat=ios) ihdr, rhdr, str
      if (ios /= 0) exit chk
      if (any(ihdr /= ihdr0) .or. any(abs(rhdr - rhdr0) > chk_eps) .or. str /= str0) exit chk

      ! the derivative method, the grids in the file, and the grid
      read (lu,iostat=ios) lfourier, llat, llfrag
      if (ios /= 0 .or. (lfourier .neqv. opts%fourierint)) exit chk
      if (dopromol .and. .not.llat) exit chk
      if (nfrag > 0 .and. .not.llfrag) exit chk
      read (lu,iostat=ios) x0, xmat, n
      if (ios /= 0 .or. any(abs(x0-opts%x0) > 1d-13) .or. any(abs(xmat-opts%xmat) > 1d-13) .or.&
         any(n /= opts%nstep)) exit chk

      ! the fragments
      read (lu,iostat=ios) i, nattot
      if (ios /= 0 .or. i /= nfrag) exit chk
      call fragment_atoms(opts,natfr0,iz0,xat0)
      if (nattot /= size(iz0)) exit chk
      allocate(natfr(nfrag),iz(nattot),xat(3,nattot))
      read (lu,iostat=ios) natfr, iz, xat
      if (ios /= 0 .or. any(natfr /= natfr0) .or. any(iz /= iz0) .or. any(abs(xat - xat0) > chk_eps)) exit chk

      ! read the grids
      res%ifield = opts%ifield
      res%fourierint = opts%fourierint
      res%x0 = opts%x0
      res%xmat = opts%xmat
      res%n = n
      allocate(res%srho(n(1),n(2),n(3)),res%rdg(n(1),n(2),n(3)))
      read (lu,iostat=ios) res%srho, res%rdg
      if (ios /= 0) exit chk
      if (llat) then
         allocate(res%rhoat(n(1),n(2),n(3)))
         read (lu,iostat=ios) res%rhoat
         if (ios /= 0) exit chk
      end if
      if (llfrag) then
         allocate(res%rhofrag(n(1),n(2),n(3),nfrag))
         read (lu,iostat=ios) res%rhofrag
         if (ios /= 0) exit chk
      end if
      ok = .true.
    end block chk
    call fclose(lu)

    ! leave no half-read result behind
    if (.not.ok) then
       res%ifield = -1
       res%n = 0
       if (allocated(res%srho)) deallocate(res%srho)
       if (allocated(res%rdg)) deallocate(res%rdg)
       if (allocated(res%rhoat)) deallocate(res%rhoat)
       if (allocated(res%rhofrag)) deallocate(res%rhofrag)
    end if

  end function nci_chk_read

  !> Write the NCI grids res, calculated with options opts in system
  !> s, to checkpoint file file. errmsg is non-empty on error.
  module subroutine nci_chk_write(file,opts,s,res,errmsg)
    use systemmod, only: system
    use tools_io, only: fopen_write, fclose
    character*(*), intent(in) :: file
    type(nci_options), intent(in) :: opts
    type(system), intent(in) :: s
    type(nci_result), intent(in) :: res
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: lu
    integer, allocatable :: ihdr(:), iz(:), natfr(:)
    real*8, allocatable :: rhdr(:), xat(:,:)
    character(len=:), allocatable :: str

    errmsg = ""
    if (.not.s%goodfield(res%ifield)) then
       errmsg = "the density field of the NCI result is not available"
       return
    end if
    lu = fopen_write(file,"unformatted",errstop=.false.)
    if (lu < 0) then
       errmsg = "could not open checkpoint file " // file
       return
    end if

    ! the format, the fingerprint of the density, the grid, and the fragments
    call s%f(res%ifield)%chk_fingerprint(ihdr,rhdr,str)
    call fragment_atoms(opts,natfr,iz,xat)
    write (lu) chk_magic, chk_version
    write (lu) size(ihdr), size(rhdr), len(str)
    write (lu) ihdr, rhdr, str
    write (lu) res%fourierint, allocated(res%rhoat), allocated(res%rhofrag)
    write (lu) res%x0, res%xmat, res%n
    write (lu) opts%nfrag, size(iz)
    write (lu) natfr, iz, xat

    ! the grids
    write (lu) res%srho, res%rdg
    if (allocated(res%rhoat)) write (lu) res%rhoat
    if (allocated(res%rhofrag)) write (lu) res%rhofrag
    call fclose(lu)

  end subroutine nci_chk_write

  !xx! private procedures

  !> Whether field ifield of system s has a core contribution (a
  !> grid field with core augmentation).
  function field_has_core(s,ifield) result(usecore)
    use systemmod, only: system
    type(system), intent(in) :: s
    integer, intent(in) :: ifield
    logical :: usecore

    usecore = s%f(ifield)%usecore
    if (usecore) usecore = any(s%f(ifield)%zpsp /= -1)

  end function field_has_core

  !> The atoms of the fragments in options opts as flat lists: the
  !> number of atoms of each fragment (natfr), and the atomic numbers
  !> (iz) and Cartesian coordinates (xat) of all of them in order
  !> (for the checkpoint fingerprint).
  subroutine fragment_atoms(opts,natfr,iz,xat)
    type(nci_options), intent(in) :: opts
    integer, allocatable, intent(out) :: natfr(:)
    integer, allocatable, intent(out) :: iz(:)
    real*8, allocatable, intent(out) :: xat(:,:)

    integer :: i, j, n, nattot

    allocate(natfr(opts%nfrag))
    nattot = 0
    do i = 1, opts%nfrag
       natfr(i) = opts%fr(i)%nat
       nattot = nattot + natfr(i)
    end do
    allocate(iz(nattot),xat(3,nattot))
    n = 0
    do i = 1, opts%nfrag
       do j = 1, opts%fr(i)%nat
          n = n + 1
          iz(n) = opts%fr(i)%spc(opts%fr(i)%at(j)%is)%z
          xat(:,n) = opts%fr(i)%at(j)%r
       end do
    end do

  end subroutine fragment_atoms

  !> The axis-aligned box around the grid of options opts (origin
  !> x0, nstep steps along each column of xmat): its lower corner xlo
  !> and upper corner xhi (Cartesian, cell frame). The atoms the
  !> output files list for a non-periodic grid are the ones inside it.
  subroutine grid_box(opts,xlo,xhi)
    type(nci_options), intent(in) :: opts
    real*8, intent(out) :: xlo(3), xhi(3)

    integer :: j
    real*8 :: v(3)

    xlo = opts%x0
    xhi = opts%x0
    do j = 1, 3
       v = real(opts%nstep(j),8) * opts%xmat(:,j)
       xlo = xlo + min(v,0d0)
       xhi = xhi + max(v,0d0)
    end do

  end subroutine grid_box

  !> The frame VMD uses for a cube with step vectors xmat (columns).
  !> VMD screw-up log, entry #1: vmd internally converts the
  !> coordinate system to a Cartesian matrix that is lower triangular;
  !> rchol transforms the coordinates to it (x_vmd = matmul(x+delta,rchol)).
  !> Entry #3: the origin of the isosurface box is translated by -1/2
  !> of a step in each direction; delta is that shift. isortho is
  !> whether the step vectors are orthogonal. ok is false if the
  !> region is too small.
  subroutine grid_frame(xmat,isortho,rchol,delta,ok)
    use tools_math, only: m_x2c_from_cellpar, matinv, cellpar_from_metric
    real*8, intent(in) :: xmat(3,3)
    logical, intent(out) :: isortho
    real*8, intent(out) :: rchol(3,3)
    real*8, intent(out) :: delta(3)
    logical, intent(out) :: ok

    real*8 :: gg(3,3), aal(3), bbl(3), xmati(3,3)

    ok = .false.
    isortho = .true.
    rchol = 0d0
    delta = 0d0
    gg = matmul(transpose(xmat),xmat)
    call cellpar_from_metric(gg,aal,bbl)
    if (any(abs(aal) < 1d-2)) return
    isortho = all(abs(bbl - 90d0) < 1d-4)
    rchol = m_x2c_from_cellpar(aal,bbl)
    xmati = xmat
    call matinv(xmati,3)
    rchol = transpose(matmul(rchol,xmati))
    delta = -0.5d0 * (xmat(:,1) + xmat(:,2) + xmat(:,3))
    ok = .true.

  end subroutine grid_frame

  !> The atoms listed in the headers of the NCI cube files for options
  !> opts in system s: atomic numbers iz and Cartesian coordinates
  !> xat of the nat atoms. They are the atoms of the unit cell
  !> (periodic grids), of the fragments, or inside the region box plus
  !> a small margin; the unit cell again if the box has none.
  subroutine header_atoms(opts,s,nat,iz,xat)
    use systemmod, only: system
    use types, only: realloc
    type(nci_options), intent(in) :: opts
    type(system), intent(in) :: s
    integer, intent(out) :: nat
    integer, allocatable, intent(out) :: iz(:)
    real*8, allocatable, intent(out) :: xat(:,:)

    integer :: i, j, i1, i2, i3, nmin(3), nmax(3)
    real*8 :: xx(3), x0x(3), x1x(3), xlo(3), xhi(3)

    real*8, parameter :: eps = 0.1d0

    nat = 0
    if (.not.opts%periodic) then
       if (opts%nfrag == 0) then
          call grid_box(opts,xlo,xhi)
          x0x = s%c%c2x(xlo)
          x1x = s%c%c2x(xhi)
          nmin = floor(x0x) - 1
          nmax = ceiling(x1x) + 1
          ! the cell of a molecule is an artifact: no lattice images
          if (s%c%ismolecule) then
             nmin = 0
             nmax = 0
          end if
          allocate(iz(s%c%ncel),xat(3,s%c%ncel))
          do i = 1, s%c%ncel
             do i1 = nmin(1), nmax(1)
                do i2 = nmin(2), nmax(2)
                   do i3 = nmin(3), nmax(3)
                      xx = s%c%x2c(s%c%atcel(i)%x + real((/i1,i2,i3/),8))
                      if (all(xx > xlo-eps) .and. all(xx < xhi+eps)) then
                         nat = nat + 1
                         if (nat > size(iz)) then
                            call realloc(iz,2*nat)
                            call realloc(xat,3,2*nat)
                         end if
                         iz(nat) = s%c%spc(s%c%atcel(i)%is)%z
                         xat(:,nat) = xx
                      end if
                   end do
                end do
             end do
          end do
       else
          do i = 1, opts%nfrag
             nat = nat + opts%fr(i)%nat
          end do
          allocate(iz(nat),xat(3,nat))
          nat = 0
          do i = 1, opts%nfrag
             do j = 1, opts%fr(i)%nat
                nat = nat + 1
                iz(nat) = s%c%spc(opts%fr(i)%at(j)%is)%z
                xat(:,nat) = opts%fr(i)%at(j)%r
             end do
          end do
       end if
    end if

    ! periodic grids, or an empty box: the unit cell
    if (nat == 0) then
       nat = s%c%ncel
       if (allocated(iz)) deallocate(iz)
       if (allocated(xat)) deallocate(xat)
       allocate(iz(nat),xat(3,nat))
       do i = 1, nat
          iz(i) = s%c%spc(s%c%atcel(i)%is)%z
          xat(:,i) = s%c%atcel(i)%r
       end do
    end if

  end subroutine header_atoms

  !> Write the header of an NCI cube file to unit lu with title lines
  !> l1 and l2, for the grid of result res and the nat atoms with
  !> atomic numbers iz and Cartesian coordinates xat. isortho is
  !> whether the grid step vectors are orthogonal (grid_frame).
  subroutine write_cube_header(lu,l1,l2,res,isortho,nat,iz,xat)
    integer, intent(in) :: lu
    character*(*), intent(in) :: l1, l2
    type(nci_result), intent(in) :: res
    logical, intent(in) :: isortho
    integer, intent(in) :: nat
    integer, intent(in) :: iz(nat)
    real*8, intent(in) :: xat(3,nat)

    integer :: i

    write(lu,*) trim(l1)
    write(lu,*) trim(l2)
    write(lu,'(I5,3(F12.6))') nat, res%x0
    write(lu,'(I5,3(F12.6))') res%n(1), res%xmat(:,1)

    ! VMD screw-up log, entry #2:
    ! vmd detects whether a cube is orthogonal or not by the (1,2) element
    ! of the Cartesian matrix only. If this value is zero, vmd will misalign
    ! the isosurface and its box, the latter being displayed correctly. To
    ! prevent this, nudge the (1,2) element a little bit if the cube is not
    ! orthogonal and that element happens to be zero.
    if (isortho .or. abs(res%xmat(1,2)) > 2.d-6) then
       write(lu,'(I5,3(F12.6))') res%n(2), res%xmat(:,2)
    else
       write(lu,'(I5,3(F12.6))') res%n(2), 1.01d-6, res%xmat(2:3,2)
    end if
    write(lu,'(I5,3(F12.6))') res%n(3), res%xmat(:,3)

    do i = 1, nat
       write(lu,'(I4,F5.1,F11.6,F11.6,F11.6)') iz(i), 0d0, xat(:,i)
    end do

  end subroutine write_cube_header

  !> Write the body of a cube file (grid c times scale, dimensions n,
  !> in grid order) to unit lu and close it. The cube lists the values
  !> in (i,j,k) order with k fastest: the grid is transposed through a
  !> buffer nblk i-planes at a time, so it is read in memory order.
  subroutine write_cube_body(lu,n,c,scale)
    use tools_io, only: fclose
    integer, intent(in) :: lu
    integer, intent(in) :: n(3)
    real*8, intent(in) :: c(n(1),n(2),n(3))
    real*8, intent(in) :: scale

    integer :: i0, ii, nb, j, k
    real*8, allocatable :: buf(:,:,:)

    allocate(buf(n(3),n(2),nblk))
    do i0 = 1, n(1), nblk
       nb = min(nblk,n(1)-i0+1)
       do k = 1, n(3)
          do j = 1, n(2)
             do ii = 1, nb
                buf(k,j,ii) = c(i0+ii-1,j,k) * scale
             end do
          end do
       end do
       do ii = 1, nb
          do j = 1, n(2)
             write (lu,'(6(" ",1p,e13.5e3))') (buf(k,j,ii),k=1,n(3))
          end do
       end do
    end do
    call fclose(lu)

  end subroutine write_cube_body

  function check_no_extra_word(line,lp,routine)
    use tools_io, only: getword, ferror, faterr
    character(len=:), allocatable :: aux2
    character*(*), intent(in) :: line, routine
    integer, intent(inout) :: lp
    logical :: check_no_extra_word

    check_no_extra_word = .true.
    aux2 = getword(line,lp)
    if (len_trim(aux2) > 0) then
       call ferror(trim(routine),'Unknown extra keyword',faterr,line,syntax=.true.)
       check_no_extra_word = .false.
    end if
  end function check_no_extra_word

  !> Read a FRAGMENT environment (atomic positions in angstrom, up to
  !> END/ENDFRAGMENT) from unit lu and build the fragment of crystal
  !> c with those atoms.
  function read_fragment(c,lu) result(fr)
    use crystalmod, only: crystal
    use global, only: eval_next
    use tools_io, only: uin, getline, lgetword, ucopy, equal, ferror, faterr
    use types, only: realloc
    use param, only: bohrtoa
    class(crystal), intent(inout) :: c
    integer, intent(in) :: lu
    type(fragment) :: fr

    character(len=:), allocatable :: line, word
    integer :: lp, n
    logical :: ok
    real*8 :: x(3)
    real*8, allocatable :: fx(:,:)

    ! gather the atoms given in the input
    n = 0
    allocate(fx(3,10))
    do while (.true.)
       if (lu == uin) then
          ok = getline(lu,line,.true.,ucopy)
       else
          ok = getline(lu,line,.true.)
       endif
       lp = 1
       word = lgetword(line,lp)
       if (equal(word,"end") .or. equal(word,"endfragment")) exit
       lp = 1
       ok = eval_next(x(1),line,lp)
       ok = ok .and. eval_next(x(2),line,lp)
       ok = ok .and. eval_next(x(3),line,lp)
       if (.not.ok) call ferror('nciplot','bad atom in fragment',faterr)
       n = n + 1
       if (n > size(fx,2)) call realloc(fx,3,2*n)
       fx(:,n) = x / bohrtoa - c%molx0
    end do

    ! identify the atoms in the crystal and build the fragment
    fr = c%identify_fragment(n,fx(:,1:n))

  end function read_fragment

  !> The checks of the promolecular part of an NCI calculation with
  !> options opts of system s: the grid, field 0, and the fragments.
  !> errmsg is non-empty if they fail.
  subroutine promol_checks(opts,s,errmsg)
    use systemmod, only: system
    use fieldmod, only: type_promol
    type(nci_options), intent(in) :: opts
    type(system), intent(in) :: s
    character(len=:), allocatable, intent(out) :: errmsg

    errmsg = ""
    if (any(opts%nstep < 1)) then
       errmsg = "invalid NCI grid dimensions"
    elseif (nci_dopromol(opts) .and. s%f(0)%type /= type_promol) then
       errmsg = "promolecular density not found"
    elseif (opts%nfrag > 0) then
       if (.not.allocated(opts%fr)) then
          errmsg = "the fragments are not allocated"
       elseif (size(opts%fr) < opts%nfrag) then
          errmsg = "the fragments are not allocated"
       end if
    end if

  end subroutine promol_checks

  !> At point x (cell-frame Cartesian bohr), grid node (i,j,k) of res:
  !> the promolecular density of system s (if dopromol, in res%rhoat)
  !> and of the first nfrag fragments of opts (in res%rhofrag): their
  !> atoms where they are or, with opts%frcell, their cell atoms and
  !> all the lattice translations.
  subroutine promol_point(opts,s,x,dopromol,nfrag,res,i,j,k)
    use systemmod, only: system
    use grid1mod, only: agrid
    use param, only: icrd_cart
    type(nci_options), intent(in) :: opts
    type(system), intent(inout) :: s
    real*8, intent(in) :: x(3)
    logical, intent(in) :: dopromol
    integer, intent(in) :: nfrag
    type(nci_result), intent(inout) :: res
    integer, intent(in) :: i, j, k

    integer :: ifr, iat, iz
    real*8 :: xd(3), dist, rrho, rrho1, rrho2, rp(3), rpp(3,3)

    if (dopromol) res%rhoat(i,j,k) = s%f(0)%grd0(x)
    if (opts%frcell) then
       ! (the cell atoms of each fragment and their lattice
       ! translations, like the promolecular density)
       do ifr = 1, nfrag
          call s%c%promolecular_atom(x,icrd_cart,res%rhofrag(i,j,k,ifr),rp,rpp,0,fr=opts%fr(ifr))
       end do
       return
    end if
    do ifr = 1, nfrag
       res%rhofrag(i,j,k,ifr) = 0d0
       do iat = 1, opts%fr(ifr)%nat
          xd = x - opts%fr(ifr)%at(iat)%r
          dist = norm2(xd)
          iz = opts%fr(ifr)%spc(opts%fr(ifr)%at(iat)%is)%z
          if (.not.agrid(iz)%isinit) cycle
          if (dist > agrid(iz)%rmax) cycle
          dist = max(dist,agrid(iz)%r(1))
          dist = max(dist,1d-14)
          call agrid(iz)%interp(dist,rrho,rrho1,rrho2)
          rrho = max(rrho,0d0)
          res%rhofrag(i,j,k,ifr) = res%rhofrag(i,j,k,ifr) + rrho
       end do
    end do

  end subroutine promol_point

end submodule proc
