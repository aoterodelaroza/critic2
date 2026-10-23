! Copyright (c) 2016-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! FHI-aims wavefunctions, read from the dgrid_aims.dat file (output dgrid).
submodule (fhiaims_private) proc
  implicit none

  !xx! private procedures
  ! subroutine check_kgrid(nkxyz,kpt,errmsg)
  ! subroutine map_atoms(f,c,lat,errmsg)
  ! subroutine build_blocks(f,c)
  ! subroutine build_dm(f,evc,occ,kpt,wk,nsv)
  ! function find_block(f,a,b,dt,trans) result(ib)
  ! subroutine radial_index(s,r,i,dt,tp)
  ! subroutine radial_eval(s,ifn,r,i,dt,tp,nder,u,up,upp)
  ! subroutine contract(d,nr,nc,both,nder,phir,phic,gr,gc,wr,wc,wur,wuc)

  ! minimum distance (bohr)
  real*8, parameter :: mindist = 1d-6
  ! occupations below this are not stored
  real*8, parameter :: occ_cutoff = 1d-14
  ! tolerance for comparing lattice vectors (bohr)
  real*8, parameter :: lat_eps = 1d-3

contains

  !> Deallocate data
  module subroutine fhiaims_end(f)
    use iso_c_binding, only: c_null_ptr
    class(fhiaimswfn), intent(inout) :: f

    if (allocated(f%spc)) deallocate(f%spc)
    if (allocated(f%ispc)) deallocate(f%ispc)
    if (allocated(f%xat)) deallocate(f%xat)
    if (allocated(f%ibas0)) deallocate(f%ibas0)
    if (allocated(f%iaims)) deallocate(f%iaims)
    if (allocated(f%ioff)) deallocate(f%ioff)
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    if (allocated(f%pair0)) deallocate(f%pair0)
    if (allocated(f%pairn)) deallocate(f%pairn)
    if (allocated(f%blkt)) deallocate(f%blkt)
    if (allocated(f%blkoff)) deallocate(f%blkoff)
    if (allocated(f%dm)) deallocate(f%dm)
    f%nblk = 0
    f%cptr = c_null_ptr

  end subroutine fhiaims_end

  !> Read an FHI-aims wavefunction from the dgrid_aims.dat file
  !> (output dgrid keyword in control.in), map its atoms to the
  !> structure in cptr, and build the real-space density matrix.
  module subroutine fhiaims_read(f,cptr,file,errmsg,ti)
    use crystalmod, only: crystal
    use iso_c_binding, only: c_f_pointer
    use tools_io, only: fopen_read, fclose, getline_raw
    class(fhiaimswfn), intent(inout) :: f !< Output field
    type(c_ptr), intent(in) :: cptr !< Crystal structure
    character*(*), intent(in) :: file !< The dgrid_aims.dat file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, k, is, nper, nfnall, idum(3), nkxyz(3), iflag, nsv
    character(len=:), allocatable :: line
    character*20 :: name
    integer, allocatable :: fsp(:), fl(:), nfnsp(:), ifn(:)
    real*8, allocatable :: frc(:), occ(:,:,:), xr(:), kpt(:,:), wk(:)
    complex*16, allocatable :: evc(:,:,:,:)
    real*8 :: rdum, lat(3,3), rinc
    type(crystal), pointer :: c

    call c_f_pointer(cptr,c)
    f%cptr = cptr
    errmsg = "Error reading file: " // trim(file)
    lu = fopen_read(file,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if

    ! header and dimensions
    if (.not.getline_raw(lu,line)) goto 999
    if (index(line,"FHI-aims - DGrid Interface") == 0) then
       errmsg = "Not a dgrid_aims.dat file (FHI-aims output dgrid): " // trim(file)
       goto 999
    end if
    read (lu,*,err=999,end=999) f%nspc, f%nat, idum(1), f%nstates, f%nbas, nfnall, nper, f%nspin
    read (lu,*,err=999,end=999) f%nkpt, nkxyz, idum(1)
    if (nper /= 0 .and. nper /= 3) then
       errmsg = "Only molecules and 3D-periodic FHI-aims calculations are supported: " // trim(file)
       goto 999
    end if
    f%ismol = (nper == 0)
    allocate(f%spc(f%nspc))
    read (lu,*,err=999,end=999) (f%spc(i)%ngrid, i=1,f%nspc)

    ! atoms
    allocate(f%ispc(f%nat),f%xat(3,f%nat))
    do i = 1, f%nat
       read (lu,*,err=999,end=999) f%ispc(i), f%xat(:,i), name
       if (f%ispc(i) < 1 .or. f%ispc(i) > f%nspc) goto 999
       f%spc(f%ispc(i))%name = name
    end do

    ! total energy and chemical potential
    if (.not.getline_raw(lu,line)) goto 999
    if (index(line,"Total energy") == 0) goto 999
    if (.not.getline_raw(lu,line)) goto 999
    if (index(line,"Chemical pot") == 0) goto 999

    ! radial functions
    allocate(fsp(nfnall),fl(nfnall),frc(nfnall))
    do j = 1, nfnall
       read (lu,*,err=999,end=999) fsp(j), fl(j), idum(1), frc(j)
       if (fsp(j) < 1 .or. fsp(j) > f%nspc) goto 999
    end do
    allocate(nfnsp(f%nspc))
    read (lu,*,err=999,end=999) nfnsp
    do is = 1, f%nspc
       f%spc(is)%nfn = count(fsp == is)
       if (f%spc(is)%nfn /= nfnsp(is)) then
          errmsg = "Inconsistent number of radial functions in: " // trim(file)
          goto 999
       end if
       allocate(f%spc(is)%fn(f%spc(is)%nfn))
    end do

    ! radial grids
    do is = 1, f%nspc
       read (lu,*,err=999,end=999) f%spc(is)%rmin, rinc
       f%spc(is)%linc = log(rinc)
    end do

    ! spline coefficients of the radial functions
    allocate(ifn(f%nspc))
    ifn = 0
    do j = 1, nfnall
       is = fsp(j)
       ifn(is) = ifn(is) + 1
       f%spc(is)%fn(ifn(is))%l = fl(j)
       f%spc(is)%fn(ifn(is))%rcut = frc(j)
       allocate(f%spc(is)%fn(ifn(is))%spl(4,f%spc(is)%ngrid))
       read (lu,*,err=999,end=999) f%spc(is)%fn(ifn(is))%spl
    end do
    f%lmax = 0
    do is = 1, f%nspc
       f%spc(is)%nbas = sum(2 * f%spc(is)%fn(:)%l + 1)
       f%spc(is)%rcut = maxval(f%spc(is)%fn(:)%rcut)
       f%spc(is)%lmax = maxval(f%spc(is)%fn(:)%l)
       f%lmax = max(f%lmax,f%spc(is)%lmax)
    end do

    ! first basis function of each atom
    allocate(f%ibas0(f%nat))
    k = 0
    do i = 1, f%nat
       f%ibas0(i) = k
       k = k + f%spc(f%ispc(i))%nbas
    end do
    if (k /= f%nbas) then
       errmsg = "Inconsistent number of basis functions in: " // trim(file)
       goto 999
    end if

    ! eigenvalues and occupations
    allocate(occ(f%nstates,f%nspin,f%nkpt))
    do k = 1, f%nkpt
       do is = 1, f%nspin
          do i = 1, f%nstates
             read (lu,*,err=999,end=999) rdum, occ(i,is,k)
          end do
       end do
    end do

    ! only the states up to the last occupied one are kept
    nsv = 0
    do i = 1, f%nstates
       if (any(abs(occ(i,:,:)) > occ_cutoff)) nsv = i
    end do

    ! eigenvectors (real ones are stored as complex with zero imaginary part)
    read (lu,*,err=999,end=999) iflag
    f%isreal = (iflag < 0)
    allocate(evc(f%nbas,nsv,f%nspin,f%nkpt))
    if (f%isreal) then
       allocate(xr(f%nbas))
    else
       allocate(xr(2*f%nbas))
    end if
    do k = 1, f%nkpt
       do is = 1, f%nspin
          do i = 1, f%nstates
             read (lu,*,err=999,end=999) xr
             if (i > nsv) cycle
             if (f%isreal) then
                evc(:,i,is,k) = xr
             else
                evc(:,i,is,k) = cmplx(xr(1::2),xr(2::2),8)
             end if
          end do
       end do
    end do

    ! lattice vectors and k-points
    allocate(kpt(3,f%nkpt),wk(f%nkpt))
    if (f%ismol) then
       lat = 0d0
       kpt = 0d0
       wk = 1d0
    else
       do i = 1, 3
          read (lu,*,err=999,end=999) lat(:,i)
       end do
       do i = 1, 3
          read (lu,*,err=999,end=999) rdum, rdum, rdum
       end do
       do k = 1, f%nkpt
          read (lu,*,err=999,end=999) kpt(:,k), wk(k)
       end do
    end if
    call fclose(lu)

    ! the k-points must cover the full grid under time reversal
    if (.not.f%ismol) then
       call check_kgrid(nkxyz,kpt,errmsg)
       if (len_trim(errmsg) > 0) return
    end if

    ! map the atoms, set up the density matrix blocks, and build them
    call map_atoms(f,c,lat,errmsg)
    if (len_trim(errmsg) > 0) return
    call build_blocks(f,c)
    call build_dm(f,evc,occ,kpt,wk,nsv)

    errmsg = ""
    return
999 continue
    call fclose(lu)

  end subroutine fhiaims_read

  !> Calculate the density and derivatives of an FHI-aims field (f)
  !> up to the nder degree (max = 2). xpos is in Cartesian
  !> coordinates. In output, the density (rho), the gradient (grad),
  !> the Hessian (h), and the G(r) kinetic energy density (gkin).
  !> This routine is thread-safe.
  module subroutine rho2(f,xpos,nder,rho,grad,h,gkin)
    use crystalmod, only: crystal
    use iso_c_binding, only: c_f_pointer
    use tools_math, only: tosphere, genylm, ylmderiv
    use param, only: icrd_cart
    class(fhiaimswfn), intent(inout) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    integer, intent(in) :: nder  !< Number of derivatives
    real*8, intent(out) :: rho !< Density
    real*8, intent(out) :: grad(3) !< Gradient
    real*8, intent(out) :: h(3,3) !< Hessian
    real*8, intent(out) :: gkin !< G(r), kinetic energy density

    integer :: nenv, ion, nenvl, ia, is, ifn, l, m, iy, ip, im, iloc, ispl
    integer :: nloc, p, q, ib, i, dt(3)
    logical :: trans
    real*8 :: xion(3), r, tp(2), u, up, upp, rr, rp, rpp, dspl, tspl, iphas
    integer, allocatable :: eid(:), lvec(:,:), ienv(:), iat(:), tvec(:,:), loc0(:)
    real*8, allocatable :: dist(:), phi(:), phip(:,:), phipp(:,:), w(:), wu(:,:)
    complex*16, allocatable :: ylm(:)
    complex*16 :: xgrad(3), xhess(6)
    type(crystal), pointer :: c

    real*8, parameter :: sqrt2 = sqrt(2d0)

    rho = 0d0
    grad = 0d0
    h = 0d0
    gkin = 0d0

    ! atoms around the input point; keep the ones within their own cutoff
    call c_f_pointer(f%cptr,c)
    call c%list_near_atoms(xpos,icrd_cart,.false.,nenv,eid,dist,lvec,up2dsp=f%spcutoff)
    if (nenv == 0) return
    allocate(ienv(nenv),iat(nenv),tvec(3,nenv),loc0(nenv+1))
    nenvl = 0
    nloc = 0
    do ion = 1, nenv
       ia = f%iaims(eid(ion))
       if (dist(ion) > f%spc(f%ispc(ia))%rcut) cycle
       nenvl = nenvl + 1
       ienv(nenvl) = ion
       iat(nenvl) = ia
       tvec(:,nenvl) = lvec(:,ion) - f%ioff(:,eid(ion))
       loc0(nenvl) = nloc
       nloc = nloc + f%spc(f%ispc(ia))%nbas
    end do
    loc0(nenvl+1) = nloc
    if (nenvl == 0) return

    ! values and derivatives of the local basis functions
    allocate(phi(nloc),ylm((f%lmax+3)**2))
    phi = 0d0
    if (nder >= 1) then
       allocate(phip(3,nloc))
       phip = 0d0
    end if
    if (nder >= 2) then
       allocate(phipp(6,nloc))
       phipp = 0d0
    end if
    do p = 1, nenvl
       ion = ienv(p)
       is = f%ispc(iat(p))
       xion = xpos - (c%atcel(eid(ion))%r + c%x2c(real(lvec(:,ion),8)))
       call tosphere(xion,r,tp)
       r = max(r,mindist)
       call genylm(f%spc(is)%lmax+2,tp,ylm)
       call radial_index(f%spc(is),r,ispl,dspl,tspl)

       iloc = loc0(p)
       do ifn = 1, f%spc(is)%nfn
          l = f%spc(is)%fn(ifn)%l
          if (r > f%spc(is)%fn(ifn)%rcut) then
             iloc = iloc + 2*l+1
             cycle
          end if

          ! radial part, R = u/r
          call radial_eval(f%spc(is),ifn,r,ispl,dspl,tspl,nder,u,up,upp)
          rr = u / r
          rp = up / r - u / r**2
          rpp = upp / r - 2d0 * up / r**2 + 2d0 * u / r**3

          ! real spherical harmonics in the FHI-aims convention
          ! (SHEval): S_l0 = Y_l0, S_lm = sqrt2 Re(Y_lm) and S_l-m =
          ! (-1)^m sqrt2 Im(Y_lm) for m > 0, with Condon-Shortley Y_lm.
          ! Local basis order is m = -l..l.
          iy = l*l + l + 1
          ip = iloc + l + 1
          call ylmderiv(ylm,r,l,0,rr,rp,rpp,nder,xgrad,xhess)
          phi(ip) = rr * real(ylm(iy),8)
          if (nder >= 1) phip(:,ip) = real(xgrad,8)
          if (nder >= 2) phipp(:,ip) = real(xhess,8)
          iphas = 1d0
          do m = 1, l
             iphas = -iphas
             ip = iloc + l + 1 + m
             im = iloc + l + 1 - m
             call ylmderiv(ylm,r,l,m,rr,rp,rpp,nder,xgrad,xhess)
             phi(ip) = sqrt2 * rr * real(ylm(iy+m),8)
             phi(im) = sqrt2 * iphas * rr * aimag(ylm(iy+m))
             if (nder < 1) cycle
             phip(:,ip) = sqrt2 * real(xgrad,8)
             phip(:,im) = sqrt2 * iphas * aimag(xgrad)
             if (nder < 2) cycle
             phipp(:,ip) = sqrt2 * real(xhess,8)
             phipp(:,im) = sqrt2 * iphas * aimag(xhess)
          end do
          iloc = iloc + 2*l+1
       end do
    end do

    ! contract with the density matrix: w = D * phi, wu = D * grad(phi)
    allocate(w(nloc))
    w = 0d0
    if (nder >= 1) then
       allocate(wu(3,nloc))
       wu = 0d0
    else
       allocate(wu(3,1),phip(3,1))
    end if
    do p = 1, nenvl
       do q = p, nenvl
          dt = tvec(:,p) - tvec(:,q)
          ib = find_block(f,iat(p),iat(q),dt,trans)
          if (ib == 0) cycle
          if (trans) then
             call contract(f%dm(f%blkoff(ib)+1),loc0(q+1)-loc0(q),loc0(p+1)-loc0(p),p /= q,nder,&
                phi(loc0(q)+1:),phi(loc0(p)+1:),phip(:,loc0(q)+1:),phip(:,loc0(p)+1:),&
                w(loc0(q)+1:),w(loc0(p)+1:),wu(:,loc0(q)+1:),wu(:,loc0(p)+1:))
          else
             call contract(f%dm(f%blkoff(ib)+1),loc0(p+1)-loc0(p),loc0(q+1)-loc0(q),p /= q,nder,&
                phi(loc0(p)+1:),phi(loc0(q)+1:),phip(:,loc0(p)+1:),phip(:,loc0(q)+1:),&
                w(loc0(p)+1:),w(loc0(q)+1:),wu(:,loc0(p)+1:),wu(:,loc0(q)+1:))
          end if
       end do
    end do

    ! accumulate properties
    rho = sum(phi * w)
    if (nder >= 1) then
       do i = 1, nloc
          grad = grad + 2d0 * phip(:,i) * w(i)
          gkin = gkin + 0.5d0 * dot_product(phip(:,i),wu(:,i))
       end do
    end if
    if (nder >= 2) then
       do i = 1, nloc
          h(1,1) = h(1,1) + 2d0 * (phipp(1,i) * w(i) + phip(1,i) * wu(1,i))
          h(1,2) = h(1,2) + 2d0 * (phipp(2,i) * w(i) + phip(1,i) * wu(2,i))
          h(1,3) = h(1,3) + 2d0 * (phipp(3,i) * w(i) + phip(1,i) * wu(3,i))
          h(2,2) = h(2,2) + 2d0 * (phipp(4,i) * w(i) + phip(2,i) * wu(2,i))
          h(2,3) = h(2,3) + 2d0 * (phipp(5,i) * w(i) + phip(2,i) * wu(3,i))
          h(3,3) = h(3,3) + 2d0 * (phipp(6,i) * w(i) + phip(3,i) * wu(3,i))
       end do
       h(2,1) = h(1,2)
       h(3,1) = h(1,3)
       h(3,2) = h(2,3)
    end if

  end subroutine rho2

  !xx! private procedures

  !> Check that the k-point list (fractional coordinates) is the
  !> time-reversal reduction of the full nkxyz grid: expanding every k
  !> to {k,-k} modulo G must give product(nkxyz) distinct points. A
  !> space-group-reduced list (symmetry_reduced_k_grid_spg) fails this
  !> and would give the unsymmetrized density. Skipped if the header
  !> does not give the grid (e.g. external k-point lists).
  subroutine check_kgrid(nkxyz,kpt,errmsg)
    use tools_io, only: string
    integer, intent(in) :: nkxyz(3)
    real*8, intent(in) :: kpt(:,:)
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: nk, n, i, j, s
    integer, allocatable :: ik(:,:)
    logical :: found

    integer, parameter :: kscale = 1000000

    errmsg = ""
    if (any(nkxyz <= 0)) return
    nk = size(kpt,2)
    allocate(ik(3,2*nk))
    n = 0
    do i = 1, nk
       do s = 1, -1, -2
          ik(:,n+1) = modulo(nint(s * kpt(:,i) * kscale),kscale)
          found = .false.
          do j = 1, n
             if (all(ik(:,j) == ik(:,n+1))) then
                found = .true.
                exit
             end if
          end do
          if (.not.found) n = n + 1
       end do
    end do
    if (n /= product(nkxyz)) then
       errmsg = "The " // string(nk) // " k-points do not cover the " // string(nkxyz(1)) // "x" //&
          string(nkxyz(2)) // "x" // string(nkxyz(3)) // " grid under time reversal (" //&
          string(n) // " points); space-group-reduced k-grids are not supported, rerun FHI-aims &
          &without symmetry_reduced_k_grid_spg"
    end if

  end subroutine check_kgrid

  !> Map the FHI-aims atoms in f to the cell atoms of crystal c. The
  !> two structures must have the same Cartesian frame: the lattice
  !> vectors (lat) must coincide for crystals, and molecules may only
  !> be translated (critic2 recenters them by c%molx0).
  subroutine map_atoms(f,c,lat,errmsg)
    use crystalmod, only: crystal
    use tools_io, only: zatguess
    use param, only: icrd_cart
    type(fhiaimswfn), intent(inout) :: f
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: lat(3,3)
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: i, ia, iz, lv(3)

    errmsg = ""
    if (c%ncel /= f%nat) then
       errmsg = "The number of atoms in the FHI-aims file does not match the structure &
          &(ghost/empty atoms are not supported)"
       return
    end if
    if (c%ismolecule .neqv. f%ismol) then
       errmsg = "The periodicity of the FHI-aims calculation does not match the structure"
       return
    end if
    if (.not.f%ismol) then
       if (any(abs(c%m_x2c - lat) > lat_eps)) then
          errmsg = "The lattice vectors of the FHI-aims calculation do not match those of the &
             &structure; load the structure from geometry.in, aims.out, or the same frame"
          return
       end if
    end if

    ! each FHI-aims atom is a cell atom, possibly in another cell
    if (allocated(f%iaims)) deallocate(f%iaims)
    if (allocated(f%ioff)) deallocate(f%ioff)
    allocate(f%iaims(c%ncel),f%ioff(3,c%ncel))
    f%iaims = 0
    f%ioff = 0
    do ia = 1, f%nat
       i = c%identify_atom(f%xat(:,ia) - c%molx0,icrd_cart,lvec=lv)
       if (i == 0) then
          errmsg = "Could not match the atoms of the FHI-aims calculation to the structure; &
             &load the structure from geometry.in, aims.out, or the same frame"
          return
       elseif (f%iaims(i) /= 0) then
          errmsg = "Two FHI-aims atoms match the same atom of the structure"
          return
       end if
       iz = zatguess(f%spc(f%ispc(ia))%name)
       if (iz > 0 .and. iz /= c%spc(c%atcel(i)%is)%z) then
          errmsg = "The atomic species of the FHI-aims calculation do not match the structure"
          return
       end if
       f%iaims(i) = ia
       f%ioff(:,i) = lv
    end do

    ! cutoffs for the atomic environments
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    allocate(f%spcutoff(c%nspc,2))
    f%spcutoff = 0d0
    do i = 1, c%ncel
       f%spcutoff(c%atcel(i)%is,2) = max(f%spcutoff(c%atcel(i)%is,2),&
          f%spc(f%ispc(f%iaims(i)))%rcut)
    end do

  end subroutine map_atoms

  !> Determine the density matrix blocks: atom pairs (a,b) and lattice
  !> vectors dT = T_a - T_b such that the basis functions on a+T_a and
  !> on b+T_b overlap. Only a <= b is stored (for a == b, dT
  !> lexicographically non-negative); the block (b,a,-dT) is the
  !> transpose.
  subroutine build_blocks(f,c)
    use crystalmod, only: crystal
    type(fhiaimswfn), intent(inout) :: f
    type(crystal), intent(in) :: c

    integer :: a, b, i, j, k, n, nbox, tmax(3), mblk
    integer*8 :: ntot
    real*8 :: rmax, dmax, v(3)
    integer, allocatable :: tbox(:,:)
    real*8, allocatable :: lbox(:,:)

    ! lattice vectors that can connect overlapping atoms: |dT_i| <=
    ! |row i of m_c2x| * (2 rcut + max interatomic distance)
    tmax = 0
    if (.not.f%ismol) then
       rmax = 2d0 * maxval(f%spc(:)%rcut)
       dmax = 0d0
       do a = 1, f%nat
          do b = a+1, f%nat
             dmax = max(dmax,norm2(f%xat(:,a)-f%xat(:,b)))
          end do
       end do
       do i = 1, 3
          tmax(i) = ceiling(c%ar(i) * (rmax + dmax))
       end do
    end if
    nbox = product(2*tmax+1)
    allocate(tbox(3,nbox),lbox(3,nbox))
    n = 0
    do k = -tmax(3), tmax(3)
       do j = -tmax(2), tmax(2)
          do i = -tmax(1), tmax(1)
             n = n + 1
             tbox(:,n) = (/i,j,k/)
             lbox(:,n) = c%x2c(real(tbox(:,n),8))
          end do
       end do
    end do

    ! find the blocks; grow the block list as needed
    if (allocated(f%pair0)) deallocate(f%pair0)
    if (allocated(f%pairn)) deallocate(f%pairn)
    if (allocated(f%blkt)) deallocate(f%blkt)
    if (allocated(f%blkoff)) deallocate(f%blkoff)
    allocate(f%pair0(f%nat,f%nat),f%pairn(f%nat,f%nat))
    f%pair0 = 0
    f%pairn = 0
    mblk = 4 * f%nat * nbox
    allocate(f%blkt(3,mblk),f%blkoff(mblk))
    f%nblk = 0
    ntot = 0
    do a = 1, f%nat
       do b = a, f%nat
          f%pair0(a,b) = f%nblk + 1
          do n = 1, nbox
             if (a == b .and. lexneg(tbox(:,n))) cycle
             v = f%xat(:,a) - f%xat(:,b) + lbox(:,n)
             if (norm2(v) > f%spc(f%ispc(a))%rcut + f%spc(f%ispc(b))%rcut) cycle
             f%nblk = f%nblk + 1
             if (f%nblk > mblk) then
                mblk = 2 * mblk
                call realloc_blocks(mblk)
             end if
             f%pairn(a,b) = f%pairn(a,b) + 1
             f%blkt(:,f%nblk) = tbox(:,n)
             f%blkoff(f%nblk) = ntot
             ntot = ntot + f%spc(f%ispc(a))%nbas * f%spc(f%ispc(b))%nbas
          end do
       end do
    end do
    call realloc_blocks(f%nblk)
    if (allocated(f%dm)) deallocate(f%dm)
    allocate(f%dm(ntot))
    f%dm = 0d0

  contains
    subroutine realloc_blocks(nnew)
      integer, intent(in) :: nnew
      integer, allocatable :: it(:,:)
      integer*8, allocatable :: io(:)

      allocate(it(3,nnew),io(nnew))
      it(:,1:min(nnew,size(f%blkt,2))) = f%blkt(:,1:min(nnew,size(f%blkt,2)))
      io(1:min(nnew,size(f%blkoff))) = f%blkoff(1:min(nnew,size(f%blkoff)))
      call move_alloc(it,f%blkt)
      call move_alloc(io,f%blkoff)
    end subroutine realloc_blocks
  end subroutine build_blocks

  !> True if the first non-zero component of the integer vector is
  !> negative.
  pure function lexneg(t)
    integer, intent(in) :: t(3)
    logical :: lexneg

    integer :: i

    lexneg = .false.
    do i = 1, 3
       if (t(i) /= 0) then
          lexneg = (t(i) < 0)
          return
       end if
    end do

  end function lexneg

  !> Find the block for atoms a, b and lattice vector dT = T_a - T_b.
  !> Returns 0 if the block does not exist. trans is true if the
  !> stored block is (b,a,-dT), the transpose of the requested one.
  function find_block(f,a,b,dt,trans) result(ib)
    type(fhiaimswfn), intent(in) :: f
    integer, intent(in) :: a, b, dt(3)
    logical, intent(out) :: trans
    integer :: ib

    integer :: aa, bb, tt(3), i

    trans = (a > b) .or. (a == b .and. lexneg(dt))
    if (trans) then
       aa = b
       bb = a
       tt = -dt
    else
       aa = a
       bb = b
       tt = dt
    end if
    do i = f%pair0(aa,bb), f%pair0(aa,bb) + f%pairn(aa,bb) - 1
       if (all(f%blkt(:,i) == tt)) then
          ib = i
          return
       end if
    end do
    ib = 0

  end function find_block

  !> Build the real-space density matrix blocks from the eigenvectors
  !> (evc, nsv states kept), occupations (occ, not k-weighted),
  !> fractional k-points (kpt) and k-point weights (wk):
  !>   D_ij(dT) = Re sum_{s,k} wk sum_n occ c_in(k) c_jn(k)^* exp(2 pi i k.dT)
  !> The real part is exact for the time-reversal-reduced k-point list.
  subroutine build_dm(f,evc,occ,kpt,wk,nsv)
    use param, only: tpi, img
    type(fhiaimswfn), intent(inout) :: f
    complex*16, intent(in) :: evc(:,:,:,:)
    real*8, intent(in) :: occ(:,:,:)
    real*8, intent(in) :: kpt(:,:)
    real*8, intent(in) :: wk(:)
    integer, intent(in) :: nsv

    integer :: a, b, na, nb, a0, b0, is, k, i, j, n, ib
    integer*8 :: off
    complex*16, allocatable :: fw(:,:), ct(:,:), pa(:,:)
    complex*16 :: phase

    if (nsv == 0) return
    allocate(fw(f%nbas,nsv),ct(nsv,f%nbas))
    do is = 1, f%nspin
       do k = 1, f%nkpt
          ! occupation- and weight-scaled eigenvectors, and the conjugate transpose
          do n = 1, nsv
             fw(:,n) = evc(:,n,is,k) * (occ(n,is,k) * wk(k))
          end do
          ct = conjg(transpose(evc(:,1:nsv,is,k)))

          do a = 1, f%nat
             if (all(f%pairn(a,a:f%nat) == 0)) cycle
             na = f%spc(f%ispc(a))%nbas
             a0 = f%ibas0(a)
             ! the a-rows of the k-space density matrix, all columns
             pa = matmul(fw(a0+1:a0+na,:),ct)
             do b = a, f%nat
                nb = f%spc(f%ispc(b))%nbas
                b0 = f%ibas0(b)
                do ib = f%pair0(a,b), f%pair0(a,b) + f%pairn(a,b) - 1
                   phase = exp(img * tpi * dot_product(kpt(:,k),real(f%blkt(:,ib),8)))
                   off = f%blkoff(ib)
                   do j = 1, nb
                      do i = 1, na
                         f%dm(off+i+(j-1)*na) = f%dm(off+i+(j-1)*na) + real(pa(i,b0+j) * phase,8)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine build_dm

  !> Position of distance r in the logarithmic grid of species s:
  !> spline interval i, offset dt within it, and dt/dr (tp).
  pure subroutine radial_index(s,r,i,dt,tp)
    type(fhiaimsspecies), intent(in) :: s
    real*8, intent(in) :: r
    integer, intent(out) :: i
    real*8, intent(out) :: dt, tp

    real*8 :: t

    t = 1d0 + log(r / s%rmin) / s%linc
    i = int(t)
    i = max(i,1)
    i = min(i,s%ngrid-1)
    dt = t - real(i,8)
    tp = 1d0 / (r * s%linc)

  end subroutine radial_index

  !> Evaluate radial function ifn of species s and its first two
  !> radial derivatives at distance r, given the spline position (i,
  !> dt, tp) from radial_index: u(r) is a cubic spline in the grid
  !> index, extrapolated from the last interval as in FHI-aims
  !> (val_spline). Below the first grid point, where the extrapolated
  !> cubic diverges, u = u(rmin) (r/rmin)^(l+1).
  pure subroutine radial_eval(s,ifn,r,i,dt,tp,nder,u,up,upp)
    type(fhiaimsspecies), intent(in) :: s
    integer, intent(in) :: ifn, i, nder
    real*8, intent(in) :: r, dt, tp
    real*8, intent(out) :: u, up, upp

    integer :: l
    real*8 :: ut, utt

    up = 0d0
    upp = 0d0
    if (r < s%rmin) then
       l = s%fn(ifn)%l
       u = s%fn(ifn)%spl(1,1) * (r / s%rmin)**(l+1)
       if (nder < 1) return
       up = (l+1) * u / r
       upp = l * (l+1) * u / r**2
       return
    end if

    associate (a => s%fn(ifn)%spl(:,i))
      u = a(1) + dt * (a(2) + dt * (a(3) + dt * a(4)))
      if (nder < 1) return
      ut = a(2) + dt * (2d0 * a(3) + dt * 3d0 * a(4))
      utt = 2d0 * a(3) + 6d0 * dt * a(4)
    end associate
    up = ut * tp
    upp = utt * tp * tp - ut * tp / r

  end subroutine radial_eval

  !> Accumulate the contraction of the (nr,nc) density matrix block d
  !> with the basis values of its row (phir, gr) and column (phic, gc)
  !> functions: wr += d phic, wur += d gc, and, if both, wc += d^T
  !> phir, wuc += d^T gr. The gradients are used only if nder >= 1.
  pure subroutine contract(d,nr,nc,both,nder,phir,phic,gr,gc,wr,wc,wur,wuc)
    integer, intent(in) :: nr, nc, nder
    real*8, intent(in) :: d(nr,nc)
    logical, intent(in) :: both
    real*8, intent(in) :: phir(*), phic(*), gr(3,*), gc(3,*)
    real*8, intent(inout) :: wr(*), wc(*), wur(3,*), wuc(3,*)

    integer :: i, j
    real*8 :: s, g(3)

    if (nder < 1) then
       do j = 1, nc
          wr(1:nr) = wr(1:nr) + d(:,j) * phic(j)
       end do
       if (both) then
          do j = 1, nc
             wc(j) = wc(j) + dot_product(d(:,j),phir(1:nr))
          end do
       end if
    else
       do j = 1, nc
          do i = 1, nr
             wr(i) = wr(i) + d(i,j) * phic(j)
             wur(:,i) = wur(:,i) + d(i,j) * gc(:,j)
          end do
       end do
       if (both) then
          do j = 1, nc
             s = 0d0
             g = 0d0
             do i = 1, nr
                s = s + d(i,j) * phir(i)
                g = g + d(i,j) * gr(:,i)
             end do
             wc(j) = wc(j) + s
             wuc(:,j) = wuc(:,j) + g
          end do
       end if
    end if

  end subroutine contract

end submodule proc
