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
  ! subroutine read_basis_and_occ(lu,f,ngrid,spcname,occ,nsv,errmsg)
  ! subroutine check_kgrid(nkxyz,kpt,errmsg)
  ! subroutine map_atoms(f,c,lat,xat,errmsg)
  ! subroutine build_blocks(f,c,xat)
  ! subroutine accumulate_dm(f,evk,occ,kpt,wk,nsv,ibas0)
  ! function lexneg(t)
  ! function block_key(f,dt) result(key)
  ! function find_block(f,a,b,dt,trans) result(ib)
  ! subroutine radial_index(s,r,i,dt,tp)
  ! subroutine radial_eval(s,ifn,r,i,dt,tp,nder,u,up,upp)

  ! first line of a dgrid file
  character(len=*), parameter :: dgrid_magic = "FHI-aims - DGrid Interface"
  ! hint appended to the structure mismatch messages
  character(len=*), parameter :: frame_hint = &
     "; load the structure from the dgrid file itself (CRYSTAL/MOLECULE file.dat)"
  ! minimum distance (bohr)
  real*8, parameter :: mindist = 1d-6
  ! occupations below this are not stored
  real*8, parameter :: occ_cutoff = 1d-14
  ! tolerance for comparing lattice vectors (bohr)
  real*8, parameter :: lat_eps = 1d-3

contains

  !> True if file is an FHI-aims dgrid_aims.dat file (output dgrid),
  !> from its first line.
  module function dgrid_is_file(file,ti) result(ok)
    use tools_io, only: fopen_read, fclose, getline_raw
    character*(*), intent(in) :: file
    type(thread_info), intent(in), optional :: ti
    logical :: ok

    integer :: lu
    character(len=:), allocatable :: line

    ok = .false.
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) return
    if (getline_raw(lu,line)) ok = (index(line,dgrid_magic) > 0)
    call fclose(lu)

  end function dgrid_is_file

  !> Read the header (first three lines) of a dgrid file open in
  !> unit lu, positioned at its start. errmsg is empty on success.
  module subroutine dgrid_read_header(lu,hdr,errmsg)
    use tools_io, only: getline_raw
    integer, intent(in) :: lu
    type(dgrid_header), intent(out) :: hdr
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: ios, idum
    character(len=:), allocatable :: line

    errmsg = "Error reading the header of the dgrid file"
    if (.not.getline_raw(lu,line)) return
    if (index(line,dgrid_magic) == 0) then
       errmsg = "Not a dgrid_aims.dat file (FHI-aims output dgrid)"
       return
    end if
    read (lu,*,iostat=ios) hdr%nspc, hdr%nat, idum, hdr%nstates, hdr%nbas, hdr%nfn, hdr%nper, hdr%nspin
    if (ios /= 0) return
    read (lu,*,iostat=ios) hdr%nkpt, hdr%nkxyz, idum
    if (ios /= 0) return
    if (hdr%nper /= 0 .and. hdr%nper /= 3) then
       errmsg = "Only molecules and 3D-periodic FHI-aims calculations are supported"
       return
    end if
    errmsg = ""

  end subroutine dgrid_read_header

  !> Read the grid sizes, the atoms (species index, Cartesian
  !> coordinates in bohr, species name) and the two energy lines that
  !> follow the header in a dgrid file. errmsg is empty on success.
  module subroutine dgrid_read_atoms(lu,hdr,ngrid,ispc,xat,spcname,errmsg)
    use tools_io, only: getline_raw
    integer, intent(in) :: lu
    type(dgrid_header), intent(in) :: hdr
    integer, allocatable, intent(inout) :: ngrid(:)
    integer, allocatable, intent(inout) :: ispc(:)
    real*8, allocatable, intent(inout) :: xat(:,:)
    character*20, allocatable, intent(inout) :: spcname(:)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, ios
    character*20 :: name
    character(len=:), allocatable :: line

    errmsg = "Error reading the atoms of the dgrid file"
    if (allocated(ngrid)) deallocate(ngrid)
    if (allocated(ispc)) deallocate(ispc)
    if (allocated(xat)) deallocate(xat)
    if (allocated(spcname)) deallocate(spcname)
    allocate(ngrid(hdr%nspc),ispc(hdr%nat),xat(3,hdr%nat),spcname(hdr%nspc))
    read (lu,*,iostat=ios) ngrid
    if (ios /= 0) return
    spcname = ""
    do i = 1, hdr%nat
       read (lu,*,iostat=ios) ispc(i), xat(:,i), name
       if (ios /= 0) return
       if (ispc(i) < 1 .or. ispc(i) > hdr%nspc) return
       spcname(ispc(i)) = name
    end do
    if (.not.getline_raw(lu,line)) return
    if (index(line,"Total energy") == 0) return
    if (.not.getline_raw(lu,line)) return
    if (index(line,"Chemical pot") == 0) return
    errmsg = ""

  end subroutine dgrid_read_atoms

  !> Read the lattice vectors (columns of lat, bohr), the fractional
  !> k-points and their weights from the end of a dgrid file: the
  !> last 6+nkpt non-blank lines, obtained with a stream read of the
  !> tail of the file (the eigenvector block before them is not
  !> parsed). The reciprocal vectors in between are checked against
  !> the lattice. errmsg is empty on success.
  module subroutine dgrid_read_tail(file,nkpt,lat,kpt,wk,errmsg)
    use tools_math, only: matinv
    use param, only: tpi
    character*(*), intent(in) :: file
    integer, intent(in) :: nkpt
    real*8, intent(out) :: lat(3,3)
    real*8, intent(out) :: kpt(3,nkpt)
    real*8, intent(out) :: wk(nkpt)
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: lu, ntail, i, j, n, ios, i0, i1, nbytes
    integer*8 :: fsize
    character(len=:), allocatable :: buf
    integer, allocatable :: lstart(:), lend(:)
    real*8 :: rinv(3,3), rec(3)
    logical :: ok

    ! read chunks from the end of the file until they hold ntail
    ! complete non-blank lines (the first line of a chunk may be cut)
    errmsg = "Error reading the lattice vectors and k-points at the end of the dgrid file (truncated?)"
    ntail = 6 + nkpt
    inquire(file=file,size=fsize)
    if (fsize <= 0) return
    open(newunit=lu,file=file,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if (ios /= 0) return
    nbytes = int(min(fsize,int(ntail+2,8) * 160_8))
    do
       if (allocated(buf)) deallocate(buf)
       allocate(character(len=nbytes) :: buf)
       read (lu,pos=fsize-int(nbytes,8)+1,iostat=ios) buf
       if (ios /= 0) then
          close(lu)
          return
       end if
       ! non-blank lines in the chunk, excluding the first (possibly cut) one
       n = 0
       i0 = 1
       do i = 1, nbytes
          if (buf(i:i) == new_line('a') .or. i == nbytes) then
             i1 = i
             if (buf(i:i) == new_line('a')) i1 = i - 1
             if (i0 > 1 .and. len_trim(buf(i0:i1)) > 0) n = n + 1
             i0 = i + 1
          end if
       end do
       ok = (n >= ntail) .or. (int(nbytes,8) == fsize)
       if (ok) exit
       nbytes = int(min(fsize,2_8*int(nbytes,8)))
    end do
    close(lu)

    ! positions of the last ntail non-blank lines
    allocate(lstart(ntail),lend(ntail))
    n = 0
    i0 = 1
    do i = 1, nbytes
       if (buf(i:i) == new_line('a') .or. i == nbytes) then
          i1 = i
          if (buf(i:i) == new_line('a')) i1 = i - 1
          if ((i0 > 1 .or. int(nbytes,8) == fsize) .and. len_trim(buf(i0:i1)) > 0) then
             n = n + 1
             lstart(mod(n-1,ntail)+1) = i0
             lend(mod(n-1,ntail)+1) = i1
          end if
          i0 = i + 1
       end if
    end do
    if (n < ntail) return

    ! the oldest of the kept lines is at position mod(n,ntail)+1
    do j = 1, 3
       read (buf(lstart(ln(j)):lend(ln(j))),*,iostat=ios) lat(:,j)
       if (ios /= 0) return
    end do
    rinv = lat
    call matinv(rinv,3)
    do j = 1, 3
       read (buf(lstart(ln(3+j)):lend(ln(3+j))),*,iostat=ios) rec
       if (ios /= 0) return
       if (any(abs(rec - tpi * rinv(j,:)) > 1d-5)) then
          errmsg = "The reciprocal vectors at the end of the dgrid file do not match the lattice (unknown layout?)"
          return
       end if
    end do
    do j = 1, nkpt
       read (buf(lstart(ln(6+j)):lend(ln(6+j))),*,iostat=ios) kpt(:,j), wk(j)
       if (ios /= 0) return
    end do
    if (abs(sum(wk) - 1d0) > 1d-6) then
       errmsg = "The k-point weights at the end of the dgrid file do not add up to one"
       return
    end if
    errmsg = ""

  contains
    ! ring position of the j-th kept line (1 = oldest)
    integer function ln(j)
      integer, intent(in) :: j
      ln = mod(n - ntail + j - 1,ntail) + 1
    end function ln
  end subroutine dgrid_read_tail

  !> Deallocate data
  module subroutine fhiaims_end(f)
    use iso_c_binding, only: c_null_ptr
    class(fhiaimswfn), intent(inout) :: f

    if (allocated(f%spc)) deallocate(f%spc)
    if (allocated(f%ispc)) deallocate(f%ispc)
    if (allocated(f%iaims)) deallocate(f%iaims)
    if (allocated(f%ioff)) deallocate(f%ioff)
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    if (allocated(f%pair0)) deallocate(f%pair0)
    if (allocated(f%pairn)) deallocate(f%pairn)
    if (allocated(f%blkt)) deallocate(f%blkt)
    if (allocated(f%blkkey)) deallocate(f%blkkey)
    if (allocated(f%blkoff)) deallocate(f%blkoff)
    if (allocated(f%dm)) deallocate(f%dm)
    f%nblk = 0
    f%cptr = c_null_ptr

  end subroutine fhiaims_end

  !> Read an FHI-aims wavefunction from the dgrid_aims.dat file
  !> (output dgrid keyword in control.in), map its atoms to the
  !> structure in cptr, and build the real-space density matrix. The
  !> lattice vectors and k-points, which the file keeps after the
  !> eigenvectors, are read first from the tail of the file; the
  !> eigenvectors are then streamed one k-point and spin at a time,
  !> so only one set of them is in memory.
  module subroutine fhiaims_read(f,cptr,file,errmsg,ti)
    use crystalmod, only: crystal
    use iso_c_binding, only: c_f_pointer
    use tools_io, only: fopen_read, fclose
    class(fhiaimswfn), intent(inout) :: f !< Output field
    type(c_ptr), intent(in) :: cptr !< Crystal structure
    character*(*), intent(in) :: file !< The dgrid_aims.dat file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, k, is, nsv, ios
    type(dgrid_header) :: hdr
    integer, allocatable :: ngrid(:), ibas0(:)
    character*20, allocatable :: spcname(:)
    real*8, allocatable :: occ(:,:,:), xat(:,:), xr(:), kpt(:,:), wk(:)
    complex*16, allocatable :: evk(:,:)
    real*8 :: lat(3,3)
    type(crystal), pointer :: c

    call c_f_pointer(cptr,c)
    f%cptr = cptr

    ! header, atoms, basis and occupations
    lu = fopen_read(file,errstop=.false.,ti=ti)
    if (lu < 0) then
       errmsg = "Error opening file: " // trim(file)
       return
    end if
    call dgrid_read_header(lu,hdr,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    call dgrid_read_atoms(lu,hdr,ngrid,f%ispc,xat,spcname,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    f%nspc = hdr%nspc
    f%nat = hdr%nat
    f%nstates = hdr%nstates
    f%nbas = hdr%nbas
    f%nspin = hdr%nspin
    f%nkpt = hdr%nkpt
    f%ismol = (hdr%nper == 0)
    call read_basis_and_occ(lu,f,hdr,ngrid,spcname,occ,nsv,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! lattice vectors and k-points, from the end of the file
    allocate(kpt(3,f%nkpt),wk(f%nkpt))
    if (f%ismol) then
       lat = 0d0
       kpt = 0d0
       wk = 1d0
    else
       call dgrid_read_tail(file,f%nkpt,lat,kpt,wk,errmsg)
       if (len_trim(errmsg) > 0) goto 999
       call check_kgrid(hdr%nkxyz,kpt,errmsg)
       if (len_trim(errmsg) > 0) goto 999
    end if

    ! map the atoms and set up the density matrix blocks
    call map_atoms(f,c,lat,xat,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    call build_blocks(f,c,xat)
    allocate(ibas0(f%nat))
    k = 0
    do i = 1, f%nat
       ibas0(i) = k
       k = k + f%spc(f%ispc(i))%nbas
    end do

    ! eigenvectors: accumulate the density matrix per k-point and spin
    errmsg = "Error reading the eigenvectors of the dgrid file: " // trim(file)
    allocate(evk(f%nbas,nsv))
    if (f%isreal) then
       allocate(xr(f%nbas))
    else
       allocate(xr(2*f%nbas))
    end if
    do k = 1, f%nkpt
       do is = 1, f%nspin
          ! only the states up to the last occupied one are kept
          do i = 1, f%nstates
             read (lu,*,iostat=ios) xr
             if (ios /= 0) goto 999
             if (i > nsv) cycle
             if (f%isreal) then
                evk(:,i) = xr
             else
                evk(:,i) = cmplx(xr(1::2),xr(2::2),8)
             end if
          end do
          call accumulate_dm(f,evk,occ(:,is,k),kpt(:,k),wk(k),nsv,ibas0)
       end do
    end do
    call fclose(lu)

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
    integer :: nloc, nd, ie, je, ib, i, ir, ic, nr, nc, dt(3)
    logical :: trans
    real*8 :: r, tp(2), u, up, upp, rr, rp, rpp, dspl, tspl, iphas
    integer, allocatable :: eid(:), lvec(:,:), iat(:), tvec(:,:), loc0(:)
    real*8, allocatable :: dist(:), xenv(:,:), rc(:), p(:,:), phipp(:,:), w(:,:)
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
    allocate(iat(nenv),tvec(3,nenv),loc0(nenv+1),xenv(3,nenv),rc(nenv))
    nenvl = 0
    nloc = 0
    do ion = 1, nenv
       ia = f%iaims(eid(ion))
       if (dist(ion) > f%spc(f%ispc(ia))%rcut) cycle
       nenvl = nenvl + 1
       iat(nenvl) = ia
       rc(nenvl) = f%spc(f%ispc(ia))%rcut
       tvec(:,nenvl) = lvec(:,ion) - f%ioff(:,eid(ion))
       xenv(:,nenvl) = xpos - (c%atcel(eid(ion))%r + c%x2c(real(lvec(:,ion),8)))
       loc0(nenvl) = nloc
       nloc = nloc + f%spc(f%ispc(ia))%nbas
    end do
    loc0(nenvl+1) = nloc
    if (nenvl == 0) return

    ! values and gradients of the local basis functions: p(:,1) =
    ! phi, p(:,2:4) = grad phi; second derivatives in phipp
    nd = 1
    if (nder >= 1) nd = 4
    allocate(p(nloc,nd),ylm((f%lmax+3)**2))
    p = 0d0
    if (nder >= 2) then
       allocate(phipp(6,nloc))
       phipp = 0d0
    end if
    do ie = 1, nenvl
       is = f%ispc(iat(ie))
       call tosphere(xenv(:,ie),r,tp)
       r = max(r,mindist)
       call genylm(f%spc(is)%lmax+2,tp,ylm)
       call radial_index(f%spc(is),r,ispl,dspl,tspl)

       iloc = loc0(ie)
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
          p(ip,1) = rr * real(ylm(iy),8)
          if (nder >= 1) p(ip,2:4) = real(xgrad,8)
          if (nder >= 2) phipp(:,ip) = real(xhess,8)
          iphas = 1d0
          do m = 1, l
             iphas = -iphas
             ip = iloc + l + 1 + m
             im = iloc + l + 1 - m
             call ylmderiv(ylm,r,l,m,rr,rp,rpp,nder,xgrad,xhess)
             p(ip,1) = sqrt2 * rr * real(ylm(iy+m),8)
             p(im,1) = sqrt2 * iphas * rr * aimag(ylm(iy+m))
             if (nder < 1) cycle
             p(ip,2:4) = sqrt2 * real(xgrad,8)
             p(im,2:4) = sqrt2 * iphas * aimag(xgrad)
             if (nder < 2) cycle
             phipp(:,ip) = sqrt2 * real(xhess,8)
             phipp(:,im) = sqrt2 * iphas * aimag(xhess)
          end do
          iloc = iloc + 2*l+1
       end do
    end do

    ! contract with the density matrix, w = D p, one stored block at
    ! a time; the transposed block gives the (q,p) contribution
    allocate(w(nloc,nd))
    w = 0d0
    do ie = 1, nenvl
       do je = ie, nenvl
          if (norm2(xenv(:,ie)-xenv(:,je)) > rc(ie) + rc(je)) cycle
          dt = tvec(:,ie) - tvec(:,je)
          ib = find_block(f,iat(ie),iat(je),dt,trans)
          if (ib == 0) cycle
          if (trans) then
             ir = je
             ic = ie
          else
             ir = ie
             ic = je
          end if
          nr = loc0(ir+1) - loc0(ir)
          nc = loc0(ic+1) - loc0(ic)
          call dgemm('N','N',nr,nd,nc,1d0,f%dm(f%blkoff(ib)+1),nr,p(loc0(ic)+1,1),nloc,1d0,w(loc0(ir)+1,1),nloc)
          if (ie /= je) &
             call dgemm('T','N',nc,nd,nr,1d0,f%dm(f%blkoff(ib)+1),nr,p(loc0(ir)+1,1),nloc,1d0,w(loc0(ic)+1,1),nloc)
       end do
    end do

    ! accumulate properties
    rho = dot_product(p(:,1),w(:,1))
    if (nder >= 1) then
       do i = 1, 3
          grad(i) = 2d0 * dot_product(p(:,1+i),w(:,1))
          gkin = gkin + 0.5d0 * dot_product(p(:,1+i),w(:,1+i))
       end do
    end if
    if (nder >= 2) then
       do i = 1, nloc
          h(1,1) = h(1,1) + 2d0 * (phipp(1,i) * w(i,1) + p(i,2) * w(i,2))
          h(1,2) = h(1,2) + 2d0 * (phipp(2,i) * w(i,1) + p(i,2) * w(i,3))
          h(1,3) = h(1,3) + 2d0 * (phipp(3,i) * w(i,1) + p(i,2) * w(i,4))
          h(2,2) = h(2,2) + 2d0 * (phipp(4,i) * w(i,1) + p(i,3) * w(i,3))
          h(2,3) = h(2,3) + 2d0 * (phipp(5,i) * w(i,1) + p(i,3) * w(i,4))
          h(3,3) = h(3,3) + 2d0 * (phipp(6,i) * w(i,1) + p(i,4) * w(i,4))
       end do
       h(2,1) = h(1,2)
       h(3,1) = h(1,3)
       h(3,2) = h(2,3)
    end if

  end subroutine rho2

  !xx! private procedures

  !> Read the part of a dgrid file between the atoms and the
  !> eigenvectors: radial functions, radial grids, spline
  !> coefficients (into f%spc), eigenvalues and occupations (occ,
  !> not k-weighted), and the real/complex flag (f%isreal). nsv is
  !> the last state with a non-zero occupation. errmsg is empty on
  !> success.
  subroutine read_basis_and_occ(lu,f,hdr,ngrid,spcname,occ,nsv,errmsg)
    integer, intent(in) :: lu
    type(fhiaimswfn), intent(inout) :: f
    type(dgrid_header), intent(in) :: hdr
    integer, intent(in) :: ngrid(:)
    character*20, intent(in) :: spcname(:)
    real*8, allocatable, intent(inout) :: occ(:,:,:)
    integer, intent(out) :: nsv
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: i, j, k, is, ios, iflag, idum
    integer, allocatable :: fsp(:), fl(:), nfnsp(:), ifn(:)
    real*8, allocatable :: frc(:)
    real*8 :: rdum, rinc

    errmsg = "Error reading the basis set and occupations of the dgrid file"
    if (allocated(f%spc)) deallocate(f%spc)
    allocate(f%spc(f%nspc))
    do is = 1, f%nspc
       f%spc(is)%name = spcname(is)
       f%spc(is)%ngrid = ngrid(is)
    end do

    ! radial functions
    allocate(fsp(hdr%nfn),fl(hdr%nfn),frc(hdr%nfn),nfnsp(f%nspc),ifn(f%nspc))
    do j = 1, hdr%nfn
       read (lu,*,iostat=ios) fsp(j), fl(j), idum, frc(j)
       if (ios /= 0) return
       if (fsp(j) < 1 .or. fsp(j) > f%nspc) return
    end do
    read (lu,*,iostat=ios) nfnsp
    if (ios /= 0) return
    do is = 1, f%nspc
       f%spc(is)%nfn = count(fsp == is)
       if (f%spc(is)%nfn /= nfnsp(is)) then
          errmsg = "Inconsistent number of radial functions in the dgrid file"
          return
       end if
       allocate(f%spc(is)%fn(f%spc(is)%nfn))
    end do

    ! radial grids and spline coefficients
    do is = 1, f%nspc
       read (lu,*,iostat=ios) f%spc(is)%rmin, rinc
       if (ios /= 0) return
       f%spc(is)%linc = log(rinc)
    end do
    ifn = 0
    do j = 1, hdr%nfn
       is = fsp(j)
       ifn(is) = ifn(is) + 1
       f%spc(is)%fn(ifn(is))%l = fl(j)
       f%spc(is)%fn(ifn(is))%rcut = frc(j)
       allocate(f%spc(is)%fn(ifn(is))%spl(4,f%spc(is)%ngrid))
       read (lu,*,iostat=ios) f%spc(is)%fn(ifn(is))%spl
       if (ios /= 0) return
    end do
    f%lmax = 0
    k = 0
    do is = 1, f%nspc
       f%spc(is)%nbas = sum(2 * f%spc(is)%fn(:)%l + 1)
       f%spc(is)%rcut = maxval(f%spc(is)%fn(:)%rcut)
       f%spc(is)%lmax = maxval(f%spc(is)%fn(:)%l)
       f%lmax = max(f%lmax,f%spc(is)%lmax)
    end do
    do i = 1, f%nat
       k = k + f%spc(f%ispc(i))%nbas
    end do
    if (k /= f%nbas) then
       errmsg = "Inconsistent number of basis functions in the dgrid file"
       return
    end if

    ! eigenvalues and occupations
    if (allocated(occ)) deallocate(occ)
    allocate(occ(f%nstates,f%nspin,f%nkpt))
    do k = 1, f%nkpt
       do is = 1, f%nspin
          do i = 1, f%nstates
             read (lu,*,iostat=ios) rdum, occ(i,is,k)
             if (ios /= 0) return
          end do
       end do
    end do
    nsv = 0
    do i = 1, f%nstates
       if (any(abs(occ(i,:,:)) > occ_cutoff)) nsv = i
    end do

    ! real (-5) or complex (5) eigenvectors
    read (lu,*,iostat=ios) iflag
    if (ios /= 0) return
    f%isreal = (iflag < 0)
    errmsg = ""

  end subroutine read_basis_and_occ

  !> Check that the k-point list (fractional coordinates) samples the
  !> full nkxyz grid: expanding every k to {k,-k} modulo G must give
  !> at least product(nkxyz) distinct points (exactly that many for
  !> the usual time-reversal reduction, more for an unreduced shifted
  !> grid). A space-group-reduced list (symmetry_reduced_k_grid_spg)
  !> gives fewer and would yield the unsymmetrized density. Skipped
  !> if the header does not give the grid (external k-point lists).
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
    if (n < product(nkxyz)) then
       errmsg = "The " // string(nk) // " k-points do not cover the " // string(nkxyz(1)) // "x" //&
          string(nkxyz(2)) // "x" // string(nkxyz(3)) // " grid under time reversal (" //&
          string(n) // " points); space-group-reduced k-grids are not supported, rerun FHI-aims &
          &without symmetry_reduced_k_grid_spg"
    end if

  end subroutine check_kgrid

  !> Map the FHI-aims atoms (Cartesian positions xat, bohr) to the
  !> cell atoms of crystal c. The two structures must have the same
  !> Cartesian frame: the lattice vectors (lat) must coincide for
  !> crystals, and molecules may only be translated (critic2
  !> recenters them by c%molx0).
  subroutine map_atoms(f,c,lat,xat,errmsg)
    use crystalmod, only: crystal
    use tools_io, only: zatguess
    use param, only: icrd_cart
    type(fhiaimswfn), intent(inout) :: f
    type(crystal), intent(inout) :: c
    real*8, intent(in) :: lat(3,3)
    real*8, intent(in) :: xat(:,:)
    character(len=:), allocatable, intent(inout) :: errmsg

    integer :: i, ia, iz, lv(3)

    errmsg = ""
    if (c%ncel /= f%nat) then
       errmsg = "The number of atoms in the FHI-aims file does not match the structure &
          &(ghost/empty atoms are not supported)" // frame_hint
       return
    end if
    if (c%ismolecule .neqv. f%ismol) then
       errmsg = "The periodicity of the FHI-aims calculation does not match the structure" // frame_hint
       return
    end if
    if (.not.f%ismol) then
       if (any(abs(c%m_x2c - lat) > lat_eps)) then
          errmsg = "The lattice vectors of the FHI-aims calculation do not match those of the &
             &structure" // frame_hint
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
       i = c%identify_atom(xat(:,ia) - c%molx0,icrd_cart,lvec=lv)
       if (i == 0) then
          errmsg = "Could not match the atoms of the FHI-aims calculation to the structure" // frame_hint
          return
       elseif (f%iaims(i) /= 0) then
          errmsg = "Two FHI-aims atoms match the same atom of the structure" // frame_hint
          return
       end if
       iz = zatguess(f%spc(f%ispc(ia))%name)
       if (iz > 0 .and. iz /= c%spc(c%atcel(i)%is)%z) then
          errmsg = "The atomic species of the FHI-aims calculation do not match the structure" // frame_hint
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
  !> on b+T_b overlap (xat: FHI-aims Cartesian positions). Only a <= b
  !> is stored (for a == b, dT lexicographically non-negative); the
  !> block (b,a,-dT) is the transpose. Allocates dm (zeroed).
  subroutine build_blocks(f,c,xat)
    use crystalmod, only: crystal
    use types, only: realloc
    type(fhiaimswfn), intent(inout) :: f
    type(crystal), intent(in) :: c
    real*8, intent(in) :: xat(:,:)

    integer :: a, b, i, j, k, n, nbox, mblk
    integer*8 :: ntot
    real*8 :: rmax, dmax, v(3)
    integer, allocatable :: tbox(:,:)
    real*8, allocatable :: lbox(:,:)

    ! lattice vectors that can connect overlapping atoms: |dT_i| <=
    ! |row i of m_c2x| * (2 rcut + max interatomic distance)
    f%tmax = 0
    if (.not.f%ismol) then
       rmax = 2d0 * maxval(f%spc(:)%rcut)
       dmax = 0d0
       do a = 1, f%nat
          do b = a+1, f%nat
             dmax = max(dmax,norm2(xat(:,a)-xat(:,b)))
          end do
       end do
       do i = 1, 3
          f%tmax(i) = ceiling(c%ar(i) * (rmax + dmax))
       end do
    end if
    nbox = product(2*f%tmax+1)
    allocate(tbox(3,nbox),lbox(3,nbox))
    n = 0
    do k = -f%tmax(3), f%tmax(3)
       do j = -f%tmax(2), f%tmax(2)
          do i = -f%tmax(1), f%tmax(1)
             n = n + 1
             tbox(:,n) = (/i,j,k/)
             lbox(:,n) = c%x2c(real(tbox(:,n),8))
          end do
       end do
    end do

    ! find the blocks; the box is traversed in increasing key order
    if (allocated(f%pair0)) deallocate(f%pair0)
    if (allocated(f%pairn)) deallocate(f%pairn)
    if (allocated(f%blkt)) deallocate(f%blkt)
    allocate(f%pair0(f%nat,f%nat),f%pairn(f%nat,f%nat))
    f%pair0 = 0
    f%pairn = 0
    mblk = 4 * f%nat * nbox
    allocate(f%blkt(3,mblk))
    f%nblk = 0
    do a = 1, f%nat
       do b = a, f%nat
          f%pair0(a,b) = f%nblk + 1
          do n = 1, nbox
             if (a == b .and. lexneg(tbox(:,n))) cycle
             v = xat(:,a) - xat(:,b) + lbox(:,n)
             if (norm2(v) > f%spc(f%ispc(a))%rcut + f%spc(f%ispc(b))%rcut) cycle
             f%nblk = f%nblk + 1
             if (f%nblk > mblk) then
                mblk = 2 * mblk
                call realloc(f%blkt,3,mblk)
             end if
             f%pairn(a,b) = f%pairn(a,b) + 1
             f%blkt(:,f%nblk) = tbox(:,n)
          end do
       end do
    end do
    call realloc(f%blkt,3,f%nblk)

    ! keys and offsets
    if (allocated(f%blkkey)) deallocate(f%blkkey)
    if (allocated(f%blkoff)) deallocate(f%blkoff)
    allocate(f%blkkey(f%nblk),f%blkoff(f%nblk))
    ntot = 0
    do a = 1, f%nat
       do b = a, f%nat
          do i = f%pair0(a,b), f%pair0(a,b) + f%pairn(a,b) - 1
             f%blkkey(i) = block_key(f,f%blkt(:,i))
             f%blkoff(i) = ntot
             ntot = ntot + f%spc(f%ispc(a))%nbas * f%spc(f%ispc(b))%nbas
          end do
       end do
    end do
    if (allocated(f%dm)) deallocate(f%dm)
    allocate(f%dm(ntot))
    f%dm = 0d0

  end subroutine build_blocks

  !> Accumulate the contribution of one spin channel and k-point
  !> (fractional kpt, weight wk) to the real-space density matrix
  !> blocks, from the eigenvectors evk (nbas,nsv) and the occupations
  !> (occ, not k-weighted):
  !>   D_ij(dT) += Re wk sum_n occ_n c_in c_jn^* exp(2 pi i k.dT)
  !> The real part is exact: the imaginary parts cancel in the
  !> symmetric sum over basis function pairs. Only the columns of the
  !> atoms b >= a that share a block with a are computed. ibas0 is
  !> the first basis function of each atom minus one.
  subroutine accumulate_dm(f,evk,occ,kpt,wk,nsv,ibas0)
    use param, only: tpi, img
    type(fhiaimswfn), intent(inout) :: f
    complex*16, intent(in) :: evk(:,:)
    real*8, intent(in) :: occ(:)
    real*8, intent(in) :: kpt(3)
    real*8, intent(in) :: wk
    integer, intent(in) :: nsv
    integer, intent(in) :: ibas0(:)

    integer :: a, b, na, nb, a0, b0, n, ib, ncol
    integer*8 :: off
    integer, allocatable :: col0(:)
    complex*16, allocatable :: fw(:,:), en(:,:), pa(:,:)
    complex*16 :: phase

    complex*16, parameter :: cone = (1d0,0d0), czero = (0d0,0d0)

    if (nsv == 0) return

    ! occupation- and weight-scaled eigenvectors
    allocate(fw(f%nbas,nsv),col0(f%nat))
    do n = 1, nsv
       fw(:,n) = evk(:,n) * (occ(n) * wk)
    end do

    do a = 1, f%nat
       na = f%spc(f%ispc(a))%nbas
       a0 = ibas0(a)

       ! gather the eigenvector rows of the neighbours of a
       ncol = 0
       do b = a, f%nat
          if (f%pairn(a,b) == 0) cycle
          col0(b) = ncol
          ncol = ncol + f%spc(f%ispc(b))%nbas
       end do
       if (allocated(en)) deallocate(en)
       if (allocated(pa)) deallocate(pa)
       allocate(en(ncol,nsv),pa(na,ncol))
       do b = a, f%nat
          if (f%pairn(a,b) == 0) cycle
          nb = f%spc(f%ispc(b))%nbas
          b0 = ibas0(b)
          en(col0(b)+1:col0(b)+nb,:) = evk(b0+1:b0+nb,1:nsv)
       end do

       ! the a-rows of the k-space density matrix, neighbour columns
       call zgemm('N','C',na,ncol,nsv,cone,fw(a0+1,1),f%nbas,en,ncol,czero,pa,na)

       do b = a, f%nat
          if (f%pairn(a,b) == 0) cycle
          nb = f%spc(f%ispc(b))%nbas
          n = na * nb
          do ib = f%pair0(a,b), f%pair0(a,b) + f%pairn(a,b) - 1
             phase = exp(img * tpi * dot_product(kpt,real(f%blkt(:,ib),8)))
             off = f%blkoff(ib)
             f%dm(off+1:off+n) = f%dm(off+1:off+n) + reshape(real(pa(:,col0(b)+1:col0(b)+nb) * phase,8),(/n/))
          end do
       end do
    end do

  end subroutine accumulate_dm

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

  !> Integer key of a block lattice vector: the position of dT in the
  !> -tmax..tmax box, in the order build_blocks traverses it.
  pure function block_key(f,dt) result(key)
    type(fhiaimswfn), intent(in) :: f
    integer, intent(in) :: dt(3)
    integer :: key

    key = (dt(1) + f%tmax(1)) + (2*f%tmax(1)+1) * ((dt(2) + f%tmax(2)) + (2*f%tmax(2)+1) * (dt(3) + f%tmax(3)))

  end function block_key

  !> Find the block for atoms a, b and lattice vector dT = T_a - T_b
  !> by bisection on the keys of the pair. Returns 0 if the block
  !> does not exist. trans is true if the stored block is (b,a,-dT),
  !> the transpose of the requested one.
  function find_block(f,a,b,dt,trans) result(ib)
    type(fhiaimswfn), intent(in) :: f
    integer, intent(in) :: a, b, dt(3)
    logical, intent(out) :: trans
    integer :: ib

    integer :: aa, bb, tt(3), key, lo, hi, mid

    ib = 0
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
    if (any(abs(tt) > f%tmax)) return
    key = block_key(f,tt)
    lo = f%pair0(aa,bb)
    hi = lo + f%pairn(aa,bb) - 1
    do while (lo <= hi)
       mid = (lo + hi) / 2
       if (f%blkkey(mid) == key) then
          ib = mid
          return
       elseif (f%blkkey(mid) < key) then
          lo = mid + 1
       else
          hi = mid - 1
       end if
    end do

  end function find_block

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

end submodule proc
