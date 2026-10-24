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
  ! subroutine sb_open(sb,file,errmsg)
  ! subroutine sb_fill(sb,ios)
  ! subroutine sb_values(sb,nval,ios,x)
  ! subroutine accumulate_dm(f,evk,evr,occ,kpt,wk,nsv,is)
  ! function max_nloc(f,c) result(nloc)
  ! subroutine block_norms(f)
  ! subroutine blockmul(op,m,n,k,a,lda,b,ldb,beta,c,ldc)
  ! subroutine accum_props(x,y,nder,rho,grad,h,gkin)
  ! subroutine spin_fill(nsw,rho,grad,h,gkin)
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
  ! a pair of atoms is skipped at a point if the bound on its
  ! contribution, |D|max * sum|p_a| * sum|p_b|, is below this
  real*8, parameter :: dm_pairthr = 1d-14

  ! blank characters in the eigenvector block (space, tab, LF, CR)
  character(len=4), parameter :: blanks = " " // char(9) // char(10) // char(13)

  !> Chunked reader of the eigenvector block of a dgrid file as a
  !> byte stream. The block is 95% of the file, and a list-directed
  !> read of it costs several times more than the conversion of the
  !> numbers (and most of the states are unoccupied and skipped).
  integer, parameter :: sb_size = 1048576
  type streambuf
     integer :: lu = -1 !< unit (stream access), -1 if closed
     integer*8 :: fsize !< size of the file (bytes)
     integer*8 :: fpos !< position of the next byte to read from the file
     integer :: n = 0 !< number of bytes in buf
     integer :: ip = 1 !< next unread byte in buf
     logical :: eof = .false. !< the whole file has been read into buf
     character(len=:), allocatable :: buf !< the buffer (sb_size, plus a blank sentinel after the data)
  end type streambuf

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
    if (allocated(f%ibas0)) deallocate(f%ibas0)
    if (allocated(f%cs)) deallocate(f%cs)
    if (allocated(f%iaims)) deallocate(f%iaims)
    if (allocated(f%ioff)) deallocate(f%ioff)
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    if (allocated(f%pair0)) deallocate(f%pair0)
    if (allocated(f%pairn)) deallocate(f%pairn)
    if (allocated(f%blkt)) deallocate(f%blkt)
    if (allocated(f%blkkey)) deallocate(f%blkkey)
    if (allocated(f%blkoff)) deallocate(f%blkoff)
    if (allocated(f%blknorm)) deallocate(f%blknorm)
    if (allocated(f%dm)) deallocate(f%dm)
    f%nocc = 0
    f%noccs = 0
    f%nblk = 0
    f%cptr = c_null_ptr

  end subroutine fhiaims_end

  !> Read an FHI-aims wavefunction from the dgrid_aims.dat file
  !> (output dgrid keyword in control.in), map its atoms to the
  !> structure in cptr, and build the real-space density matrix or,
  !> when that is cheaper, keep the occupied orbitals (real
  !> eigenvectors at Gamma only, non-negative occupations, and fewer
  !> occupied orbitals than local basis functions at the atoms: the
  !> orbital route costs nocc*nloc per point, the density matrix
  !> ~nloc^2/2). The lattice vectors and k-points, which the file
  !> keeps after the eigenvectors, are read first from the tail of
  !> the file; the eigenvectors are then streamed one k-point and
  !> spin at a time, so only one set of them is in memory.
  module subroutine fhiaims_read(f,cptr,file,errmsg,ti)
    use crystalmod, only: crystal
    use iso_c_binding, only: c_f_pointer
    use tools_io, only: fopen_read, fclose
    class(fhiaimswfn), intent(inout) :: f !< Output field
    type(c_ptr), intent(in) :: cptr !< Crystal structure
    character*(*), intent(in) :: file !< The dgrid_aims.dat file
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, k, is, nsv, ios, icol, nocc, nval
    type(dgrid_header) :: hdr
    integer, allocatable :: ngrid(:), iperm(:)
    character*20, allocatable :: spcname(:)
    real*8, allocatable :: occ(:,:,:), xat(:,:), xr(:), kpt(:,:), wk(:), evr(:,:)
    complex*16, allocatable :: evk(:,:)
    real*8 :: lat(3,3)
    type(crystal), pointer :: c
    type(streambuf) :: sb

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
    call fclose(lu)
    lu = -1

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

    ! map the atoms
    call map_atoms(f,c,lat,xat,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    if (allocated(f%ibas0)) deallocate(f%ibas0)
    allocate(f%ibas0(f%nat),iperm(f%nbas))
    k = 0
    do i = 1, f%nat
       f%ibas0(i) = k
       ! file-order index of each local (sorted) basis function
       iperm(k+1:k+f%spc(f%ispc(i))%nbas) = k + f%spc(f%ispc(i))%ifile
       k = k + f%spc(f%ispc(i))%nbas
    end do

    ! the eigenvectors, as a byte stream
    call sb_open(sb,file,errmsg)
    if (len_trim(errmsg) > 0) goto 999
    errmsg = "Error reading the eigenvectors of the dgrid file: " // trim(file)
    nocc = count(occ(:,:,1) > occ_cutoff)
    if (f%isreal .and. f%nkpt == 1 .and. all(abs(kpt(:,1)) < 1d-10) .and. all(occ >= -occ_cutoff) .and.&
       nocc <= max_nloc(f,c)) then
       ! occupied orbitals, scaled by sqrt(occ), rows in local order
       f%nocc = nocc
       if (allocated(f%cs)) deallocate(f%cs)
       allocate(f%cs(f%nbas,f%nocc),xr(f%nbas))
       icol = 0
       f%noccs = 0
       do is = 1, f%nspin
          do i = 1, f%nstates
             if (occ(i,is,1) <= occ_cutoff) then
                call sb_values(sb,f%nbas,ios)
                if (ios /= 0) goto 999
                cycle
             end if
             call sb_values(sb,f%nbas,ios,xr)
             if (ios /= 0) goto 999
             icol = icol + 1
             f%noccs(is) = f%noccs(is) + 1
             f%cs(:,icol) = xr(iperm) * sqrt(occ(i,is,1) * wk(1))
          end do
       end do
    else
       ! density matrix, accumulated per k-point and spin
       call build_blocks(f,c,xat)
       if (f%isreal) then
          allocate(evr(f%nbas,nsv),evk(0,0))
          nval = f%nbas
       else
          allocate(evk(f%nbas,nsv),evr(0,0),xr(2*f%nbas))
          nval = 2 * f%nbas
       end if
       do k = 1, f%nkpt
          do is = 1, f%nspin
             ! only the states up to the last occupied one are kept
             do i = 1, nsv
                if (f%isreal) then
                   call sb_values(sb,nval,ios,evr(:,i))
                else
                   call sb_values(sb,nval,ios,xr)
                   evk(:,i) = cmplx(xr(1::2),xr(2::2),8)
                end if
                if (ios /= 0) goto 999
             end do
             call sb_values(sb,(f%nstates-nsv)*nval,ios)
             if (ios /= 0) goto 999
             call accumulate_dm(f,evk,evr,occ(:,is,k),kpt(:,k),wk(k),nsv,is)
          end do
       end do
       call block_norms(f)
    end if
    close(sb%lu)

    errmsg = ""
    return
999 continue
    if (lu >= 0) call fclose(lu)
    if (sb%lu >= 0) close(sb%lu)

  end subroutine fhiaims_read

  !> Calculate the density and derivatives of an FHI-aims field (f)
  !> up to the nder degree (max = 2). xpos is in Cartesian
  !> coordinates. In output, the density (rho), the gradient (grad),
  !> the Hessian (h), and the G(r) kinetic energy density (gkin); the
  !> last index is 1 for the total, 2 for spin up and 3 for spin down.
  !> The spin components are filled only if dospin and the
  !> calculation is spin-polarized (at twice the cost of the total
  !> alone on the density matrix route). This routine is thread-safe.
  module subroutine rho2(f,xpos,nder,dospin,rho,grad,h,gkin)
    use crystalmod, only: crystal
    use iso_c_binding, only: c_f_pointer
    use tools_math, only: tosphere, genylm, ylmderiv_real
    use tools, only: qcksort
    use param, only: icrd_cart
    class(fhiaimswfn), intent(inout) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    integer, intent(in) :: nder  !< Number of derivatives
    logical, intent(in) :: dospin !< Compute the spin densities of a spin-polarized calculation
    real*8, intent(out) :: rho(3) !< Density (total, spin up, spin down)
    real*8, intent(out) :: grad(3,3) !< Gradient (:,1:3 = total, spin up, spin down)
    real*8, intent(out) :: h(3,3,3) !< Hessian (:,:,1:3 = total, spin up, spin down)
    real*8, intent(out) :: gkin(3) !< G(r), kinetic energy density (total, spin up, spin down)

    integer :: nenv, ion, nenvl, ia, is, ifn, l, m, ip, im, iloc, ispl
    integer :: nloc, nd, ndw, ie, je, ib, k, k1, k2, ir, ic, nr, mr, mc, n, n2, np, nbmax, maxrun, dt(3)
    integer :: nsw, js
    logical :: trans
    real*8 :: r, tp(2), u, up, upp, rr, rp, rpp, dspl, tspl, iphas
    integer, allocatable :: eid(:), lvec(:,:), iat(:), tvec(:,:), loc0(:), nact(:)
    integer, allocatable :: pib(:), pir(:), pic(:), iord(:)
    logical, allocatable :: pboth(:)
    real*8, allocatable :: dist(:), xenv(:,:), rc(:), psum(:), p(:,:), w(:,:,:), g(:,:), wk(:,:), psi(:,:)
    complex*16, allocatable :: ylm(:)
    real*8 :: s(2), sg(3,2), sh(6,2)
    type(crystal), pointer :: c

    rho = 0d0
    grad = 0d0
    h = 0d0
    gkin = 0d0
    nsw = 1
    if (dospin .and. f%nspin == 2) nsw = 2

    ! atoms around the input point; keep the ones within their own cutoff
    call c_f_pointer(f%cptr,c)
    call c%list_near_atoms(xpos,icrd_cart,.false.,nenv,eid,dist,lvec,up2dsp=f%spcutoff)
    if (nenv == 0) return
    allocate(iat(nenv),tvec(3,nenv),loc0(nenv+1),xenv(3,nenv),rc(nenv),nact(nenv))
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

    ! values and derivatives of the local basis functions: p(:,1) =
    ! phi, p(:,2:4) = grad phi, p(:,5:10) = second derivatives
    ! (xx,xy,xz,yy,yz,zz)
    nd = 1
    if (nder >= 1) nd = 4
    if (nder >= 2) nd = 10
    allocate(p(nloc,nd),ylm((f%lmax+nder+1)**2))
    p = 0d0
    do ie = 1, nenvl
       is = f%ispc(iat(ie))
       call tosphere(xenv(:,ie),r,tp)
       r = max(r,mindist)
       call genylm(f%spc(is)%lmax+nder,tp,ylm)
       call radial_index(f%spc(is),r,ispl,dspl,tspl)

       iloc = loc0(ie)
       do ifn = 1, f%spc(is)%nfn
          l = f%spc(is)%fn(ifn)%l
          ! by decreasing outer radius: this and the rest are zero here
          if (r > f%spc(is)%fn(ifn)%rcut) exit

          ! radial part, R = u/r
          call radial_eval(f%spc(is),ifn,r,ispl,dspl,tspl,nder,u,up,upp)
          rr = u / r
          rp = up / r - u / r**2
          rpp = upp / r - 2d0 * up / r**2 + 2d0 * u / r**3

          ! real spherical harmonics in the FHI-aims convention
          ! (SHEval): S_l0 = Y_l0, S_lm = sqrt2 Re(Y_lm) and S_l-m =
          ! (-1)^m sqrt2 Im(Y_lm) for m > 0, with Condon-Shortley Y_lm.
          ! Local basis order is m = -l..l.
          do m = 0, l
             iphas = real((-1)**m,8)
             ip = iloc + l + 1 + m
             call ylmderiv_real(ylm,r,l,m,rr,rp,rpp,nder,s,sg,sh)
             p(ip,1) = s(1)
             if (nder >= 1) p(ip,2:4) = sg(:,1)
             if (nder >= 2) p(ip,5:10) = sh(:,1)
             if (m == 0) cycle
             im = iloc + l + 1 - m
             p(im,1) = iphas * s(2)
             if (nder >= 1) p(im,2:4) = iphas * sg(:,2)
             if (nder >= 2) p(im,5:10) = iphas * sh(:,2)
          end do
          iloc = iloc + 2*l+1
       end do
       ! the functions of this atom that are non-zero at the point
       nact(ie) = iloc - loc0(ie)
    end do

    if (allocated(f%cs)) then
       ! orbital route: psi(n,:) = sqrt(f_n) sum_i c_in p(i,:)
       if (f%nocc == 0) return
       allocate(psi(f%nocc,nd))
       psi = 0d0
       do ie = 1, nenvl
          call blockmul('T',f%nocc,nd,nact(ie),f%cs(f%ibas0(iat(ie))+1,1),f%nbas,&
             p(loc0(ie)+1,1),nloc,1d0,psi,f%nocc)
       end do
       if (nsw == 2) then
          ! spin up and down from their own orbitals; the total is the sum
          n = f%noccs(1)
          call accum_props(psi(1:n,:),psi(1:n,:),nder,rho(2),grad(:,2),h(:,:,2),gkin(2))
          call accum_props(psi(n+1:,:),psi(n+1:,:),nder,rho(3),grad(:,3),h(:,:,3),gkin(3))
          rho(1) = rho(2) + rho(3)
          grad(:,1) = grad(:,2) + grad(:,3)
          h(:,:,1) = h(:,:,2) + h(:,:,3)
          gkin(1) = gkin(2) + gkin(3)
       else
          call accum_props(psi,psi,nder,rho(1),grad(:,1),h(:,:,1),gkin(1))
       end if
       return
    end if

    ! the pairs of environment atoms that share a block: the pair
    ! (ir,ic) reads the stored block as rows of ir and columns of ic;
    ! pairs whose contribution, bounded by |D|max sum|p_a| sum|p_b|,
    ! is below threshold are skipped (only the values and gradients
    ! enter the products)
    ndw = min(nd,4)
    allocate(w(nloc,ndw,nsw),psum(nenvl))
    w = 0d0
    do ie = 1, nenvl
       psum(ie) = sum(abs(p(loc0(ie)+1:loc0(ie)+nact(ie),1:ndw)))
    end do
    allocate(pib(nenvl*(nenvl+1)/2),pir(nenvl*(nenvl+1)/2),pic(nenvl*(nenvl+1)/2),pboth(nenvl*(nenvl+1)/2))
    np = 0
    do ie = 1, nenvl
       do je = ie, nenvl
          if (norm2(xenv(:,ie)-xenv(:,je)) > rc(ie) + rc(je)) cycle
          dt = tvec(:,ie) - tvec(:,je)
          ib = find_block(f,iat(ie),iat(je),dt,trans)
          if (ib == 0) cycle
          if (maxval(f%blknorm(ib,1:nsw)) * psum(ie) * psum(je) < dm_pairthr) cycle
          np = np + 1
          pib(np) = ib
          if (trans) then
             pir(np) = je
             pic(np) = ie
          else
             pir(np) = ie
             pic(np) = je
          end if
          pboth(np) = (ie /= je)
       end do
    end do
    if (np == 0) return

    ! contract with the density matrix, w = D p, one stored block at a
    ! time for all the pairs that use it (in a small cell, several),
    ! on the leading sub-block of the functions that are non-zero at
    ! the point; the transposed product gives the (ic,ir) contribution.
    ! A block used by one pair multiplies p into w directly; several
    ! pairs gather their columns in g and scatter the product from wk.
    allocate(iord(np))
    do k = 1, np
       iord(k) = k
    end do
    call qcksort(pib,iord,1,np)
    maxrun = 1
    k1 = 1
    do k = 2, np
       if (pib(iord(k)) /= pib(iord(k1))) k1 = k
       maxrun = max(maxrun,k - k1 + 1)
    end do
    nbmax = maxval(f%spc(:)%nbas)
    if (maxrun > 1) allocate(g(nbmax,maxrun*ndw),wk(nbmax,maxrun*ndw))
    k1 = 1
    do while (k1 <= np)
       ib = pib(iord(k1))
       k2 = k1
       do while (k2 < np)
          if (pib(iord(k2+1)) /= ib) exit
          k2 = k2 + 1
       end do
       ir = pir(iord(k1))
       ic = pic(iord(k1))
       nr = loc0(ir+1) - loc0(ir)
       if (k2 == k1) then
          mr = nact(ir)
          mc = nact(ic)
          do js = 1, nsw
             call blockmul('N',mr,ndw,mc,f%dm(f%blkoff(ib)+1,js),nr,p(loc0(ic)+1,1),nloc,1d0,&
                w(loc0(ir)+1,1,js),nloc)
             if (pboth(iord(k1))) &
                call blockmul('T',mc,ndw,mr,f%dm(f%blkoff(ib)+1,js),nr,p(loc0(ir)+1,1),nloc,1d0,&
                w(loc0(ic)+1,1,js),nloc)
          end do
          k1 = k2 + 1
          cycle
       end if
       mr = 0
       mc = 0
       do k = k1, k2
          mr = max(mr,nact(pir(iord(k))))
          mc = max(mc,nact(pic(iord(k))))
       end do

       ! rows of ir: D * p(ic)
       n = 0
       do k = k1, k2
          ic = pic(iord(k))
          g(1:mc,n+1:n+ndw) = p(loc0(ic)+1:loc0(ic)+mc,1:ndw)
          n = n + ndw
       end do
       do js = 1, nsw
          call blockmul('N',mr,n,mc,f%dm(f%blkoff(ib)+1,js),nr,g,nbmax,0d0,wk,nbmax)
          n2 = 0
          do k = k1, k2
             ir = pir(iord(k))
             w(loc0(ir)+1:loc0(ir)+mr,:,js) = w(loc0(ir)+1:loc0(ir)+mr,:,js) + wk(1:mr,n2+1:n2+ndw)
             n2 = n2 + ndw
          end do
       end do

       ! columns of ic: D^T * p(ir), for the pairs of two different atoms
       n = 0
       do k = k1, k2
          if (.not.pboth(iord(k))) cycle
          ir = pir(iord(k))
          g(1:mr,n+1:n+ndw) = p(loc0(ir)+1:loc0(ir)+mr,1:ndw)
          n = n + ndw
       end do
       if (n > 0) then
          do js = 1, nsw
             call blockmul('T',mc,n,mr,f%dm(f%blkoff(ib)+1,js),nr,g,nbmax,0d0,wk,nbmax)
             n2 = 0
             do k = k1, k2
                if (.not.pboth(iord(k))) cycle
                ic = pic(iord(k))
                w(loc0(ic)+1:loc0(ic)+mc,:,js) = w(loc0(ic)+1:loc0(ic)+mc,:,js) + wk(1:mc,n2+1:n2+ndw)
                n2 = n2 + ndw
             end do
          end do
       end if
       k1 = k2 + 1
    end do

    do js = 1, nsw
       call accum_props(p,w(:,:,js),nder,rho(js),grad(:,js),h(:,:,js),gkin(js))
    end do
    call spin_fill(nsw,rho,grad,h,gkin)

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
    integer, allocatable :: fsp(:), fl(:), nfnsp(:), ifn(:), iord(:), off(:)
    real*8, allocatable :: frc(:)
    real*8 :: rdum, rinc
    type(fhiaimsradial), allocatable :: fsort(:)

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

    ! the radial functions of each species by decreasing outer radius
    ! (ties in file order), so the functions that are non-zero at a
    ! given distance from the atom are a prefix; ifile maps the local
    ! basis functions back to the order of the file (the eigenvectors)
    do is = 1, f%nspc
       associate (s => f%spc(is))
         allocate(iord(s%nfn),off(s%nfn))
         off(1) = 0
         do j = 2, s%nfn
            off(j) = off(j-1) + 2 * s%fn(j-1)%l + 1
         end do
         do j = 1, s%nfn
            k = j
            do i = j-1, 1, -1
               if (s%fn(iord(i))%rcut >= s%fn(j)%rcut) exit
               iord(i+1) = iord(i)
               k = i
            end do
            iord(k) = j
         end do
         fsort = s%fn(iord)
         call move_alloc(fsort,s%fn)
         allocate(s%ifile(sum(2 * s%fn(:)%l + 1)))
         i = 0
         do k = 1, s%nfn
            do j = 1, 2 * s%fn(k)%l + 1
               i = i + 1
               s%ifile(i) = off(iord(k)) + j
            end do
         end do
         deallocate(iord,off)
       end associate
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
    allocate(f%dm(ntot,f%nspin))
    f%dm = 0d0

  end subroutine build_blocks

  !> Open file as a byte stream and position the buffer at the start
  !> of the eigenvector block: the line after the first one that
  !> contains "eigenvectors" (the real/complex flag line). errmsg is
  !> empty on success.
  subroutine sb_open(sb,file,errmsg)
    type(streambuf), intent(inout) :: sb
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: ios, m, j, k
    integer*8 :: fpos

    errmsg = "Error locating the eigenvectors in the dgrid file: " // trim(file)
    if (sb%lu >= 0) close(sb%lu)
    if (.not.allocated(sb%buf)) allocate(character(len=sb_size+1) :: sb%buf)
    inquire(file=file,size=sb%fsize)
    if (sb%fsize < 0) return
    open(newunit=sb%lu,file=file,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if (ios /= 0) then
       sb%lu = -1
       return
    end if

    ! find the marker; chunks overlap by 64 bytes so it is never split
    fpos = 1
    do while (fpos <= sb%fsize)
       m = int(min(int(sb_size,8),sb%fsize - fpos + 1))
       read (sb%lu,pos=fpos,iostat=ios) sb%buf(1:m)
       if (ios /= 0) return
       j = index(sb%buf(1:m),"eigenvectors")
       if (j > 0) then
          k = index(sb%buf(j:m),char(10))
          if (k == 0) then
             ! the line runs past the chunk: reread from the marker
             if (m < sb_size) return
             fpos = fpos + j - 1
             cycle
          end if
          sb%fpos = fpos + (j + k - 1)
          sb%n = 0
          sb%ip = 1
          sb%eof = .false.
          call sb_fill(sb,ios)
          if (ios /= 0) return
          errmsg = ""
          return
       end if
       if (m < sb_size) exit
       fpos = fpos + m - 64
    end do

  end subroutine sb_open

  !> Move the unread part of the buffer to its start and read more
  !> bytes from the file after it. ios is non-zero on error.
  subroutine sb_fill(sb,ios)
    type(streambuf), intent(inout) :: sb
    integer, intent(out) :: ios

    integer :: m

    ios = 0
    m = sb%n - sb%ip + 1
    if (m > 0) sb%buf(1:m) = sb%buf(sb%ip:sb%n)
    sb%n = max(m,0)
    sb%ip = 1
    sb%buf(sb%n+1:sb%n+1) = " "
    m = int(min(int(sb_size - sb%n,8),sb%fsize - sb%fpos + 1))
    if (m <= 0) then
       sb%eof = .true.
       return
    end if
    read (sb%lu,pos=sb%fpos,iostat=ios) sb%buf(sb%n+1:sb%n+m)
    if (ios /= 0) return
    sb%n = sb%n + m
    sb%fpos = sb%fpos + m
    sb%eof = (sb%fpos > sb%fsize)
    sb%buf(sb%n+1:sb%n+1) = " "

  end subroutine sb_fill

  !> Read the next nval blank-separated real values of the stream into
  !> x or, if x is absent, skip them. The values must end at a record
  !> boundary (FHI-aims writes one state per write statement). ios is
  !> non-zero on a malformed value, a missing record boundary or end
  !> of file. The conversion is C strtod on the buffer (a token is
  !> always followed by a blank, the sentinel at worst), which assumes
  !> the C numeric locale (nothing in critic2 calls setlocale); a
  !> token strtod does not consume whole (a D exponent, a three-digit
  !> exponent without the letter) falls back to a list-directed read.
  subroutine sb_values(sb,nval,ios,x)
    use iso_c_binding, only: c_ptr, c_loc, c_intptr_t, c_char, c_double
    type(streambuf), intent(inout), target :: sb
    integer, intent(in) :: nval
    integer, intent(out) :: ios
    real*8, intent(out), optional :: x(:)

    interface
       function strtod(str,endptr) bind(C,name="strtod")
         import c_char, c_ptr, c_double
         character(kind=c_char), intent(in) :: str(*)
         type(c_ptr), intent(out) :: endptr
         real(c_double) :: strtod
       end function strtod
    end interface

    integer, parameter :: toklen = 63

    integer :: nv, i0, i1, j, k
    type(c_ptr) :: endp

    ios = 0
    if (nval == 0) return
    nv = 0
    do while (nv < nval)
       ! skip blanks
       do
          if (sb%ip > sb%n) then
             call sb_fill(sb,ios)
             if (ios /= 0) return
             if (sb%n == 0) goto 999
          end if
          k = verify(sb%buf(sb%ip:sb%n),blanks)
          if (k > 0) exit
          sb%ip = sb%n + 1
       end do
       sb%ip = sb%ip + k - 1

       ! the token, whole in the buffer
       if (sb%n - sb%ip < toklen .and. .not.sb%eof) then
          call sb_fill(sb,ios)
          if (ios /= 0) return
       end if
       i0 = sb%ip
       k = scan(sb%buf(i0:min(i0+toklen,sb%n)),blanks)
       if (k > 0) then
          i1 = i0 + k - 2
       elseif (i0 + toklen > sb%n) then
          i1 = sb%n
       else
          goto 999
       end if
       sb%ip = i1 + 1
       nv = nv + 1
       if (present(x)) then
          x(nv) = strtod(sb%buf(i0:sb%n+1),endp)
          if (transfer(endp,0_c_intptr_t) - transfer(c_loc(sb%buf(i0:i0)),0_c_intptr_t) /= i1 - i0 + 1) then
             read (sb%buf(i0:i1),*,iostat=ios) x(nv)
             if (ios /= 0) goto 999
          end if
       end if
    end do

    ! only blanks up to the end of the record (or of the file)
    do
       if (sb%ip > sb%n) then
          call sb_fill(sb,ios)
          if (ios /= 0) return
          if (sb%n == 0) return
       end if
       j = index(sb%buf(sb%ip:sb%n),char(10))
       k = verify(sb%buf(sb%ip:sb%n),blanks)
       if (k > 0 .and. (j == 0 .or. k < j)) goto 999
       if (j > 0) then
          sb%ip = sb%ip + j
          return
       end if
       sb%ip = sb%n + 1
    end do

999 continue
    ios = 1

  end subroutine sb_values

  !> Accumulate the contribution of one spin channel and k-point
  !> (fractional kpt, weight wk) to the real-space density matrix
  !> blocks, from the eigenvectors (evk (nbas,nsv) if complex, evr if
  !> real) and the occupations (occ, not k-weighted):
  !>   D_ij(dT) += Re wk sum_n occ_n c_in c_jn^* exp(2 pi i k.dT)
  !> The real part is exact: the imaginary parts cancel in the
  !> symmetric sum over basis function pairs. Only the columns of the
  !> atoms b >= a that share a block with a are computed. The
  !> contribution goes to the total matrix and, if the calculation is
  !> spin-polarized and is == 1, to the spin-up matrix as well.
  subroutine accumulate_dm(f,evk,evr,occ,kpt,wk,nsv,is)
    use param, only: tpi, img
    type(fhiaimswfn), intent(inout) :: f
    complex*16, intent(in) :: evk(:,:)
    real*8, intent(in) :: evr(:,:)
    real*8, intent(in) :: occ(:)
    real*8, intent(in) :: kpt(3)
    real*8, intent(in) :: wk
    integer, intent(in) :: nsv
    integer, intent(in) :: is

    integer :: a, b, na, nb, a0, b0, n, ib, ncol, js, nj
    integer*8 :: off
    integer, allocatable :: col0(:)
    complex*16, allocatable :: fw(:,:), fwa(:,:), en(:,:), pa(:,:)
    real*8, allocatable :: fwr(:,:), fwar(:,:), enr(:,:), par(:,:), blk(:)
    complex*16 :: phase
    real*8 :: arg

    complex*16, parameter :: cone = (1d0,0d0), czero = (0d0,0d0)

    if (nsv == 0) return
    nj = 1
    if (f%nspin == 2 .and. is == 1) nj = 2

    ! occupation- and weight-scaled eigenvectors
    allocate(col0(f%nat),blk(maxval(f%spc(:)%nbas)**2))
    if (f%isreal) then
       allocate(fwr(f%nbas,nsv))
       do n = 1, nsv
          fwr(:,n) = evr(:,n) * (occ(n) * wk)
       end do
    else
       allocate(fw(f%nbas,nsv))
       do n = 1, nsv
          fw(:,n) = evk(:,n) * (occ(n) * wk)
       end do
    end if

    do a = 1, f%nat
       na = f%spc(f%ispc(a))%nbas
       a0 = f%ibas0(a)

       ! gather the eigenvector rows of a and of its neighbours, in the
       ! local (sorted) order of the basis functions
       ncol = 0
       do b = a, f%nat
          if (f%pairn(a,b) == 0) cycle
          col0(b) = ncol
          ncol = ncol + f%spc(f%ispc(b))%nbas
       end do
       if (f%isreal) then
          if (allocated(fwar)) deallocate(fwar,enr,par)
          allocate(fwar(na,nsv),enr(ncol,nsv),par(na,ncol))
          fwar = fwr(a0 + f%spc(f%ispc(a))%ifile,:)
       else
          if (allocated(fwa)) deallocate(fwa,en,pa)
          allocate(fwa(na,nsv),en(ncol,nsv),pa(na,ncol))
          fwa = fw(a0 + f%spc(f%ispc(a))%ifile,:)
       end if
       do b = a, f%nat
          if (f%pairn(a,b) == 0) cycle
          nb = f%spc(f%ispc(b))%nbas
          b0 = f%ibas0(b)
          if (f%isreal) then
             enr(col0(b)+1:col0(b)+nb,:) = evr(b0 + f%spc(f%ispc(b))%ifile,1:nsv)
          else
             en(col0(b)+1:col0(b)+nb,:) = evk(b0 + f%spc(f%ispc(b))%ifile,1:nsv)
          end if
       end do

       ! the a-rows of the k-space density matrix, neighbour columns
       if (f%isreal) then
          call dgemm('N','T',na,ncol,nsv,1d0,fwar,na,enr,ncol,0d0,par,na)
       else
          call zgemm('N','C',na,ncol,nsv,cone,fwa,na,en,ncol,czero,pa,na)
       end if

       do b = a, f%nat
          if (f%pairn(a,b) == 0) cycle
          nb = f%spc(f%ispc(b))%nbas
          n = na * nb
          do ib = f%pair0(a,b), f%pair0(a,b) + f%pairn(a,b) - 1
             arg = tpi * dot_product(kpt,real(f%blkt(:,ib),8))
             off = f%blkoff(ib)
             if (f%isreal) then
                blk(1:n) = cos(arg) * reshape(par(:,col0(b)+1:col0(b)+nb),(/n/))
             else
                phase = exp(img * arg)
                blk(1:n) = reshape(real(pa(:,col0(b)+1:col0(b)+nb) * phase,8),(/n/))
             end if
             do js = 1, nj
                f%dm(off+1:off+n,js) = f%dm(off+1:off+n,js) + blk(1:n)
             end do
          end do
       end do
    end do

  end subroutine accumulate_dm

  !> Largest number of local basis functions (those of the atoms
  !> within their cutoff) at an atomic position, the densest points:
  !> an estimate of nloc in rho2.
  function max_nloc(f,c) result(nloc)
    use crystalmod, only: crystal
    use param, only: icrd_cart
    type(fhiaimswfn), intent(in) :: f
    type(crystal), intent(inout) :: c
    integer :: nloc

    integer :: e, ion, nenv, ia, n
    integer, allocatable :: eid(:), lvec(:,:)
    real*8, allocatable :: dist(:)

    nloc = 0
    do e = 1, c%ncel
       call c%list_near_atoms(c%atcel(e)%r,icrd_cart,.false.,nenv,eid,dist,lvec,up2dsp=f%spcutoff)
       n = 0
       do ion = 1, nenv
          ia = f%iaims(eid(ion))
          if (dist(ion) <= f%spc(f%ispc(ia))%rcut) n = n + f%spc(f%ispc(ia))%nbas
       end do
       nloc = max(nloc,n)
    end do

  end function max_nloc

  !> Largest element of each density matrix block, for the pair
  !> screening in rho2.
  subroutine block_norms(f)
    type(fhiaimswfn), intent(inout) :: f

    integer :: a, b, ib, n, js
    integer*8 :: off

    if (allocated(f%blknorm)) deallocate(f%blknorm)
    allocate(f%blknorm(f%nblk,f%nspin))
    do a = 1, f%nat
       do b = a, f%nat
          n = f%spc(f%ispc(a))%nbas * f%spc(f%ispc(b))%nbas
          do ib = f%pair0(a,b), f%pair0(a,b) + f%pairn(a,b) - 1
             off = f%blkoff(ib)
             do js = 1, f%nspin
                f%blknorm(ib,js) = maxval(abs(f%dm(off+1:off+n,js)))
             end do
          end do
       end do
    end do

  end subroutine block_norms

  !> C = op(A) B + beta C, with op(A) = A (op = 'N', A is m x k) or
  !> A^T (op = 'T', A is k x m) and B k x n. A single column goes
  !> through dgemv: dgemm spends most of its time packing A when n
  !> is small (per-column dgemv is slower again from n = 4).
  subroutine blockmul(op,m,n,k,a,lda,b,ldb,beta,c,ldc)
    character*1, intent(in) :: op
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real*8, intent(in) :: a(lda,*), b(ldb,*)
    real*8, intent(in) :: beta
    real*8, intent(inout) :: c(ldc,*)

    if (n > 1) then
       call dgemm(op,'N',m,n,k,1d0,a,lda,b,ldb,beta,c,ldc)
    elseif (op == 'N') then
       call dgemv('N',m,k,1d0,a,lda,b,1,beta,c,1)
    else
       call dgemv('T',k,m,1d0,a,lda,b,1,beta,c,1)
    end if

  end subroutine blockmul

  !> Density and derivatives from two sets of columns x and y, where
  !> column 1 is the value, 2:4 the gradient and 5:10 the second
  !> derivatives (xx,xy,xz,yy,yz,zz), of basis functions and their
  !> density matrix contraction (x = p, y = D p) or of the orbitals
  !> (x = y = psi): rho = x1.y1, grad = 2 x_i.y1, G = 1/2 x_i.y_i,
  !> H_ij = 2 (x_ij.y1 + x_i.y_j).
  subroutine accum_props(x,y,nder,rho,grad,h,gkin)
    real*8, intent(in) :: x(:,:), y(:,:)
    integer, intent(in) :: nder
    real*8, intent(out) :: rho, grad(3), h(3,3), gkin

    integer :: i, j, k

    rho = dot_product(x(:,1),y(:,1))
    grad = 0d0
    h = 0d0
    gkin = 0d0
    if (nder >= 1) then
       do i = 1, 3
          grad(i) = 2d0 * dot_product(x(:,1+i),y(:,1))
          gkin = gkin + 0.5d0 * dot_product(x(:,1+i),y(:,1+i))
       end do
    end if
    if (nder >= 2) then
       k = 4
       do i = 1, 3
          do j = i, 3
             k = k + 1
             h(i,j) = 2d0 * (dot_product(x(:,k),y(:,1)) + dot_product(x(:,1+i),y(:,1+j)))
             h(j,i) = h(i,j)
          end do
       end do
    end if

  end subroutine accum_props

  !> Complete the spin components of the rho2 output when the spin-up
  !> part was computed (nsw == 2): spin down (index 3) is the total
  !> minus spin up.
  subroutine spin_fill(nsw,rho,grad,h,gkin)
    integer, intent(in) :: nsw
    real*8, intent(inout) :: rho(3), grad(3,3), h(3,3,3), gkin(3)

    if (nsw == 2) then
       rho(3) = rho(1) - rho(2)
       grad(:,3) = grad(:,1) - grad(:,2)
       h(:,:,3) = h(:,:,1) - h(:,:,2)
       gkin(3) = gkin(1) - gkin(2)
    end if

  end subroutine spin_fill

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
