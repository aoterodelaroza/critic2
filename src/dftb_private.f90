! Copyright (c) 2016 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

! Interface to DFTB+ wavefunctions.
module dftb_private
  implicit none
  
  private

  public :: dftb_read
  public :: dftb_rho2
  public :: dftb_register_struct
  private :: next_logical
  private :: next_integer
  private :: read_kpointsandweights
  private :: read_occupations
  private :: dftb_read_reals1
  private :: next_hsd_atom
  private :: build_interpolation_grid1
  private :: calculate_rl

  ! private structural info
  real*8, allocatable :: renv(:,:)
  integer, allocatable :: lenv(:,:), idxenv(:), zenv(:)
  integer :: nenv
  real*8 :: globalcutoff = 0d0

  ! minimum distance (bohr)
  real*8, parameter :: mindist = 1d-6

contains

  !> Read the information for a DFTB+ field from the detailed.xml,
  !> eigenvec.bin, and the basis set definition in HSD format.
  subroutine dftb_read(f,filexml,filebin,filehsd,zcel)
    use types, only: field, dftbatom, realloc
    use tools_io, only: fopen_read, getline_raw, lower, ferror, faterr, string, fclose
    use param, only: tpi, maxzat0
    
    type(field), intent(out) :: f !< Output field
    character*(*), intent(in) :: filexml !< The detailed.xml file
    character*(*), intent(in) :: filebin !< The eigenvec.bin file
    character*(*), intent(in) :: filehsd !< The definition of the basis set in hsd format
    integer, intent(in) :: zcel(:) !< atomic numbers for the atoms in the complete cell (in order)

    integer :: lu, i, j, k, idum, n, nat, id
    character(len=:), allocatable :: line
    logical :: ok, iread(5)
    type(dftbatom) :: at
    real*8, allocatable :: dw(:)

    ! number of atoms in the cell
    nat = size(zcel)

    ! detailed.xml, first pass
    lu = fopen_read(filexml)
    iread = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = adjustl(lower(line))
       if (index(line,"<real>") > 0) then
          f%isreal = next_logical(lu,line,"real")
          iread(1) = .true.
       elseif (index(line,"<nrofkpoints>") > 0) then
          f%nkpt = next_integer(lu,line,"nrofkpoints")
          iread(2) = .true.
       elseif (index(line,"<nrofspins>") > 0) then
          f%nspin = next_integer(lu,line,"nrofspins")
          iread(3) = .true.
       elseif (index(line,"<nrofstates>") > 0) then
          f%nstates = next_integer(lu,line,"nrofstates")
          iread(4) = .true.
       elseif (index(line,"<nroforbitals>") > 0) then
          f%norb = next_integer(lu,line,"nroforbitals")
          iread(5) = .true.
       end if
    end do
    if (.not.all(iread)) &
       call ferror('dftb_read','missing information in xml',faterr)

    ! detailed.xml, second pass
    if (allocated(f%dkpt)) deallocate(f%dkpt)
    if (allocated(f%docc)) deallocate(f%docc)
    allocate(f%dkpt(3,f%nkpt),dw(f%nkpt),f%docc(f%nstates,f%nkpt,f%nspin))
    rewind(lu)
    call read_kpointsandweights(lu,f%dkpt,dw)
    call read_occupations(lu,f%docc)
    call fclose(lu)

    ! use occupations times weights; scale the kpts by 2pi
    do i = 1, f%nkpt
       f%docc(:,i,:) = f%docc(:,i,:) * dw(i)
    end do
    f%dkpt = f%dkpt * tpi
    deallocate(dw)

    ! eigenvec.bin
    lu = fopen_read(filebin,"unformatted")
    read (lu) idum ! this is the identity number

    ! read the eigenvectors
    if (f%isreal) then
       if (allocated(f%evecr)) deallocate(f%evecr)
       allocate(f%evecr(f%norb,f%nstates,f%nspin))
       do i = 1, f%nspin
          do k = 1, f%nstates
             read (lu) f%evecr(:,k,i)
          end do
       end do
    else
       if (allocated(f%evecc)) deallocate(f%evecc)
       allocate(f%evecc(f%norb,f%nstates,f%nkpt,f%nspin))
       do i = 1, f%nspin
          do j = 1, f%nkpt
             do k = 1, f%nstates
                read (lu) f%evecc(:,k,j,i)
             end do
          end do
       end do
    end if
    call fclose(lu)

    ! open the hsd file with the basis definition
    lu = fopen_read(filehsd)
    if (allocated(f%bas)) deallocate(f%bas)
    allocate(f%bas(10))
    n = 0
    do while(next_hsd_atom(lu,at))
       if (.not.any(zcel == at%z)) cycle
       n = n + 1 
       if (n > size(f%bas)) call realloc(f%bas,2*n)
       f%bas(n) = at
    end do
    call realloc(f%bas,n)
    call fclose(lu)

    ! tie the atomic numbers to the basis types
    if (allocated(f%ispec)) deallocate(f%ispec)
    allocate(f%ispec(maxzat0))
    f%ispec = 0
    do i = 1, maxzat0
       do j = 1, n
          if (f%bas(j)%z == i) then
             f%ispec(i) = j
             exit
          end if
       end do
    end do

    ! indices for the atomic orbitals, for array sizes later on
    if (allocated(f%idxorb)) deallocate(f%idxorb)
    allocate(f%idxorb(nat))
    n = 0
    f%maxnorb = 0
    f%maxlm = 0
    do i = 1, nat
       id = f%ispec(zcel(i))
       if (id == 0) call ferror('dftb_read','basis missing for atomic number ' // string(zcel(i)),faterr)

       f%maxnorb = max(f%maxnorb,f%bas(id)%norb)
       f%idxorb(i) = n + 1
       do j = 1, f%bas(id)%norb
          n = n + 2*f%bas(id)%l(j) + 1
          f%maxlm = max(f%maxlm,f%bas(id)%l(j))
       end do
    end do
    f%midxorb = n
    f%maxlm = (f%maxlm+3)*(f%maxlm+3)

    call build_interpolation_grid1(f)

    ! finished
    f%init = .true.

  end subroutine dftb_read

  !> Calculate the density and derivatives of a DFTB+ field (f) up to
  !> the nder degree (max = 2). xpos is in Cartesian coordinates. In
  !> output, the density (rho), the gradient (grad), the Hessian (h),
  !> and the G(r) kinetic energy density (gkin).  This routine is
  !> thread-safe.
  subroutine dftb_rho2(f,xpos,nder,rho,grad,h,gkin)
    use grid1_tools
    use tools_math
    use tools_io
    use types
    use param

    type(field), intent(inout) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    integer, intent(in) :: nder  !< Number of derivatives
    real*8, intent(out) :: rho !< Density
    real*8, intent(out) :: grad(3) !< Gradient
    real*8, intent(out) :: h(3,3) !< Hessian 
    real*8, intent(out) :: gkin !< G(r), kinetic energy density

    integer, parameter :: maxenvl = 50

    integer :: ion, it, is, istate, ik, iorb, i, l, m, lmax
    integer :: ixorb, ixorb0
    real*8 :: xion(3), rcut, dist, r, tp(2)
    complex*16 :: xao(f%midxorb), xaol(f%midxorb), xmo, xmop(3), xmopp(6)
    complex*16 :: xaolp(3,f%midxorb), xaolpp(6,f%midxorb)
    complex*16 :: xaop(3,f%midxorb), xaopp(6,f%midxorb)
    real*8 :: rao(f%midxorb), raol(f%midxorb), rmo, rmop(3), rmopp(6)
    real*8 :: raolp(3,f%midxorb), raolpp(6,f%midxorb)
    real*8 :: raop(3,f%midxorb), raopp(6,f%midxorb)
    integer :: nenvl, idxion(maxenvl), ionl
    real*8 :: rl(f%maxnorb,maxenvl), rlp(f%maxnorb,maxenvl), rlpp(f%maxnorb,maxenvl)
    complex*16 :: phase(maxenvl,f%nkpt), ylm(f%maxlm)
    integer :: imin, ip, im, iphas
    real*8 :: phi(f%maxlm,f%maxnorb,maxenvl)
    real*8 :: phip(3,f%maxlm,f%maxnorb,maxenvl)
    real*8 :: phipp(6,f%maxlm,f%maxnorb,maxenvl)
    complex*16 :: xgrad1(3), xgrad2(3), xhess1(6), xhess2(6)

    ! precalculate the quantities that depend only on the environment
    nenvl = 0
    phi = 0d0
    phip = 0d0
    phipp = 0d0
    do ion = 1, nenv
       xion = xpos - renv(:,ion)
       it = f%ispec(zenv(ion))

       ! apply the distance cutoff
       rcut = maxval(f%bas(it)%cutoff(1:f%bas(it)%norb))
       if (any(abs(xion) > rcut)) cycle
       dist = sqrt(dot_product(xion,xion))
       if (dist > rcut) cycle

       ! write down this atom
       nenvl = nenvl + 1
       if (nenvl > maxenvl) call ferror('dftb_rho2','local environment exceeded array size',faterr)
       idxion(nenvl) = ion

       ! calculate the spherical harmonics contributions for this atom
       lmax = maxval(f%bas(it)%l(1:f%bas(it)%norb))
       call tosphere(xion,r,tp)
       r = max(r,mindist)
       call genylm(lmax+2,tp,ylm)

       ! calculate the radial contributions for this atom
       do iorb = 1, f%bas(it)%norb
          if (dist > f%bas(it)%cutoff(iorb)) cycle
          if (f%exact) then
             call calculate_rl(f,it,iorb,dist,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl))
          else
             call grid1_interp(f%bas(it)%orb(iorb),dist,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl))
          end if
       end do

       ! populate the phi array and calculate the derivatives
       do iorb = 1, f%bas(it)%norb
          l = f%bas(it)%l(iorb)
          imin = l*l+1

          ! m = 0
          ip = imin+l
          call ylmderiv(ylm,r,l,0,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl),xgrad1,xhess1)
          phi(ip,iorb,nenvl) = rl(iorb,nenvl) * real(ylm(ip),8)
          phip(:,ip,iorb,nenvl) = real(xgrad1,8)
          phipp(:,ip,iorb,nenvl) = real(xhess1,8)

          ! |m| > 0
          do m = 1, l
             ip = imin + l + m
             im = imin + l - m
             iphas = (-1)**m
             call ylmderiv(ylm,r,l, m,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl),xgrad1,xhess1)
             call ylmderiv(ylm,r,l,-m,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl),xgrad2,xhess2)

             phi(ip,iorb,nenvl) = rl(iorb,nenvl) * real(iphas*ylm(ip)+ylm(im),8) / sqrt(2d0)
             phip(:,ip,iorb,nenvl) = real(iphas*xgrad1+xgrad2,8) / sqrt(2d0)
             phipp(:,ip,iorb,nenvl) = real(iphas*xhess1+xhess2,8) / sqrt(2d0)

             phi(im,iorb,nenvl) = rl(iorb,nenvl) * real(-iphas*img*ylm(ip)+img*ylm(im),8) / sqrt(2d0)
             phip(:,im,iorb,nenvl) = real(-iphas*img*xgrad1+img*xgrad2,8) / sqrt(2d0)
             phipp(:,im,iorb,nenvl) = real(-iphas*img*xhess1+img*xhess2,8) / sqrt(2d0)
          end do
       end do

       ! calculate the phases
       if (.not.f%isreal) then
          do ik = 1, f%nkpt
             phase(nenvl,ik) = exp(img * dot_product(lenv(:,ion),f%dkpt(:,ik)))
          end do
       end if
    end do

    ! initialize
    rho = 0d0
    grad = 0d0
    h = 0d0
    gkin = 0d0

    if (.not.f%isreal) then
       ! complex version; for crystals with non-gamma k-points
       ! run over spins
       do is = 1, f%nspin
          ! run over k-points
          do ik = 1, f%nkpt
             ! run over the states
             do istate = 1, f%nstates

                ! determine the atomic contributions. The unit cell atoms
                ! in critic are in the same order as in dftb+, which is
                ! why indexing the evecc/evecr works.
                xao = 0d0
                xaop = 0d0
                xaopp = 0d0
                ! run over atoms
                do ionl = 1, nenvl
                   ion = idxion(ionl)
                   it = f%ispec(zenv(ion))

                   ! run over atomic orbitals
                   ixorb0 = f%idxorb(idxenv(ion))
                   ixorb = ixorb0 - 1
                   do iorb = 1, f%bas(it)%norb
                      ! run over ms for the same l
                      do i = f%bas(it)%l(iorb)*f%bas(it)%l(iorb)+1, (f%bas(it)%l(iorb)+1)*(f%bas(it)%l(iorb)+1)
                         ixorb = ixorb + 1
                         xaol(ixorb) = phi(i,iorb,ionl)
                         xaolp(:,ixorb) = phip(:,i,iorb,ionl)
                         xaolpp(:,ixorb) = phipp(:,i,iorb,ionl)
                      end do
                   end do ! iorb

                   xao(ixorb0:ixorb) = xao(ixorb0:ixorb) + xaol(ixorb0:ixorb) * phase(ionl,ik)
                   xaop(:,ixorb0:ixorb) = xaop(:,ixorb0:ixorb) + xaolp(:,ixorb0:ixorb) * phase(ionl,ik)
                   xaopp(:,ixorb0:ixorb) = xaopp(:,ixorb0:ixorb) + xaolpp(:,ixorb0:ixorb) * phase(ionl,ik)
                end do ! ion

                ! calculate the value of this extended orbital and its derivatives
                xmo = 0d0
                xmop = 0d0
                xmopp = 0d0
                do i = 1, f%midxorb
                   xmo = xmo + conjg(xao(i))*f%evecc(i,istate,ik,is)
                   xmop = xmop + conjg(xaop(:,i))*f%evecc(i,istate,ik,is)
                   xmopp = xmopp + conjg(xaopp(:,i))*f%evecc(i,istate,ik,is)
                end do

                ! accumulate properties
                rho = rho + (conjg(xmo)*xmo) * f%docc(istate,ik,is)
                grad = grad + (conjg(xmop)*xmo+conjg(xmo)*xmop) * f%docc(istate,ik,is)
                h(1,1) = h(1,1) + (conjg(xmopp(1))*xmo+conjg(xmop(1))*xmop(1)+conjg(xmop(1))*xmop(1)+conjg(xmo)*xmopp(1)) * f%docc(istate,ik,is)
                h(1,2) = h(1,2) + (conjg(xmopp(2))*xmo+conjg(xmop(1))*xmop(2)+conjg(xmop(2))*xmop(1)+conjg(xmo)*xmopp(2)) * f%docc(istate,ik,is)
                h(1,3) = h(1,3) + (conjg(xmopp(3))*xmo+conjg(xmop(1))*xmop(3)+conjg(xmop(3))*xmop(1)+conjg(xmo)*xmopp(3)) * f%docc(istate,ik,is)
                h(2,2) = h(2,2) + (conjg(xmopp(4))*xmo+conjg(xmop(2))*xmop(2)+conjg(xmop(2))*xmop(2)+conjg(xmo)*xmopp(4)) * f%docc(istate,ik,is)
                h(2,3) = h(2,3) + (conjg(xmopp(5))*xmo+conjg(xmop(2))*xmop(3)+conjg(xmop(3))*xmop(2)+conjg(xmo)*xmopp(5)) * f%docc(istate,ik,is)
                h(3,3) = h(3,3) + (conjg(xmopp(6))*xmo+conjg(xmop(3))*xmop(3)+conjg(xmop(3))*xmop(3)+conjg(xmo)*xmopp(6)) * f%docc(istate,ik,is)
                gkin = gkin + (conjg(xmop(1))*xmop(1)+conjg(xmop(2))*xmop(2)+conjg(xmop(3))*xmop(3)) * f%docc(istate,ik,is)
             end do ! states
          end do ! k-points
       end do ! spins
    else
       ! real, for molecules or crystals with a single k-point at gamma
       ! run over spins
       do is = 1, f%nspin
          ! run over the states
          do istate = 1, f%nstates

             ! determine the atomic contributions. The unit cell atoms
             ! in critic are in the same order as in dftb+, which is
             ! why indexing the evecc/evecr works.
             rao = 0d0
             raop = 0d0
             raopp = 0d0
             ! run over atoms
             do ionl = 1, nenvl
                ion = idxion(ionl)
                it = f%ispec(zenv(ion))

                ! run over atomic orbitals
                ixorb0 = f%idxorb(idxenv(ion))
                ixorb = ixorb0 - 1
                do iorb = 1, f%bas(it)%norb
                   ! run over ms for the same l
                   do i = f%bas(it)%l(iorb)*f%bas(it)%l(iorb)+1, (f%bas(it)%l(iorb)+1)*(f%bas(it)%l(iorb)+1)
                      ixorb = ixorb + 1
                      raol(ixorb) = phi(i,iorb,ionl)
                      raolp(:,ixorb) = phip(:,i,iorb,ionl)
                      raolpp(:,ixorb) = phipp(:,i,iorb,ionl)
                   end do
                end do ! iorb

                rao(ixorb0:ixorb) = rao(ixorb0:ixorb) + raol(ixorb0:ixorb)
                raop(:,ixorb0:ixorb) = raop(:,ixorb0:ixorb) + raolp(:,ixorb0:ixorb)
                raopp(:,ixorb0:ixorb) = raopp(:,ixorb0:ixorb) + raolpp(:,ixorb0:ixorb)
             end do ! ion

             ! calculate the value of this extended orbital and its derivatives
             rmo = 0d0
             rmop = 0d0
             rmopp = 0d0
             do i = 1, f%midxorb
                rmo = rmo + rao(i)*f%evecr(i,istate,is)
                rmop = rmop + raop(:,i)*f%evecr(i,istate,is)
                rmopp = rmopp + raopp(:,i)*f%evecr(i,istate,is)
             end do

             ! accumulate properties
             rho = rho + (rmo*rmo) * f%docc(istate,1,is)
             grad = grad + (rmop*rmo+rmo*rmop) * f%docc(istate,1,is)
             h(1,1) = h(1,1) + (rmopp(1)*rmo+rmop(1)*rmop(1)+rmop(1)*rmop(1)+rmo*rmopp(1)) * f%docc(istate,1,is)
             h(1,2) = h(1,2) + (rmopp(2)*rmo+rmop(1)*rmop(2)+rmop(2)*rmop(1)+rmo*rmopp(2)) * f%docc(istate,1,is)
             h(1,3) = h(1,3) + (rmopp(3)*rmo+rmop(1)*rmop(3)+rmop(3)*rmop(1)+rmo*rmopp(3)) * f%docc(istate,1,is)
             h(2,2) = h(2,2) + (rmopp(4)*rmo+rmop(2)*rmop(2)+rmop(2)*rmop(2)+rmo*rmopp(4)) * f%docc(istate,1,is)
             h(2,3) = h(2,3) + (rmopp(5)*rmo+rmop(2)*rmop(3)+rmop(3)*rmop(2)+rmo*rmopp(5)) * f%docc(istate,1,is)
             h(3,3) = h(3,3) + (rmopp(6)*rmo+rmop(3)*rmop(3)+rmop(3)*rmop(3)+rmo*rmopp(6)) * f%docc(istate,1,is)
             gkin = gkin + (rmop(1)*rmop(1)+rmop(2)*rmop(2)+rmop(3)*rmop(3)) * f%docc(istate,1,is)
          end do ! states
       end do ! spins
    end if
    ! clean up
    h(2,1) = h(1,2)
    h(3,1) = h(1,3)
    h(3,2) = h(2,3)
    gkin = 0.50 * gkin

  end subroutine dftb_rho2

  !> Register structural information. rmat is the crys2car matrix,
  !> maxcutoff is the maximum orbital cutoff, nenv, renv, lenv, idx,
  !> and zenv is the environment information (number, position,
  !> lattice vector, index in the complete list, and atomic number.
  subroutine dftb_register_struct(rmat,maxcutoff,nenv0,renv0,lenv0,idx0,zenv0)
    use types
    use tools_math

    real*8, intent(in) :: rmat(3,3)
    real*8, intent(in) :: maxcutoff
    integer, intent(in) :: nenv0
    real*8, intent(in) :: renv0(:,:)
    integer, intent(in) :: lenv0(:,:), idx0(:), zenv0(:)

    real*8 :: sphmax, x0(3), dist
    integer :: i

    ! calculate the sphere radius that encompasses the unit cell
    if (maxcutoff > globalcutoff) then
       sphmax = norm(matmul(rmat,(/0d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/)))
       sphmax = max(sphmax,norm(matmul(rmat,(/1d0,0d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
       sphmax = max(sphmax,norm(matmul(rmat,(/0d0,1d0,0d0/) - (/0.5d0,0.5d0,0.5d0/))))
       sphmax = max(sphmax,norm(matmul(rmat,(/0d0,0d0,1d0/) - (/0.5d0,0.5d0,0.5d0/))))

       if (allocated(renv)) deallocate(renv)
       allocate(renv(3,size(renv0,2)))
       if (allocated(lenv)) deallocate(lenv)
       allocate(lenv(3,size(lenv0,2)))
       if (allocated(idxenv)) deallocate(idxenv)
       allocate(idxenv(size(idx0)))
       if (allocated(zenv)) deallocate(zenv)
       allocate(zenv(size(zenv0)))

       ! save the atomic environment
       nenv = 0
       x0 = matmul(rmat,(/0.5d0,0.5d0,0.5d0/))
       do i = 1, nenv0
          dist = norm(renv0(:,i) - x0)
          if (dist <= sphmax+maxcutoff) then
             nenv = nenv + 1
             renv(:,nenv) = renv0(:,i)
             lenv(:,nenv) = lenv0(:,i)
             idxenv(nenv) = idx0(i)
             zenv(nenv) = zenv0(i)
          end if
       end do
       globalcutoff = maxcutoff

       ! reallocate
       call realloc(renv,3,nenv)
       call realloc(lenv,3,nenv)
       call realloc(idxenv,nenv)
       call realloc(zenv,nenv)
    end if

  end subroutine dftb_register_struct

  !xx! private !xx! 

  !> Read the next logical value in xml format.
  function next_logical(lu,line0,key0) result(next)
    use tools_io
    integer, intent(in) :: lu
    character*(*), intent(in) :: line0, key0
    logical :: next

    character(len=:), allocatable :: line, key, aux
    integer :: idx, idxn, idxy
    logical :: found, ok, lastpass
    
    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    ! lowercase and build key
    key = "<" // trim(adjustl(key0)) // ">"
    line = lower(line0)

    ! delete up to the key in the current line, leave the rest
    idx = index(line,key)
    aux = line(idx+len_trim(key):)
    line = trim(adjustl(aux))

    ! parse all lines until </key> is found
    key = "</" // trim(adjustl(key0)) // ">"
    found = .false.
    lastpass = .false.
    next = .false.
    do while (.true.)
       ! are there any nos or yes? if so, which is first?
       idxn = index(line,"no")
       idxy = index(line,"yes")
       if (idxn > 0 .and. idxy == 0) then
          found = .true.
          next = .false.
       elseif (idxy > 0 .and. idxn == 0) then
          found = .true.
          next = .true.
       elseif (idxy > 0 .and. idxn > 0) then
          found = .true.
          next = (idxy < idxn)
       end if
       
       ! exit if found or if this was the last pass (after a </key> was read)
       if (found .or. lastpass) exit

       ! get a new line
       ok = getline_raw(lu,line)
       if (.not.ok) exit

       ! clean up and lowercase
       aux = lower(line)
       line = trim(adjustl(aux))

       ! if the </key> is detected, read up to the </key> and activate the last pass flag
       idx = index(line,key)
       if (idx > 0) then
          aux = line(1:idx-1)
          line = aux
          lastpass = .true.
       end if
    end do

    if (.not.found) &
       call ferror("next_logical","could not parse xml file",faterr)

  end function next_logical

  !> Read the next integer value in xml format.
  function next_integer(lu,line0,key0) result(next)
    use tools_io
    integer, intent(in) :: lu
    character*(*), intent(in) :: line0, key0
    integer :: next

    character(len=:), allocatable :: line, key, aux
    integer :: idx, lp 
    logical :: found, ok, lastpass
    
    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    ! lowercase and build key
    key = "<" // trim(adjustl(key0)) // ">"
    line = lower(line0)

    ! delete up to the key in the current line, leave the rest
    idx = index(line,key)
    aux = line(idx+len_trim(key):)
    line = trim(adjustl(aux))
    lp = 1

    ! parse all lines until </key> is found
    key = "</" // trim(adjustl(key0)) // ">"
    found = .false.
    lastpass = .false.
    next = 0
    do while (.true.)
       ! are there any integers in this line? 
       found = isinteger(next,line,lp)
       
       ! exit if found or if this was the last pass (after a </key> was read)
       if (found .or. lastpass) exit

       ! get a new line
       ok = getline_raw(lu,line)
       if (.not.ok) exit

       ! clean up and lowercase
       aux = lower(line)
       line = trim(adjustl(aux))

       ! if the </key> is detected, read up to the </key> and activate the last pass flag
       idx = index(line,key)
       if (idx > 0) then
          aux = line(1:idx-1)
          line = aux
          lastpass = .true.
       end if
    end do

    if (.not.found) &
       call ferror("next_integer","could not parse xml file",faterr)

  end function next_integer

  !> Read the kpointsandweights entry from the xml.
  subroutine read_kpointsandweights(lu,kpts,w)
    use tools_io
    integer, intent(in) :: lu
    real*8, intent(out) :: kpts(:,:)
    real*8, intent(out) :: w(:)

    logical :: ok, found
    character(len=:), allocatable :: line, key
    integer :: i

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    key = "<kpointsandweights>"
    found = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = trim(adjustl(lower(line)))
       if (index(line,key) > 0) then
          found = .true.
          do i = 1, size(w)
             read (lu,*) kpts(:,i), w(i)
          end do
          exit
       end if
    end do

    if (.not.found) &
       call ferror("read_kpointsandweights","could not parse xml file",faterr)
    
  end subroutine read_kpointsandweights

  !> Read the occupations from the xml.
  subroutine read_occupations(lu,occ)
    use tools_io
    integer, intent(in) :: lu
    real*8, intent(out) :: occ(:,:,:)

    logical :: ok, found
    character(len=:), allocatable :: line, key, aux
    integer :: is, ik

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    key = "<occupations>"
    found = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       aux = lower(line)
       line = trim(adjustl(aux))
       if (index(line,key) > 0) then
          found = .true.

          do is = 1, size(occ,3)
             do ik = 1, size(occ,2)
                ! advance to the k-point with this name
                key = "<k" // string(ik) // ">"
                do while (.true.)
                   ok = getline_raw(lu,line,.true.)
                   if (index(line,key) > 0) exit
                end do
                ! read nstate real numbers
                occ(:,ik,is) = dftb_read_reals1(lu,size(occ,1))
             end do
          end do
          exit

       end if
    end do

    if (.not.found) &
       call ferror("read_occupations","could not parse xml file",faterr)
    
  end subroutine read_occupations

  !> Read a list of n reals from a logical unit.
  function dftb_read_reals1(lu,n) result(x)
    use tools_io
    integer, intent(in) :: lu, n
    real*8 :: x(n)

    integer :: kk, lp
    real*8 :: rdum
    character(len=:), allocatable :: line, aux
    logical :: ok

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    kk = 0
    lp = 1
    ok = getline_raw(lu,aux,.true.)
    line = trim(adjustl(aux))
    do while(.true.)
       if (.not.isreal(rdum,line,lp)) then
          lp = 1
          ok = getline_raw(lu,aux)
          line = trim(adjustl(aux))
          if (.not.ok .or. line(1:1) == "<") exit
       else
          kk = kk + 1
          if (kk > n) call ferror("dftb_read_reals1","exceeded size of the array",2)
          x(kk) = rdum
       endif
    enddo

  endfunction dftb_read_reals1

  !> Read the next atom from the hsd wfc file.
  function next_hsd_atom(lu,at) result(ok)
    use tools_io
    use types
    integer, intent(in) :: lu
    type(dftbatom), intent(out) :: at
    logical :: ok

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    character(len=:), allocatable :: line, word, aux
    integer :: idx, nb, lp, i, n
    real*8 :: rdum

    ok = .false.
    at%norb = 0
    at%nexp = 0
    at%ncoef = 0
    nb = 0
    do while(getline(lu,line,.false.))
       if (len_trim(line) > 0) then
          ! must be an atom?
          idx = index(line,"{")
          at%name = line(1:idx-1)
          at%name = trim(adjustl(at%name))
          nb = nb + 1

          ! read the keywords at this level: atomicnumber and orbital
          do while(getline(lu,line,.true.))
             lp = 1
             word = lgetword(line,lp)
             idx = index(line,"{")
             aux = word
             if (idx > 0) aux = word(1:idx-1)
             word = trim(adjustl(aux))

             if (equal(word,"atomicnumber")) then
                idx = index(line,"=")
                word = line(idx+1:)
                read (word,*) at%z
             elseif (equal(word,"orbital")) then
                nb = nb + 1
                at%norb = at%norb + 1
                ! read the keywords at this level: 
                ! angularmomentum, occupation, cutoff, exponents, coefficients
                do while(getline(lu,line,.true.))
                   lp = 1
                   word = lgetword(line,lp)
                   idx = index(line,"{")
                   aux = word
                   if (idx > 0) aux = word(1:idx-1)
                   word = trim(adjustl(aux))
                   if (equal(word,"angularmomentum")) then
                      idx = index(line,"=")
                      word = line(idx+1:)
                      read (word,*) at%l(at%norb)
                   elseif (equal(word,"occupation")) then
                      idx = index(line,"=")
                      word = line(idx+1:)
                      read (word,*) at%occ(at%norb)
                   elseif (equal(word,"cutoff")) then
                      idx = index(line,"=")
                      word = line(idx+1:)
                      read (word,*) at%cutoff(at%norb)
                   elseif (equal(word,"exponents")) then
                      nb = nb + 1
                      ok = getline(lu,line,.true.)
                      at%nexp(at%norb) = 0
                      lp = 1
                      do while (isreal(rdum,line,lp))
                         at%nexp(at%norb) = at%nexp(at%norb) + 1
                         at%eexp(at%nexp(at%norb),at%norb) = rdum
                      end do
                   elseif (equal(word,"coefficients")) then
                      nb = nb + 1
                      if (at%nexp(at%norb) == 0) &
                         call ferror('next_hsd_atom','coefficients must come after exponents',faterr)
                      do i = 1, at%nexp(at%norb)
                         ok = getline(lu,line,.true.)
                         n = 0
                         lp = 1
                         do while (isreal(rdum,line,lp))
                            n = n + 1
                            at%coef(n,i,at%norb) = rdum
                         end do
                         at%ncoef(i,at%norb) = n
                      end do
                   elseif (equal(word,"}")) then
                      nb = nb - 1
                      if (nb == 1) exit
                   else
                      call ferror("next_hsd_atom","unknown keyword (2): " // string(word),faterr)
                   end if
                end do
             elseif (equal(word,"}")) then
                nb = nb - 1
                if (nb == 0) exit
                exit
             else
                call ferror("next_hsd_atom","unknown keyword (1): " // string(word),faterr)
             end if
          end do

          ! succesfully read an atom
          ok = .true.
          exit
       end if
    end do

  end function next_hsd_atom

  !> Build the interpolation grids for the radial parts of the orbitals.
  subroutine build_interpolation_grid1(ff)
    use types
    use tools_io
    type(field), intent(inout) :: ff

    integer :: it, iorb, istat, ipt

    integer, parameter :: npt = 2001

    do it = 1, size(ff%bas)
       if (allocated(ff%bas(it)%orb)) deallocate(ff%bas(it)%orb)
       allocate(ff%bas(it)%orb(ff%bas(it)%norb),stat=istat)
       if (istat /= 0) call ferror('build_interpolation_grid1',&
          'could not allocate memory for orbitals',faterr)
       do iorb = 1, ff%bas(it)%norb
          ff%bas(it)%orb(iorb)%init = .true.
          ff%bas(it)%orb(iorb)%a = mindist
          ff%bas(it)%orb(iorb)%rmax = ff%bas(it)%cutoff(iorb)
          ff%bas(it)%orb(iorb)%rmax2 = ff%bas(it)%orb(iorb)%rmax * ff%bas(it)%orb(iorb)%rmax
          ff%bas(it)%orb(iorb)%ngrid = npt
          ff%bas(it)%orb(iorb)%b = log(ff%bas(it)%cutoff(iorb) / mindist) / real(npt-1,8)
          ff%bas(it)%orb(iorb)%z = ff%bas(it)%z
          ff%bas(it)%orb(iorb)%qat = 0

          allocate(ff%bas(it)%orb(iorb)%r(npt),ff%bas(it)%orb(iorb)%f(npt),ff%bas(it)%orb(iorb)%fp(npt),ff%bas(it)%orb(iorb)%fpp(npt),stat=istat)
          if (istat /= 0) call ferror('build_interpolation_grid1',&
             'could not allocate memory for orbital radial grids',faterr)
          do ipt = 1, npt
             ff%bas(it)%orb(iorb)%r(ipt) = ff%bas(it)%orb(iorb)%a * exp(ff%bas(it)%orb(iorb)%b * real(ipt-1,8))
             call calculate_rl(ff,it,iorb,ff%bas(it)%orb(iorb)%r(ipt),ff%bas(it)%orb(iorb)%f(ipt),&
                ff%bas(it)%orb(iorb)%fp(ipt),ff%bas(it)%orb(iorb)%fpp(ipt))
          end do
       end do
    end do
  end subroutine build_interpolation_grid1

  !> Calculate the radial part of an orbital and its first and second
  !> derivatives exactly.
  subroutine calculate_rl(ff,it,iorb,r0,f,fp,fpp)
    use types
    type(field), intent(in) :: ff
    integer, intent(in) :: it, iorb
    real*8, intent(in) :: r0
    real*8, intent(out) :: f, fp, fpp

    real*8 :: sumf, sumfp, sumfpp, rx(-1:5), ee, r
    integer :: i, j, l

    ! initialize
    l = ff%bas(it)%l(iorb)
    r = max(r0,mindist)

    ! powers of the distance
    rx = 0d0
    rx(1) = r**l
    do i = 2, maxval(ff%bas(it)%ncoef(:,iorb))
       rx(i) = rx(i-1) * r
    end do
    if (l > 0) rx(0) = rx(1) / r
    if (l > 1) rx(-1) = rx(0) / r

    ! calculate the Rl and its derivatives
    f = 0d0
    fp = 0d0
    fpp = 0d0
    do i = 1, ff%bas(it)%nexp(iorb)
       ee = exp(-ff%bas(it)%eexp(i,iorb) * r)

       sumf = 0d0
       sumfp = 0d0
       sumfpp = 0d0
       do j = 1, ff%bas(it)%ncoef(i,iorb)
          sumf = sumf + ff%bas(it)%coef(j,i,iorb) * rx(j)
          sumfp = sumfp + ff%bas(it)%coef(j,i,iorb) * rx(j-1) * real(l+j-1,8)
          sumfpp = sumfpp + ff%bas(it)%coef(j,i,iorb) * rx(j-2) * real(l+j-1,8) * real(l+j-2,8)
       end do
       f = f + sumf * ee
       fp = fp - ff%bas(it)%eexp(i,iorb) * ee * sumf + sumfp * ee
       fpp = fpp + ff%bas(it)%eexp(i,iorb)**2 * ee * sumf + sumfpp * ee - 2 * ff%bas(it)%eexp(i,iorb) * ee * sumfp
    end do

  end subroutine calculate_rl

end module dftb_private
