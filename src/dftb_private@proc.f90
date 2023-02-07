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

! DFTB+ wavefunctions and readers
submodule (dftb_private) proc
  implicit none

  !xx! private procedures
  ! function next_logical(lu,line0,key0) result(next)
  ! function next_integer(lu,line0,key0) result(next)
  ! subroutine read_kpointsandweights(lu,kpts,w)
  ! subroutine read_occupations(lu,occ)
  ! function dftb_read_reals1(lu,n) result(x)
  ! function next_hsd_atom(lu,at) result(ok)
  ! subroutine build_interpolation_grid1(ff)
  ! subroutine calculate_rl(ff,it,iorb,r0,f,fp,fpp)
  ! subroutine realloc_dftbbasis(a,nnew)

  ! minimum distance (bohr)
  real*8, parameter :: mindist = 1d-6

contains

  !> Deallocate data
  module subroutine dftb_end(f)
    class(dftbwfn), intent(inout) :: f

    if (allocated(f%docc)) deallocate(f%docc)
    if (allocated(f%dkpt)) deallocate(f%dkpt)
    if (allocated(f%evecr)) deallocate(f%evecr)
    if (allocated(f%evecc)) deallocate(f%evecc)
    if (allocated(f%ispec)) deallocate(f%ispec)
    if (allocated(f%idxorb)) deallocate(f%idxorb)
    if (allocated(f%bas)) deallocate(f%bas)
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    if (f%isealloc) then
       if (associated(f%e)) deallocate(f%e)
    end if
    nullify(f%e)
    f%isealloc = .false.

  end subroutine dftb_end

  !> Read the information for a DFTB+ field from the detailed.xml,
  !> eigenvec.bin, and the basis set definition in HSD format.
  module subroutine dftb_read(f,filexml,filebin,filehsd,env,errmsg,ti)
    use types, only: anyatom, species
    use tools_io, only: fopen_read, getline_raw, lower, string, fclose
    use param, only: tpi
    class(dftbwfn), intent(inout) :: f !< Output field
    character*(*), intent(in) :: filexml !< The detailed.xml file
    character*(*), intent(in) :: filebin !< The eigenvec.bin file
    character*(*), intent(in) :: filehsd !< The definition of the basis set in hsd format
    type(environ), intent(in), target :: env !< environment of the cell
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: lu, i, j, k, idum, n, id
    character(len=:), allocatable :: line
    logical :: ok, iread(5), iserr
    type(dftbbasis) :: at
    real*8, allocatable :: dw(:)

    errmsg = "Error reading file"

    ! detailed.xml, first pass
    lu = fopen_read(filexml,ti=ti)
    if (lu < 0) goto 999
    iread = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = adjustl(lower(line))
       if (index(line,"<real>") > 0) then
          f%isreal = next_logical(lu,line,"real",iserr)
          iread(1) = .true.
       elseif (index(line,"<nrofkpoints>") > 0) then
          f%nkpt = next_integer(lu,line,"nrofkpoints",iserr)
          iread(2) = .true.
       elseif (index(line,"<nrofspins>") > 0) then
          f%nspin = next_integer(lu,line,"nrofspins",iserr)
          iread(3) = .true.
       elseif (index(line,"<nrofstates>") > 0) then
          f%nstates = next_integer(lu,line,"nrofstates",iserr)
          iread(4) = .true.
       elseif (index(line,"<nroforbitals>") > 0) then
          f%norb = next_integer(lu,line,"nroforbitals",iserr)
          iread(5) = .true.
       end if
    end do
    if (iserr) then
       errmsg = "could not parse xml"
       goto 999
    end if
    if (.not.all(iread)) then
       errmsg = "missing information in xml"
       goto 999
    end if

    ! detailed.xml, second pass
    if (allocated(f%dkpt)) deallocate(f%dkpt)
    if (allocated(f%docc)) deallocate(f%docc)
    allocate(f%dkpt(3,f%nkpt),dw(f%nkpt),f%docc(f%nstates,f%nkpt,f%nspin))
    rewind(lu)
    call read_kpointsandweights(lu,f%dkpt,dw,iserr)
    if (iserr) then
       errmsg ="error reading kpoints and weights"
       goto 999
    end if
    call read_occupations(lu,f%docc,iserr)
    if (iserr) then
       errmsg ="error reading occupations"
       goto 999
    end if
    call fclose(lu)

    ! use occupations times weights; scale the kpts by 2pi
    do i = 1, f%nkpt
       f%docc(:,i,:) = f%docc(:,i,:) * dw(i)
    end do
    f%dkpt = f%dkpt * tpi
    deallocate(dw)

    ! eigenvec.bin
    lu = fopen_read(filebin,"unformatted",ti=ti)
    if (lu < 0) goto 999
    read (lu,err=999,end=999) idum ! this is the identity number

    ! read the eigenvectors
    if (f%isreal) then
       if (allocated(f%evecr)) deallocate(f%evecr)
       allocate(f%evecr(f%norb,f%nstates,f%nspin))
       do i = 1, f%nspin
          do k = 1, f%nstates
             read (lu,err=999,end=999) f%evecr(:,k,i)
          end do
       end do
    else
       if (allocated(f%evecc)) deallocate(f%evecc)
       allocate(f%evecc(f%norb,f%nstates,f%nkpt,f%nspin))
       do i = 1, f%nspin
          do j = 1, f%nkpt
             do k = 1, f%nstates
                read (lu,err=999,end=999) f%evecc(:,k,j,i)
             end do
          end do
       end do
    end if
    call fclose(lu)

    ! open the hsd file with the basis definition
    lu = fopen_read(filehsd,ti=ti)
    if (lu < 0) goto 999
    if (allocated(f%bas)) deallocate(f%bas)
    allocate(f%bas(10))
    n = 0
    do while(next_hsd_atom(lu,at,iserr))
       if (.not.any(env%spc(1:env%nspc)%z == at%z)) cycle
       n = n + 1
       if (n > size(f%bas)) call realloc_dftbbasis(f%bas,2*n)
       f%bas(n) = at
    end do
    if (iserr) then
       errmsg = "error reading hsd atoms"
       goto 999
    end if
    call realloc_dftbbasis(f%bas,n)
    call fclose(lu)

    ! tie the atomic numbers to the basis types
    if (allocated(f%ispec)) deallocate(f%ispec)
    allocate(f%ispec(env%nspc))
    f%ispec = 0
    do i = 1, env%nspc
       do j = 1, n
          if (f%bas(j)%z == env%spc(i)%z) then
             f%ispec(i) = j
             exit
          end if
       end do
    end do

    ! indices for the atomic orbitals, for array sizes later on
    if (allocated(f%idxorb)) deallocate(f%idxorb)
    allocate(f%idxorb(env%ncell))
    n = 0
    f%maxnorb = 0
    f%maxlm = 0
    do i = 1, env%ncell
       id = f%ispec(env%at(i)%is)
       if (id == 0) then
          lu = -1
          errmsg = "basis missing for atomic number"
          goto 999
       end if

       f%maxnorb = max(f%maxnorb,f%bas(id)%norb)
       f%idxorb(i) = n + 1
       do j = 1, f%bas(id)%norb
          n = n + 2*f%bas(id)%l(j) + 1
          f%maxlm = max(f%maxlm,f%bas(id)%l(j))
       end do
    end do
    f%midxorb = n
    f%maxlm = (f%maxlm+3)*(f%maxlm+3)

    call build_interpolation_grid1(f,errmsg)
    if (len_trim(errmsg) > 0) then
       lu = -1
       goto 999
    end if

    ! find the individual species cutoffs and maximum cutoff
    if (allocated(f%spcutoff)) deallocate(f%spcutoff)
    allocate(f%spcutoff(env%nspc,2))
    f%spcutoff = 0d0
    f%globalcutoff = -1d0
    do i = 1, env%nspc
       id = f%ispec(i)
       do j = 1, f%bas(id)%norb
          f%spcutoff(i,2) = max(f%spcutoff(i,2),f%bas(id)%cutoff(j))
       end do
       f%globalcutoff = max(f%globalcutoff,f%spcutoff(i,2))
    end do

    if (f%isealloc) then
       if (associated(f%e)) deallocate(f%e)
    end if
    nullify(f%e)
    if (f%globalcutoff >= env%dmax0 .and..not.env%ismolecule) then
       ! Create a new environment to satisfy all searches.
       ! The environment contains all the atoms in molecules anyway.
       f%isealloc = .true.
       nullify(f%e)
       allocate(f%e)
       call f%e%extend(env,f%globalcutoff)
    else
       ! keep a pointer to the environment
       f%isealloc = .false.
       f%e => env
    end if

    errmsg = ""
    return
999 continue
    if (lu > 0) call fclose(lu)

  end subroutine dftb_read

  !> Calculate the density and derivatives of a DFTB+ field (f) up to
  !> the nder degree (max = 2). xpos is in Cartesian coordinates. In
  !> output, the density (rho), the gradient (grad), the Hessian (h),
  !> and the G(r) kinetic energy density (gkin).  This routine is
  !> thread-safe.
  module subroutine rho2(f,xpos,exact,nder,rho,grad,h,gkin)
    use tools_math, only: tosphere, genylm, ylmderiv
    use param, only: img, icrd_cart
    class(dftbwfn), intent(inout) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    logical, intent(in) :: exact !< exact or approximate calculation
    integer, intent(in) :: nder  !< Number of derivatives
    real*8, intent(out) :: rho !< Density
    real*8, intent(out) :: grad(3) !< Gradient
    real*8, intent(out) :: h(3,3) !< Hessian
    real*8, intent(out) :: gkin !< G(r), kinetic energy density

    integer :: ion, it, is, istate, ik, iorb, i, l, m, lmax
    integer :: ixorb, ixorb0
    real*8 :: xion(3), rcut, r, tp(2)
    complex*16, allocatable :: xao(:), xaol(:), xaolp(:,:), xaolpp(:,:), xaop(:,:), xaopp(:,:)
    complex*16, allocatable :: phase(:,:), ylm(:)
    complex*16 :: xmo, xmop(3), xmopp(6)
    real*8, allocatable :: rao(:), raol(:), raolp(:,:), raolpp(:,:), raop(:,:), raopp(:,:)
    real*8, allocatable :: rl(:,:), rlp(:,:), rlpp(:,:)
    real*8, allocatable :: phi(:,:,:), phip(:,:,:,:), phipp(:,:,:,:)
    integer, allocatable :: idxion(:)
    real*8 :: rmo, rmop(3), rmopp(6), lvecc(3)
    integer :: nenvl, ionl
    integer :: imin, ip, im, iphas
    complex*16 :: xgrad1(3), xgrad2(3), xhess1(6), xhess2(6)
    integer :: nenv, lvec(3), ierr, lenv(3)
    integer, allocatable :: eid(:)
    real*8, allocatable :: dist(:)

    real*8, parameter :: docc_cutoff = 1d-20

    ! calculate the environment of the input point
    rho = 0d0
    grad = 0d0
    h = 0d0
    gkin = 0d0
    call f%e%list_near_atoms(xpos,icrd_cart,.false.,nenv,ierr,eid,dist,lvec,up2dsp=f%spcutoff)
    lvecc = f%e%x2c(real(lvec,8))
    if (ierr > 0) return ! could happen if in a molecule and very far -> zero

    ! precalculate the quantities that depend only on the environment
    allocate(rl(f%maxnorb,nenv),rlp(f%maxnorb,nenv),rlpp(f%maxnorb,nenv),ylm(f%maxlm))
    allocate(idxion(nenv))
    if (.not.f%isreal) then
       allocate(phase(nenv,f%nkpt))
    end if
    allocate(phi(f%maxlm,f%maxnorb,nenv))
    phi = 0d0
    if (nder >= 1) then
       allocate(phip(3,f%maxlm,f%maxnorb,nenv))
       phip = 0d0
    end if
    if (nder >= 2) then
       allocate(phipp(6,f%maxlm,f%maxnorb,nenv))
       phipp = 0d0
    end if

    nenvl = 0
    do ion = 1, nenv
       xion = xpos - (f%e%at(eid(ion))%r + lvecc)
       it = f%ispec(f%e%at(eid(ion))%is)

       ! apply the distance cutoff
       rcut = maxval(f%bas(it)%cutoff(1:f%bas(it)%norb))
       if (dist(ion) > rcut) cycle

       ! write down this atom
       nenvl = nenvl + 1
       idxion(nenvl) = eid(ion)

       ! calculate the spherical harmonics contributions for this atom
       lmax = maxval(f%bas(it)%l(1:f%bas(it)%norb))
       call tosphere(xion,r,tp)
       r = max(r,mindist)
       call genylm(lmax+2,tp,ylm)

       ! calculate the radial contributions for this atom
       do iorb = 1, f%bas(it)%norb
          if (dist(ion) > f%bas(it)%cutoff(iorb)) cycle
          if (exact) then
             call calculate_rl(f,it,iorb,dist(ion),rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl))
          else
             call f%bas(it)%orb(iorb)%interp(dist(ion),rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl))
          end if
       end do

       ! populate the phi array and calculate the derivatives
       do iorb = 1, f%bas(it)%norb
          l = f%bas(it)%l(iorb)
          imin = l*l+1

          ! m = 0
          ip = imin+l
          call ylmderiv(ylm,r,l,0,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl),nder,xgrad1,xhess1)
          phi(ip,iorb,nenvl) = rl(iorb,nenvl) * real(ylm(ip),8)
          if (nder >= 1) &
             phip(:,ip,iorb,nenvl) = real(xgrad1,8)
          if (nder >= 2) &
             phipp(:,ip,iorb,nenvl) = real(xhess1,8)

          ! |m| > 0
          do m = 1, l
             ip = imin + l + m
             im = imin + l - m
             iphas = (-1)**m
             call ylmderiv(ylm,r,l, m,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl),nder,xgrad1,xhess1)
             call ylmderiv(ylm,r,l,-m,rl(iorb,nenvl),rlp(iorb,nenvl),rlpp(iorb,nenvl),nder,xgrad2,xhess2)

             phi(ip,iorb,nenvl) = rl(iorb,nenvl) * real(iphas*ylm(ip)+ylm(im),8) / sqrt(2d0)
             phi(im,iorb,nenvl) = rl(iorb,nenvl) * real(-iphas*img*ylm(ip)+img*ylm(im),8) / sqrt(2d0)
             if (nder < 1) cycle
             phip(:,ip,iorb,nenvl) = real(iphas*xgrad1+xgrad2,8) / sqrt(2d0)
             phip(:,im,iorb,nenvl) = real(-iphas*img*xgrad1+img*xgrad2,8) / sqrt(2d0)
             if (nder < 2) cycle
             phipp(:,ip,iorb,nenvl) = real(iphas*xhess1+xhess2,8) / sqrt(2d0)
             phipp(:,im,iorb,nenvl) = real(-iphas*img*xhess1+img*xhess2,8) / sqrt(2d0)
          end do
       end do

       lenv = floor(f%e%xr2x(f%e%at(eid(ion))%x) + lvec)
       ! calculate the phases
       if (.not.f%isreal) then
          do ik = 1, f%nkpt
             phase(nenvl,ik) = exp(img * dot_product(lenv,f%dkpt(:,ik)))
          end do
       end if
    end do
    deallocate(rl,rlp,rlpp,ylm)

    if (.not.f%isreal) then
       allocate(xao(f%midxorb),xaol(f%midxorb))
       if (nder >= 1) &
          allocate(xaolp(3,f%midxorb),xaop(3,f%midxorb))
       if (nder >= 2) &
          allocate(xaolpp(6,f%midxorb),xaopp(6,f%midxorb))

       ! complex version; for crystals with non-gamma k-points
       ! run over spins
       do is = 1, f%nspin
          ! run over k-points
          do ik = 1, f%nkpt
             ! run over the states
             do istate = 1, f%nstates
                if (abs(f%docc(istate,ik,is)) < docc_cutoff) cycle

                ! determine the atomic contributions. The unit cell atoms
                ! in critic are in the same order as in dftb+, which is
                ! why indexing the evecc/evecr works.
                xao = 0d0
                if (nder >= 1) &
                   xaop = 0d0
                if (nder >= 2) &
                   xaopp = 0d0
                ! run over atoms
                do ionl = 1, nenvl
                   ion = idxion(ionl)
                   it = f%ispec(f%e%at(ion)%is)

                   ! run over atomic orbitals
                   ixorb0 = f%idxorb(f%e%at(ion)%cidx)
                   ixorb = ixorb0 - 1
                   do iorb = 1, f%bas(it)%norb
                      ! run over ms for the same l
                      do i = f%bas(it)%l(iorb)*f%bas(it)%l(iorb)+1, (f%bas(it)%l(iorb)+1)*(f%bas(it)%l(iorb)+1)
                         ixorb = ixorb + 1
                         xaol(ixorb) = phi(i,iorb,ionl)
                         if (nder >= 1) &
                            xaolp(:,ixorb) = phip(:,i,iorb,ionl)
                         if (nder >= 2) &
                            xaolpp(:,ixorb) = phipp(:,i,iorb,ionl)
                      end do
                   end do ! iorb

                   xao(ixorb0:ixorb) = xao(ixorb0:ixorb) + xaol(ixorb0:ixorb) * phase(ionl,ik)
                   if (nder >= 1) &
                      xaop(:,ixorb0:ixorb) = xaop(:,ixorb0:ixorb) + xaolp(:,ixorb0:ixorb) * phase(ionl,ik)
                   if (nder >= 2) &
                      xaopp(:,ixorb0:ixorb) = xaopp(:,ixorb0:ixorb) + xaolpp(:,ixorb0:ixorb) * phase(ionl,ik)
                end do ! ion

                ! calculate the value of this extended orbital and its derivatives
                xmo = 0d0
                xmop = 0d0
                xmopp = 0d0
                do i = 1, f%midxorb
                   xmo = xmo + conjg(xao(i))*f%evecc(i,istate,ik,is)
                   if (nder >= 1) &
                      xmop = xmop + conjg(xaop(:,i))*f%evecc(i,istate,ik,is)
                   if (nder >= 2) &
                      xmopp = xmopp + conjg(xaopp(:,i))*f%evecc(i,istate,ik,is)
                end do

                ! accumulate properties
                rho = rho + real(conjg(xmo)*xmo,8) * f%docc(istate,ik,is)
                if (nder < 1) cycle
                grad = grad + real(conjg(xmop)*xmo+conjg(xmo)*xmop,8) * f%docc(istate,ik,is)
                gkin = gkin + real(conjg(xmop(1))*xmop(1)+conjg(xmop(2))*xmop(2)+conjg(xmop(3))*xmop(3),8) * f%docc(istate,ik,is)
                if (nder < 2) cycle
                h(1,1) = h(1,1) + real(conjg(xmopp(1))*xmo+conjg(xmop(1))*xmop(1)&
                                + conjg(xmop(1))*xmop(1)+conjg(xmo)*xmopp(1),8) * f%docc(istate,ik,is)
                h(1,2) = h(1,2) + real(conjg(xmopp(2))*xmo+conjg(xmop(1))*xmop(2)&
                                + conjg(xmop(2))*xmop(1)+conjg(xmo)*xmopp(2),8) * f%docc(istate,ik,is)
                h(1,3) = h(1,3) + real(conjg(xmopp(3))*xmo+conjg(xmop(1))*xmop(3)&
                                + conjg(xmop(3))*xmop(1)+conjg(xmo)*xmopp(3),8) * f%docc(istate,ik,is)
                h(2,2) = h(2,2) + real(conjg(xmopp(4))*xmo+conjg(xmop(2))*xmop(2)&
                                + conjg(xmop(2))*xmop(2)+conjg(xmo)*xmopp(4),8) * f%docc(istate,ik,is)
                h(2,3) = h(2,3) + real(conjg(xmopp(5))*xmo+conjg(xmop(2))*xmop(3)&
                                + conjg(xmop(3))*xmop(2)+conjg(xmo)*xmopp(5),8) * f%docc(istate,ik,is)
                h(3,3) = h(3,3) + real(conjg(xmopp(6))*xmo+conjg(xmop(3))*xmop(3)&
                                + conjg(xmop(3))*xmop(3)+conjg(xmo)*xmopp(6),8) * f%docc(istate,ik,is)
             end do ! states
          end do ! k-points
       end do ! spins
       deallocate(xao,xaol)
       if (nder >= 1) &
          deallocate(xaolp,xaop)
       if (nder >= 2) &
          deallocate(xaolpp,xaopp)
       deallocate(phase)
    else
       allocate(rao(f%midxorb),raol(f%midxorb))
       if (nder >= 1) &
          allocate(raolp(3,f%midxorb),raop(3,f%midxorb))
       if (nder >= 2) &
          allocate(raolpp(6,f%midxorb),raopp(6,f%midxorb))
       ! real, for molecules or crystals with a single k-point at gamma
       ! run over spins
       do is = 1, f%nspin
          ! run over the states
          do istate = 1, f%nstates

             ! determine the atomic contributions. The unit cell atoms
             ! in critic are in the same order as in dftb+, which is
             ! why indexing the evecc/evecr works.
             rao = 0d0
             if (nder >= 1) &
                raop = 0d0
             if (nder >= 2) &
                raopp = 0d0
             ! run over atoms
             do ionl = 1, nenvl
                ion = idxion(ionl)
                it = f%ispec(f%e%at(ion)%is)

                ! run over atomic orbitals
                ixorb0 = f%idxorb(f%e%at(ion)%cidx)
                ixorb = ixorb0 - 1
                do iorb = 1, f%bas(it)%norb
                   ! run over ms for the same l
                   do i = f%bas(it)%l(iorb)*f%bas(it)%l(iorb)+1, (f%bas(it)%l(iorb)+1)*(f%bas(it)%l(iorb)+1)
                      ixorb = ixorb + 1
                      raol(ixorb) = phi(i,iorb,ionl)
                      if (nder >= 1) &
                         raolp(:,ixorb) = phip(:,i,iorb,ionl)
                      if (nder >= 2) &
                         raolpp(:,ixorb) = phipp(:,i,iorb,ionl)
                   end do
                end do ! iorb

                rao(ixorb0:ixorb) = rao(ixorb0:ixorb) + raol(ixorb0:ixorb)
                if (nder >= 1) &
                   raop(:,ixorb0:ixorb) = raop(:,ixorb0:ixorb) + raolp(:,ixorb0:ixorb)
                if (nder >= 2) &
                   raopp(:,ixorb0:ixorb) = raopp(:,ixorb0:ixorb) + raolpp(:,ixorb0:ixorb)
             end do ! ion

             ! calculate the value of this extended orbital and its derivatives
             rmo = 0d0
             rmop = 0d0
             rmopp = 0d0
             do i = 1, f%midxorb
                rmo = rmo + rao(i)*f%evecr(i,istate,is)
                if (nder >= 1) &
                   rmop = rmop + raop(:,i)*f%evecr(i,istate,is)
                if (nder >= 2) &
                   rmopp = rmopp + raopp(:,i)*f%evecr(i,istate,is)
             end do

             ! accumulate properties
             rho = rho + (rmo*rmo) * f%docc(istate,1,is)
             if (nder < 1) cycle
             grad = grad + (rmop*rmo+rmo*rmop) * f%docc(istate,1,is)
             gkin = gkin + (rmop(1)*rmop(1)+rmop(2)*rmop(2)+rmop(3)*rmop(3)) * f%docc(istate,1,is)
             if (nder < 2) cycle
             h(1,1) = h(1,1) + (rmopp(1)*rmo+rmop(1)*rmop(1)+rmop(1)*rmop(1)+rmo*rmopp(1)) * f%docc(istate,1,is)
             h(1,2) = h(1,2) + (rmopp(2)*rmo+rmop(1)*rmop(2)+rmop(2)*rmop(1)+rmo*rmopp(2)) * f%docc(istate,1,is)
             h(1,3) = h(1,3) + (rmopp(3)*rmo+rmop(1)*rmop(3)+rmop(3)*rmop(1)+rmo*rmopp(3)) * f%docc(istate,1,is)
             h(2,2) = h(2,2) + (rmopp(4)*rmo+rmop(2)*rmop(2)+rmop(2)*rmop(2)+rmo*rmopp(4)) * f%docc(istate,1,is)
             h(2,3) = h(2,3) + (rmopp(5)*rmo+rmop(2)*rmop(3)+rmop(3)*rmop(2)+rmo*rmopp(5)) * f%docc(istate,1,is)
             h(3,3) = h(3,3) + (rmopp(6)*rmo+rmop(3)*rmop(3)+rmop(3)*rmop(3)+rmo*rmopp(6)) * f%docc(istate,1,is)
          end do ! states
       end do ! spins
       deallocate(rao,raol)
       if (nder >= 1) deallocate(raolp,raop)
       if (nder >= 2) deallocate(raolpp,raopp)
    end if
    deallocate(idxion,phi)
    if (allocated(phip)) deallocate(phip)
    if (allocated(phipp)) deallocate(phipp)
    ! clean up
    h(2,1) = h(1,2)
    h(3,1) = h(1,3)
    h(3,2) = h(2,3)
    gkin = 0.50 * gkin

  end subroutine rho2

  !xx! private procedures

  !> Read the next logical value in xml format.
  function next_logical(lu,line0,key0,iserr) result(next)
    use tools_io, only: lower, getline_raw
    integer, intent(in) :: lu
    character*(*), intent(in) :: line0, key0
    logical, intent(out) :: iserr
    logical :: next

    character(len=:), allocatable :: line, key, aux
    integer :: idx, idxn, idxy
    logical :: found, ok, lastpass

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    iserr = .false.
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

    if (.not.found) iserr = .true.

  end function next_logical

  !> Read the next integer value in xml format.
  function next_integer(lu,line0,key0,iserr) result(next)
    use tools_io, only: lower, isinteger, getline_raw
    integer, intent(in) :: lu
    character*(*), intent(in) :: line0, key0
    logical, intent(out) :: iserr
    integer :: next

    character(len=:), allocatable :: line, key, aux
    integer :: idx, lp
    logical :: found, ok, lastpass

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.
    iserr = .false.

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

    if (.not.found) iserr = .true.

  end function next_integer

  !> Read the kpointsandweights entry from the xml.
  subroutine read_kpointsandweights(lu,kpts,w,iserr)
    use tools_io, only: getline_raw, lower
    integer, intent(in) :: lu
    real*8, intent(out) :: kpts(:,:)
    real*8, intent(out) :: w(:)
    logical, intent(out) :: iserr

    logical :: ok, found
    character(len=:), allocatable :: line, key
    integer :: i

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    iserr = .false.
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

    if (.not.found) iserr = .true.

  end subroutine read_kpointsandweights

  !> Read the occupations from the xml.
  subroutine read_occupations(lu,occ,iserr)
    use tools_io, only: getline_raw, lower, string
    integer, intent(in) :: lu
    real*8, intent(out) :: occ(:,:,:)
    logical, intent(out) :: iserr

    logical :: ok, found
    character(len=:), allocatable :: line, key, aux
    integer :: is, ik

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    iserr = .true.
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
                occ(:,ik,is) = dftb_read_reals1(lu,size(occ,1),iserr)
                if (iserr) return
             end do
          end do
          exit

       end if
    end do

    if (.not.found) iserr = .true.

  end subroutine read_occupations

  !> Read a list of n reals from a logical unit.
  function dftb_read_reals1(lu,n,iserr) result(x)
    use tools_io, only: getline_raw, isreal
    integer, intent(in) :: lu, n
    logical, intent(out) :: iserr
    real*8 :: x(n)

    integer :: kk, lp
    real*8 :: rdum
    character(len=:), allocatable :: line, aux
    logical :: ok

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    iserr = .false.
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
          if (kk > n) then
             iserr = .true.
             return
          end if
          x(kk) = rdum
       endif
    enddo

  endfunction dftb_read_reals1

  !> Read the next atom from the hsd wfc file.
  function next_hsd_atom(lu,at,iserr) result(ok)
    use types, only: realloc
    use tools_io, only: lgetword, equal, getline, isreal, string
    integer, intent(in) :: lu
    type(dftbbasis), intent(out) :: at
    logical, intent(out) :: iserr
    logical :: ok

    character(len=:), allocatable :: line, word, aux
    integer :: idx, nb, lp, i, j, k, n
    real*8 :: rdum
    real*8, allocatable :: caux(:)

    ! see xx(note1)xx, tools_io.f90 for the use of line and aux.

    iserr = .false.
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
             idx = index(word,"{")
             if (idx > 0) then
                aux = word(1:idx-1)
                word = trim(adjustl(aux))
             end if

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
                   idx = index(word,"{")
                   if (idx > 0) then
                      aux = word(1:idx-1)
                      word = trim(adjustl(aux))
                   end if
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
                      idx = index(line,"{")
                      aux = line(idx+1:)
                      line = aux

                      at%nexp(at%norb) = 0
                      do while (.true.)
                         lp = 1
                         do while (isreal(rdum,line,lp))
                            at%nexp(at%norb) = at%nexp(at%norb) + 1
                            at%eexp(at%nexp(at%norb),at%norb) = rdum
                         end do
                         if (index(line,"}") > 0) exit
                         ok = getline(lu,line,.true.)
                         if (.not.ok) exit
                      end do
                   elseif (equal(word,"coefficients")) then
                      if (at%nexp(at%norb) == 0) then
                         iserr = .true.
                         return
                      end if
                      idx = index(line,"{")
                      aux = line(idx+1:)
                      line = aux

                      allocate(caux(10))
                      n = 0
                      do while(.true.)
                         lp = 1
                         do while (isreal(rdum,line,lp))
                            n = n + 1
                            if (n > size(caux,1)) call realloc(caux,2*n)
                            caux(n) = rdum
                         end do
                         if (index(line,"}") > 0) exit
                         ok = getline(lu,line,.true.)
                         if (.not.ok) exit
                      end do
                      if (mod(n,at%nexp(at%norb)) /= 0) then
                         iserr = .true.
                         return
                      end if

                      k = 0
                      do i = 1, at%nexp(at%norb)
                         at%ncoef(i,at%norb) = n / at%nexp(at%norb)
                         do j = 1, at%ncoef(i,at%norb)
                            k = k + 1
                            at%coef(j,i,at%norb) = caux(k)
                         end do
                      end do
                      deallocate(caux)
                   elseif (equal(word,"}")) then
                      nb = nb - 1
                      if (nb == 1) exit
                   else
                      iserr = .true.
                      return
                   end if
                end do
             elseif (equal(word,"}")) then
                nb = nb - 1
                if (nb == 0) exit
                exit
             else
                iserr = .true.
                return
             end if
          end do

          ! successfully read an atom
          ok = .true.
          exit
       end if
    end do

  end function next_hsd_atom

  !> Build the interpolation grids for the radial parts of the orbitals.
  subroutine build_interpolation_grid1(ff,errmsg)
    class(dftbwfn), intent(inout) :: ff
    character(len=:), allocatable, intent(out) :: errmsg

    integer :: it, iorb, istat, ipt

    integer, parameter :: npt = 2001

    errmsg = ""
    do it = 1, size(ff%bas)
       if (allocated(ff%bas(it)%orb)) deallocate(ff%bas(it)%orb)
       allocate(ff%bas(it)%orb(ff%bas(it)%norb),stat=istat)
       if (istat /= 0) then
          errmsg = "could not allocate memory for orbitals"
          return
       end if
       do iorb = 1, ff%bas(it)%norb
          ff%bas(it)%orb(iorb)%isinit = .true.
          ff%bas(it)%orb(iorb)%a = mindist
          ff%bas(it)%orb(iorb)%rmax = ff%bas(it)%cutoff(iorb)
          ff%bas(it)%orb(iorb)%rmax2 = ff%bas(it)%orb(iorb)%rmax * ff%bas(it)%orb(iorb)%rmax
          ff%bas(it)%orb(iorb)%ngrid = npt
          ff%bas(it)%orb(iorb)%b = log(ff%bas(it)%cutoff(iorb) / mindist) / real(npt-1,8)
          ff%bas(it)%orb(iorb)%z = ff%bas(it)%z
          ff%bas(it)%orb(iorb)%qat = 0

          allocate(ff%bas(it)%orb(iorb)%r(npt),ff%bas(it)%orb(iorb)%f(npt),&
                   ff%bas(it)%orb(iorb)%fp(npt),ff%bas(it)%orb(iorb)%fpp(npt),stat=istat)
          if (istat /= 0) then
             errmsg = "could not allocate memory for orbital radial grids"
             return
          end if
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
    class(dftbwfn), intent(in) :: ff
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

  !> Adapt the size of an allocatable 1D type(atom) array
  subroutine realloc_dftbbasis(a,nnew)
    type(dftbbasis), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(dftbbasis), allocatable :: temp(:)
    integer :: nold

    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))

    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc_dftbbasis

end submodule proc
