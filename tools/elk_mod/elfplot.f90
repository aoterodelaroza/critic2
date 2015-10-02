
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: elfplot
! !INTERFACE:
subroutine elfplot
! !USES:
use modmain
! !DESCRIPTION:
!   Outputs the electron localisation function (ELF) for 1D, 2D or 3D plotting.
!   The spin-averaged ELF is given by
!   $$ f_{\rm ELF}({\bf r})=\frac{1}{1+[D({\bf r})/D^0({\bf r})]^2}, $$
!   where
!   $$ D({\bf r})=\frac{1}{2}\left(\tau({\bf r})-\frac{1}{4}
!    \frac{[\nabla n({\bf r})]^2}{n({\bf r})}\right) $$
!   and
!   $$ \tau({\bf r})=\sum_{i=1}^N \left|\nabla\Psi_i({\bf r})
!    \right|^2 $$
!   is the spin-averaged kinetic energy density from the spinor wavefunctions.
!   The function $D^0$ is the kinetic energy density for the homogeneous
!   electron gas evaluated for $n({\bf r})$:
!   $$ D^0({\bf r})=\frac{3}{5}(6\pi^2)^{2/3}\left(\frac{n({\bf r})}{2}
!    \right)^{5/3}. $$
!   The ELF is useful for the topological classification of bonding. See for
!   example T. Burnus, M. A. L. Marques and E. K. U. Gross [Phys. Rev. A 71,
!   10501 (2005)].
!
! !REVISION HISTORY:
!   Created September 2003 (JKD)
!   Fixed bug found by F. Wagner (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,is,ia,ias
integer ir,i,itp,ig,ifg
real(8) t1,t2
! allocatable arrays
real(8), allocatable :: rftp1(:),rftp2(:),rftp3(:)
real(8), allocatable :: grfmt(:,:,:)
real(8), allocatable :: grfir(:)
real(8), allocatable :: gwf2mt(:,:,:)
real(8), allocatable :: gwf2ir(:)
real(8), allocatable :: elfmt(:,:,:)
real(8), allocatable :: elfir(:)
complex(8), allocatable :: zfft1(:),zfft2(:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
! allocate local arrays
allocate(rftp1(lmmaxvr),rftp2(lmmaxvr),rftp3(lmmaxvr))
allocate(grfmt(lmmaxvr,nrmtmax,3))
allocate(grfir(ngrtot))
allocate(gwf2mt(lmmaxvr,nrmtmax,natmtot))
allocate(gwf2ir(ngrtot))
allocate(elfmt(lmmaxvr,nrmtmax,natmtot))
allocate(elfir(ngrtot))
allocate(zfft1(ngrtot),zfft2(ngrtot))
! allocate first-variational eigenvector array
allocate(evecfv(nmatmax,nstfv))
! allocate second-variational eigenvector array
allocate(evecsv(nstsv,nstsv))
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! generate the core wavefunctions and densities
call gencore
! set the gradient squared to zero
gwf2mt(:,:,:)=0.d0
gwf2ir(:)=0.d0
do ik=1,nkpt
! get the eigenvectors and occupancies from file
  call getoccsv(vkl(:,ik),occsv(:,ik))
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! add the valence wavefunction gradient squared
  call gradwf2(ik,evecfv,evecsv,gwf2mt,gwf2ir)
end do
! add core wavefunction gradient squared
call gradwfcr2(gwf2mt)
!------------------------!
!     muffin-tin ELF     !
!------------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the gradient of the density
    call gradrfmt(lmaxvr,nrmt(is),spr(:,is),lmmaxvr,nrmtmax,rhomt(:,:,ias), &
     grfmt)
    do ir=1,nrmt(is)
! convert rho from spherical harmonics to spherical coordinates
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,rhomt(:,ir,ias),1, &
       0.d0,rftp1,1)
      rftp2(:)=0.d0
! compute the square of the gradient of rho
      do i=1,3
        call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtvr,lmmaxvr,grfmt(:,ir,i),1, &
         0.d0,rftp3,1)
        do itp=1,lmmaxvr
          rftp2(itp)=rftp2(itp)+rftp3(itp)**2
        end do
      end do
      do itp=1,lmmaxvr
! D for inhomogeneous density
        t1=(1.d0/2.d0)*(gwf2mt(itp,ir,ias)-(1.d0/4.d0)*rftp2(itp)/rftp1(itp))
! D0 for uniform electron gas
        t2=(3.d0/5.d0)*((6.d0*pi**2)**(2.d0/3.d0)) &
         *(abs(rftp1(itp))/2.d0)**(5.d0/3.d0)
! ELF function
        rftp3(itp)=1.d0/(1.d0+(t1/t2)**2)
      end do
! convert from spherical coordinates to spherical harmonics
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rfshtvr,lmmaxvr,rftp3,1,0.d0, &
       elfmt(:,ir,ias),1)
    end do
  end do
end do
!--------------------------!
!     interstitial ELF     !
!--------------------------!
! Fourier transform density to G-space
zfft1(:)=rhoir(:)
call zfftifc(3,ngrid,-1,zfft1)
grfir(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
! take the gradient
    zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
  end do
! Fourier transform gradient to real-space
  call zfftifc(3,ngrid,1,zfft2)
  do ir=1,ngrtot
    grfir(ir)=grfir(ir)+dble(zfft2(ir))**2
  end do
end do
do ir=1,ngrtot
! D for inhomogeneous density
  t1=(1.d0/2.d0)*(gwf2ir(ir)-(1.d0/4.d0)*grfir(ir)/rhoir(ir))
! D0 for homogeneous electron gas
  t2=(3.d0/5.d0)*((6.d0*pi**2)**(2.d0/3.d0))*(abs(rhoir(ir))/2.d0)**(5.d0/3.d0)
! ELF function
  elfir(ir)=1.d0/(1.d0+(t1/t2)**2)
end do
! symmetrise the ELF
call symrf(1,elfmt,elfir)
! plot the ELF to file
select case(task)
case(51)
  open(50,file='ELF1D.OUT',action='WRITE',form='FORMATTED')
  open(51,file='ELFLINES.OUT',action='WRITE',form='FORMATTED')
  call plot1d(50,51,1,lmaxvr,lmmaxvr,elfmt,elfir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(elfplot):")')
  write(*,'(" 1D ELF plot written to ELF1D.OUT")')
  write(*,'(" vertex location lines written to ELFLINES.OUT")')
case(52)
  open(50,file='ELF2D.OUT',action='WRITE',form='FORMATTED')
  call plot2d(50,1,lmaxvr,lmmaxvr,elfmt,elfir)
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): 2D ELF plot written to ELF2D.OUT")')
case(53)
  open(50,file='ELF3D.OUT',action='WRITE',form='FORMATTED')
  call plot3d(50,1,lmaxvr,lmmaxvr,elfmt,elfir)
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): 3D ELF plot written to ELF3D.OUT")')
case(54)
  open(50,file='ELF.OUT',action='WRITE',form='UNFORMATTED')
  write(50) lmmaxvr, nrmtmax, natmtot, ngrtot
  write(50) elfmt, elfir
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): ELF mt+kpts description written to ELF.OUT")')
end select
write(*,*)
deallocate(rftp1,rftp2,rftp3)
deallocate(grfmt,grfir,gwf2mt,gwf2ir,elfmt,elfir)
deallocate(zfft1,zfft2,evecfv,evecsv)
return
end subroutine
!EOC

