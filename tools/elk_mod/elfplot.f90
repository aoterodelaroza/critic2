
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
!   Outputs the electron localisation function (ELF) for 1D, 2D or 3D plotting
!   or a MT+kpts representation of the ELF (ELF.OUT).
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
integer ik,is,ias
integer nr,nri,ir
integer ig,ifg,i
real(8) r,t1,t2,t3
! allocatable arrays
real(8), allocatable :: gwf2mt(:,:),gwf2ir(:)
real(8), allocatable :: rfmt1(:),rfmt2(:),grfir(:)
real(8), allocatable :: grfmt1(:,:),grfmt2(:,:)
real(8), allocatable :: elfmt(:,:),elfir(:)
complex(8), allocatable :: evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: zfft1(:),zfft2(:)
! initialise universal variables
call init0
call init1
! allocate local arrays
allocate(gwf2mt(npmtmax,natmtot),gwf2ir(ngtot))
allocate(rfmt1(npmtmax),rfmt2(npmtmax),grfir(ngtot))
allocate(grfmt1(npmtmax,3),grfmt2(npmtmax,3))
allocate(elfmt(npmtmax,natmtot),elfir(ngtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(zfft1(ngtot),zfft2(ngtot))
! read density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! set the gradient squared to zero
gwf2mt(:,:)=0.d0
gwf2ir(:)=0.d0
do ik=1,nkpt
! get the eigenvectors and occupancies from file
  call getoccsv(filext,ik,vkl(:,ik),occsv(:,ik))
  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! add the valence wavefunction gradient squared
  call gradwf2(ik,evecfv,evecsv,gwf2mt,gwf2ir)
end do
! add core wavefunction gradient squared
call gradwfcr2(gwf2mt)
!------------------------!
!     muffin-tin ELF     !
!------------------------!
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
! convert rho from spherical harmonics to spherical coordinates
  call rbsht(nr,nri,rhomt(:,ias),rfmt1)
! compute the gradient of the density
  call gradrfmt(nr,nri,rlmt(:,1,is),rlmt(:,-1,is),rhomt(:,ias),npmtmax,grfmt1)
! convert gradient to spherical coordinates
  do i=1,3
    call rbsht(nr,nri,grfmt1(:,i),grfmt2(:,i))
  end do
  do i=1,npmt(is)
    r=abs(rfmt1(i))
! square of gradient of rho
    t1=grfmt2(i,1)**2+grfmt2(i,2)**2+grfmt2(i,3)**2
! D for inhomogeneous density
    t2=(1.d0/2.d0)*(gwf2mt(i,ias)-(1.d0/4.d0)*t1/r)
! D0 for uniform electron gas
    t3=(3.d0/5.d0)*((6.d0*pi**2)**(2.d0/3.d0))*(r/2.d0)**(5.d0/3.d0)
! ELF function
    rfmt2(i)=1.d0/(1.d0+(t2/t3)**2)
  end do
! convert ELF from spherical coordinates to spherical harmonics
  call rfsht(nr,nri,rfmt2,elfmt(:,ias))
end do
!--------------------------!
!     interstitial ELF     !
!--------------------------!
! Fourier transform density to G-space
zfft1(:)=rhoir(:)
call zfftifc(3,ngridg,-1,zfft1)
grfir(:)=0.d0
do i=1,3
  zfft2(:)=0.d0
  do ig=1,ngvec
    ifg=igfft(ig)
! take the gradient
    zfft2(ifg)=zi*vgc(i,ig)*zfft1(ifg)
  end do
! Fourier transform gradient to real-space
  call zfftifc(3,ngridg,1,zfft2)
  do ir=1,ngtot
    grfir(ir)=grfir(ir)+dble(zfft2(ir))**2
  end do
end do
do ir=1,ngtot
  r=abs(rhoir(ir))
! D for inhomogeneous density
  t1=(1.d0/2.d0)*(gwf2ir(ir)-(1.d0/4.d0)*grfir(ir)/r)
! D0 for homogeneous electron gas
  t2=(3.d0/5.d0)*((6.d0*pi**2)**(2.d0/3.d0))*(r/2.d0)**(5.d0/3.d0)
! ELF function
  elfir(ir)=1.d0/(1.d0+(t1/t2)**2)
end do
! symmetrise the ELF
call symrf(nrmt,nrmti,npmt,npmtmax,elfmt,elfir)
! plot the ELF to file
select case(task)
case(51)
  open(50,file='ELF1D.OUT',form='FORMATTED')
  open(51,file='ELFLINES.OUT',form='FORMATTED')
  call plot1d(50,51,1,elfmt,elfir)
  close(50)
  close(51)
  write(*,*)
  write(*,'("Info(elfplot):")')
  write(*,'(" 1D ELF plot written to ELF1D.OUT")')
  write(*,'(" vertex location lines written to ELFLINES.OUT")')
case(52)
  open(50,file='ELF2D.OUT',form='FORMATTED')
  call plot2d(.false.,50,1,elfmt,elfir)
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): 2D ELF plot written to ELF2D.OUT")')
case(53)
  open(50,file='ELF3D.OUT',form='FORMATTED')
  call plot3d(50,1,elfmt,elfir)
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): 3D ELF plot written to ELF3D.OUT")')
case(54)
  open(50,file='ELF.OUT',action='WRITE',form='UNFORMATTED')
  write(50) lmmaxi, lmmaxo, nrmtmax, npmtmax, natmtot, ngtot, maxspecies
  write(50) nrmti, nrmt
  write(50) elfmt, elfir
  close(50)
  write(*,*)
  write(*,'("Info(elfplot): ELF mt+kpts description written to ELF.OUT")')
end select
deallocate(gwf2mt,gwf2ir,rfmt1,rfmt2,grfir)
deallocate(grfmt1,grfmt2,elfmt,elfir)
deallocate(evecfv,evecsv,zfft1,zfft2)
return
end subroutine
!EOC

