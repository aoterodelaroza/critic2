    ! xxxx
    integer :: lu
    integer :: nk, ikk, jkk, irr, jrr
    complex*16, allocatable, dimension(:,:,:) :: unk1, unk2, unkp, phase
    complex*16 :: ph
    character(len=:), allocatable :: oname
          nk = f(fid)%nwan(1)*f(fid)%nwan(2)*f(fid)%nwan(3)
          allocate(unk1(f(fid)%n(1),f(fid)%n(2),f(fid)%n(3)))
          allocate(unk2(f(fid)%n(1),f(fid)%n(2),f(fid)%n(3)))
          allocate(unkp(f(fid)%n(1),f(fid)%n(2),f(fid)%n(3)))
          allocate(phase(f(fid)%n(1),f(fid)%n(2),f(fid)%n(3)))
          ! run over k 
          do ikk = 1, nk
             ! run over k'
             do jkk = 1, nk
                ! calculate the exp(i * (k'-k) * r) phase
                do k = 1, f(fid)%n(3)
                   do j = 1, f(fid)%n(2)
                      do i = 1, f(fid)%n(1)
                         phase(i,j,k) = exp(tpi*img*(&
                            (f(fid)%wan_kpt(1,jkk)-f(fid)%wan_kpt(1,ikk))*real(i-1,8)/real(f(fid)%n(1),8)+&
                            (f(fid)%wan_kpt(2,jkk)-f(fid)%wan_kpt(2,ikk))*real(j-1,8)/real(f(fid)%n(2),8)+&
                            (f(fid)%wan_kpt(3,jkk)-f(fid)%wan_kpt(3,ikk))*real(k-1,8)/real(f(fid)%n(3),8)))
                      end do
                   end do
                end do

                ! run over spins
                do is = 1, nspin

                   ! run over bands
                   do ibnd1 = 1, nbnd
                      oname = "WNK." // string(ikk) // "." // string(ibnd1) // "." // string(is)
                      lu = fopen_read(oname,form="unformatted")
                      read(lu) unk1
                      call fclose(lu)
                      do ibnd2 = 1, nbnd
                         oname = "WNK." // string(jkk) // "." // string(ibnd2) // "." // string(is)
                         lu = fopen_read(oname,form="unformatted")
                         read(lu) unk2
                         call fclose(lu)

                         ! conjg(unk) * un'k' * exp(i(k'-k)*r)
                         unkp = conjg(unk1)*unk2*phase

                         ! run over attractors
                         do i = 1, natt1
                            padd = sum(unkp,idg1==i)

                            ! R'
                            do ia = 0, f(fid)%nwan(1)-1
                               do ja = 0, f(fid)%nwan(2)-1
                                  do ka = 0, f(fid)%nwan(3)-1
                                     ! R
                                     do ib = 0, f(fid)%nwan(1)-1
                                        do jb = 0, f(fid)%nwan(2)-1
                                           do kb = 0, f(fid)%nwan(3)-1
                                              ! exp (2*pi*i * (k' * (R' + I) - k * (R + I)))
                                              ph = exp(tpi*img*(&
                                                 f(fid)%wan_kpt(1,jkk) * (ia + ilvec(1,i)) + f(fid)%wan_kpt(2,jkk) * (ja + ilvec(2,i)) + &
                                                 f(fid)%wan_kpt(3,jkk) * (ka + ilvec(3,i)) - f(fid)%wan_kpt(1,ikk) * (ib + ilvec(1,i)) - &
                                                 f(fid)%wan_kpt(2,ikk) * (jb + ilvec(2,i)) - f(fid)%wan_kpt(3,ikk) * (kb + ilvec(3,i))))
                                              ! pack and add
                                              call packidx(ia+ilvec(1,i),ja+ilvec(2,i),ka+ilvec(3,i),ibnd1,imo1,nmo,nbnd,nwan)
                                              call packidx(ib+ilvec(1,i),jb+ilvec(2,i),kb+ilvec(3,i),ibnd2,jmo1,nmo,nbnd,nwan)
                                              sij(imo1,jmo1,iatt(i),is,ndeloc) = sij(imo1,jmo1,iatt(i),is,ndeloc) + padd * ph
                                           end do
                                        end do
                                     end do
                                  end do
                               end do
                            end do

                         end do
                      end do
                   end do
                end do
             end do
          end do
          deallocate(unk1,unk2,unkp,phase)

          ! wannier function normalization factors
          sij = sij / real(nk * nk,8)

          ! write (*,*) "bleh!"
          ! stop 1
          
