       ! ! check -- are the wannier functions orthonormal?
       ! if (allocated(f(fid)%fwan)) then
       !    do is = 1, nspin
       !       write (*,*) "spin = ", is
       !       do imo = 1, nmo
       !          call unpackidx(imo,ia,ja,ka,iba)
       !          write (uout,'(99(A))') "+ Overlaps of: b=", string(iba), " Rx=", string(ia), &
       !             " Ry=", string(ja), " Rz=", string(ka)
       !          do jmo = 1, nmo
       !             call unpackidx(jmo,ib,jb,kb,ibb)
       !             asum = 0d0
       !             do m1 = 0, f(fid)%nwan(1)-1
       !                do m2 = 0, f(fid)%nwan(2)-1
       !                   do m3 = 0, f(fid)%nwan(3)-1
       !                      idx1 = (/m1,m2,m3/) + (/ia,ja,ka/)
       !                      idx1 = modulo(idx1,f(fid)%nwan)
       !                      idx2 = (/m1,m2,m3/) + (/ib,jb,kb/)
       !                      idx2 = modulo(idx2,f(fid)%nwan)
       
       !                      asum = asum + sum(f(fid)%fwan(idx1(1)*n(1)+1:idx1(1)*n(1)+n(1),idx1(2)*n(2)+1:idx1(2)*n(2)+n(2),idx1(3)*n(3)+1:idx1(3)*n(3)+n(3),iba,is) * &
       !                         f(fid)%fwan(idx2(1)*n(1)+1:idx2(1)*n(1)+n(1),idx2(2)*n(2)+1:idx2(2)*n(2)+n(2),idx2(3)*n(3)+1:idx2(3)*n(3)+n(3),ibb,is))
       !                   end do
       !                end do
       !             end do
       !             asum = asum / (n(1)*n(2)*n(3))
       !             write (uout,'(99(A))') "b=", string(ibb), " Rx=", string(ib), &
       !                " Ry=", string(jb), " Rz=", string(kb), " ", string(asum,'e')
       !          end do
       !       end do
       !    end do
       ! else
       !    write (uout,*) "field 2 does not have wannier functions"
       !    stop 1
       ! end if

       ! ! check -- is the density correct?
       ! if (refden == 1 .and. fid == 2) then
       !    write (uout,'("+ Density check +")')
       !    write (uout,'("Maximum rho diff: ",A)') string(maxval(f(fid)%f - f(refden)%f),'e')
       !    write (uout,'("Minimum rho diff: ",A)') string(minval(f(fid)%f - f(refden)%f),'e')
       !    write (uout,'("Average rho diff: ",A)') string(sum(f(fid)%f - f(refden)%f) / (f(fid)%n(1)*f(fid)%n(2)*f(fid)%n(3)),'e')
       !    write (uout,'("Integral 1: ",A)') string(sum(f(fid)%f) * cr%omega / (f(fid)%n(1)*f(fid)%n(2)*f(fid)%n(3)),'e')
       !    write (uout,'("Integral 2: ",A)') string(sum(f(refden)%f) * cr%omega / (f(refden)%n(1)*f(refden)%n(2)*f(refden)%n(3)),'e')
       ! else
       !    write (uout,*) "rho needs to be field 1 and xsf field 2"
       !    stop 1
       ! endif


       ! write (uout,'("+ Charge check using the overlap matrices")')
       ! do is = 1, nspin
       !    write (uout,'(" Spin = ",A)') string(is)
       !    asum = 0d0
       !    do imo = 1, nmo
       !       asum = asum + sum(sij(imo,imo,:,:,ndeloc)) * fspin
       !    end do
       !    write (uout,'(" N(1) -- sum_Rn wRn^2 for atom ",I2,X,F12.6)') 1, asum 
       ! end do
       
       ! allocate(oij(nmo,nmo,nspin))
       ! oij = 0d0
       ! write (uout,'("+ Orthonormality check")')
       ! do is = 1, nspin
       !    write (uout,'(" Spin = ",A)') string(is)
       !    do imo = 1, nmo
       !       call unpackidx(imo,ia,ja,ka,iba)
       !       do jmo = 1, nmo
       !          call unpackidx(jmo,ib,jb,kb,ibb)
       !          do ic = i0, i1
       !             do jc = j0, j1
       !                do kc = k0, k1
       !                   idx1 = (/ic+ia-ib,jc+ja-jb,kc+ka-kb/)
       !                   idx1 = modulo(idx1,f(fid)%nwan)
       !                   id = idx1(1)
       !                   jd = idx1(2)
       !                   kd = idx1(3)
       !                   call packidx(id,jd,kd,iba,kmo)
       !                   call packidx(ic,jc,kc,ibb,lmo)
       !                   oij(kmo,lmo,is) = oij(kmo,lmo,is) + sum(sij(imo,jmo,:,is,ndeloc))
       !                end do
       !             end do
       !          end do
       !       end do
       !    end do
       !    do imo = 1, nmo
       !       do jmo = 1, nmo
       !          write (uout,'("Orb. ",I2,X,I2,X,F12.6)') imo, jmo, oij(imo,jmo,1)
       !       end do
       !    end do
       ! end do
       ! write (uout,*)

