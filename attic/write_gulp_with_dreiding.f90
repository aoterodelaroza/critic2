  !> Write a gulp input script
  module subroutine write_gulp(c,file,dodreiding)
    use tools_io, only: fopen_write, faterr, ferror, nameguess, fclose, string
    use param, only: bohrtoa, atmcov, pi
    class(crystal), intent(inout) :: c
    character*(*), intent(in) :: file
    logical :: dodreiding

    integer, parameter :: maxneigh = 20
    real*8, parameter :: rfac = 1.4d0
    real*8, parameter :: hbmin = 1.6 / bohrtoa
    real*8, parameter :: hbmax = 3.2 / bohrtoa

    character*(5) :: lbl
    integer :: lu, i, j, iz, jz, n, k, kz
    integer :: idx
    integer :: nneigh(maxneigh), ineigh(maxneigh,c%nneq)
    real*8 :: dneigh(maxneigh,c%nneq), dhb(maxneigh,c%nneq), avgang(c%nneq)
    integer :: nhb(maxneigh), ihb(maxneigh,c%nneq)
    real*8 :: d, x1(3), x2(3), ang, xi(3)
    logical :: ok, isat

    ! check that we have an environment
    call c%checkflags(.true.,env0=.true.)

    lu = fopen_write(file)
    if (.not. dodreiding) then
       write (lu,'("eem")')
       write (lu,'("cell ",6(A,X))') (string(c%aa(j) * bohrtoa,'f',13,9),j=1,3), &
          (string(c%bb(j),'f',10,5),j=1,3)
       write (lu,'("fractional")')
       do i = 1, c%ncel
          write (lu,'(A5,X,3(A,X))') trim(c%spc(c%atcel(i)%is)%name),&
             (string(c%atcel(i)%x(j),'f',15,9),j=1,3)
       end do
    else
       ! calculate bonded neighbors
       nneigh = 0
       nhb = 0
       do i = 1, c%nneq
          xi = c%x2xr(c%at(i)%x)
          xi = xi - floor(xi)
          xi = c%xr2c(xi)

          iz = c%spc(c%at(i)%is)%z
          n = 0
          ! determine covalent bonds
          do j = 1, c%env%n
             jz = c%spc(c%env%at(j)%is)%z
             d = norm2(c%env%at(j)%r-xi)
             if (d < 1d-10) cycle
             if (d < (atmcov(iz) + atmcov(jz)) * rfac) then
                n = n + 1
                if (n > maxneigh) call ferror("write_gulp","too many neighbors",faterr)
                ineigh(n,i) = j
                dneigh(n,i) = d
                continue
             endif
          end do
          nneigh(i) = n
          ! determine average angles with bonded neighbors
          avgang(i) = 0d0
          n = 0
          do j = 1, nneigh(i)
             do k = j+1, nneigh(i)
                n = n + 1
                x1 = c%env%at(ineigh(j,i))%r - xi
                x2 = c%env%at(ineigh(k,i))%r - xi
                ang = abs(acos(dot_product(x1,x2) / norm2(x1) / norm2(x2)) * 180d0 / pi)
                avgang(i) = avgang(i) + ang
             end do
          end do
          if (n > 0) avgang(i) = avgang(i) / n

          ! determine hydrogen bonds, only for hydrogen
          if (iz == 1) then
             n = 0
             do j = 1, c%env%n
                jz = c%spc(c%env%at(j)%is)%z
                ! only with N, O, and S
                if (jz==7 .or. jz==8 .or. jz==9 .or. jz==16 .or. jz==17 .or. jz==35 .or. jz==53) then
                   d = norm2(c%env%at(j)%r-xi)
                   ! only in the correct distance range
                   if (d > hbmin .and. d < hbmax) then
                      ! only if the angles to all other neighbor atoms is more than 145
                      ok = .true.
                      do k = 1, nneigh(i)
                         x1 = c%env%at(ineigh(k,i))%r - xi
                         x2 = c%env%at(j)%r - xi
                         ang = abs(acos(dot_product(x1,x2) / norm2(x1) / norm2(x2)) * 180d0 / pi)
                         kz = c%spc(c%env%at(ineigh(k,i))%is)%z
                         isat = (kz==7 .or. kz==8 .or. kz==9 .or. kz==16 .or. kz==17 .or. kz==35 .or. kz==53)
                         if (ang < 145d0 .or..not.isat) then
                            ok = .false.
                            exit
                         end if
                      end do
                      ! oh, sure, fine... you're a hydrogen bond
                      if (ok) then
                         n = n + 1
                         ihb(n,i) = j
                         dhb(n,i) = d
                      endif
                   end if
                end if
             end do
             nhb(i) = n
          endif
       end do

       write (lu,'("eem")')
       write (lu,'("cell ",6(A,X))') (string(c%aa(j) * bohrtoa,'f',13,9),j=1,3), &
          (string(c%bb(j),'f',10,5),j=1,3)
       write (lu,'("fractional")')
       do i = 1, c%ncel
          idx = c%atcel(i)%idx
          iz = c%spc(c%at(idx)%is)%z
          ang = avgang(idx)
          ! the first two letters is the atomic symbol
          lbl = adjustl(nameguess(iz))
          ! hydrogen types: H_ (normal), H___A (hydrogen-bonded to N, O, or S), H___b (bridging)
          if (iz == 1 .and. nneigh(idx) > 1) lbl = "H___b"
          if (iz == 1 .and. nhb(idx) > 0) lbl = "H___A"
          ! boron: sp3 (109.47 angles) and sp2 (120 angles)
          if (iz == 5 .and. abs(ang-109.47d0)<abs(ang-120d0)) lbl = "B_3"
          if (iz == 5 .and. abs(ang-109.47d0)>abs(ang-120d0)) lbl = "B_2"
          ! carbon: sp3 (109.47 angles), sp2 (120 angles), sp (180 angles)
          if (iz == 6 .and. abs(ang-109.47d0)<abs(ang-120d0) .and. abs(ang-109.47d0)<abs(ang-180d0)) lbl = "C_3"
          if (iz == 6 .and. abs(ang-120d0)<abs(ang-109.47d0) .and. abs(ang-120d0)<abs(ang-180d0)) lbl = "C_2"
          if (iz == 6 .and. abs(ang-180d0)<abs(ang-109.47d0) .and. abs(ang-180d0)<abs(ang-120d0)) lbl = "C_1"
          ! nitrogen: sp3 (109.47 angles), sp2 (120 angles), sp (180 angles)
          if (iz == 7 .and. abs(ang-109.47d0)<abs(ang-120d0) .and. abs(ang-109.47d0)<abs(ang-180d0)) lbl = "N_3"
          if (iz == 7 .and. abs(ang-120d0)<abs(ang-109.47d0) .and. abs(ang-120d0)<abs(ang-180d0)) lbl = "N_2"
          if (iz == 7 .and. (abs(ang-180d0)<abs(ang-109.47d0) .and. abs(ang-180d0)<abs(ang-120d0) .or. nneigh(idx) == 1)) lbl = "N_1"
          ! oxygen: sp3 (109.47 angles), sp2 (120 angles), sp (180 angles)
          if (iz == 8 .and. abs(ang-109.47d0)<abs(ang-120d0) .and. abs(ang-109.47d0)<abs(ang-180d0)) lbl = "O_3"
          if (iz == 8 .and. (abs(ang-120d0)<abs(ang-109.47d0) .and. abs(ang-120d0)<abs(ang-180d0) .or. nneigh(idx) == 1)) lbl = "O_2"
          if (iz == 8 .and. abs(ang-180d0)<abs(ang-109.47d0) .and. abs(ang-180d0)<abs(ang-120d0)) lbl = "O_1"
          ! Al, Si, P, S, Ga, Ge, As, Se, In, Sn, Sb, Te -> only sp3 is known
          if (iz == 13 .or. iz == 14 .or. iz == 15 .or. iz == 16 .or. iz == 31 .or. iz == 32 .or. iz == 33 .or. iz == 34 .or.&
             iz == 49 .or. iz == 50 .or. iz == 51 .or. iz == 52) then
             lbl(3:3) = "3"
          end if
          write (lu,'(A5,X,3(A,X))') adjustl(trim(lbl)), (string(c%atcel(i)%x(j),'f',15,9),j=1,3)
       end do
    end if
    call fclose(lu)

  end subroutine write_gulp
