    ! attic !

    ! ! Direct space crystal point group (deactivated)
    ! vec = 0d0
    ! call lattpg(transpose(c%crys2car),1,vec)
    ! c%pointgroup = trim(point_group)

    ! classify the lattice (deactivated)
    ! if (c%havesym >= 1) call c%classify()

     !! Attic !!
     ! Symmetry classification (classify is deactivated)
     ! character*2 :: delaunay !< Delaunay symbol (table 9.1.8.1, ITC)
     ! character*2 :: bravais_type !< Bravais type (table 9.1.8.1, ITC)
     ! character*1 :: cfam !< Crystal family 
     ! integer :: delsort !< Delaunay sort (row in table 9.1.8.1, ITC)
     ! character*3 :: pointgroup !< Crystal point group
     ! real*8 :: rmat_conventional(3,3) !< Transformation to conventional cell
     ! 1=1bar, 2=2/m, 3=mmm, 4=4/m, 5=4/mmm, 6=3bar, 7=3bar/m, 8=6/m, 
     ! 9=6/mmm, 10=m3bar, 11=m3barm
     ! ws cell neighbor information

  !> Classify the lattice using Delaunay reduction. (deactivated, 
  !> see 002_read_crinput6)
  ! subroutine classify(c)
  !   class(crystal), intent(inout) :: c
  ! 
  !   real*8 :: rmat(3,3), dmat(3,4), sc(4,4), scv(6), xn(4)
  !   integer :: nzero, ndiff, nzero123, nsame(6)
  ! 
  !   integer :: i
  !   logical :: done(6)
  ! 
  !   real*8, parameter :: eps = 1d-10
  ! 
  !   ! run the delaunay reduction of the primitive cell
  !   if (c%havesym >= 1) then
  !      call c%primitive_any(.false.,rmat)
  !      call c%delaunay_reduction(dmat,rmat,sc)
  !   else
  !      c%delaunay = "??"
  !      c%bravais_type = "??"
  !      c%cfam = "?"
  !      c%delsort = 0
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      return
  !   end if
  ! 
  !   ! unpack the scalar product matrix
  !   scv(1) = sc(1,2)
  !   scv(2) = sc(1,3)
  !   scv(3) = sc(1,4)
  !   scv(4) = sc(2,3)
  !   scv(5) = sc(2,4)
  !   scv(6) = sc(3,4)
  ! 
  !   ! implementation of the classification in table 9.1.8.1 ITC
  !   ! count number of zeros and number of distinct elements
  !   nzero = count(abs(scv) < eps)
  !   nzero123 = count(abs(scv((/1,2,4/))) < eps)
  !   done = (abs(scv) < eps)
  !   ndiff = 0
  !   do i = 1, 6
  !      if (.not.done(i)) then
  !         ndiff = ndiff + 1
  !         nsame(ndiff) = count(abs(scv - scv(i)) < eps)
  !         done = done .or. (abs(scv - scv(i)) < eps)
  !      end if
  !   end do
  !   
  !   ! Classify, from the top of the table. I'm having trouble because
  !   ! the vectors do not come in order out of the Delaunay reduction
  !   ! procedure, and it's not clear in ITC how to choose the reduced
  !   ! basis.
  !   if (nzero == 0 .and. ndiff == 1) then
  !      ! row=1, K1, cI: 12 12 12 12 12 12 
  !      c%delaunay = "K1"
  !      c%bravais_type = "cI"
  !      c%cfam = "c"
  !      c%delsort = 1
  !      c%rmat_conventional(:,1) = (/ 0d0, 1d0, 1d0 /)
  !      c%rmat_conventional(:,2) = (/ 1d0, 0d0, 1d0 /)
  !      c%rmat_conventional(:,3) = (/ 1d0, 1d0, 0d0 /)
  !   elseif (nzero == 2 .and. ndiff == 1) then
  !      ! row=2, K2, cF: 0 13 13 13 13 0
  !      c%delaunay = "K2"
  !      c%bravais_type = "cF"
  !      c%cfam = "c"
  !      c%delsort = 2
  !      c%rmat_conventional(:,1) = (/ 1d0, -1d0, 1d0 /)
  !      c%rmat_conventional(:,2) = (/ 1d0, 1d0, 1d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 2d0 /)
  !   elseif (nzero == 3 .and. ndiff == 1) then
  !      if (abs(scv(4)) > eps) then
  !         ! row=3, K3, cP: 0 0 14 14 14 0
  !         c%delaunay = "K3"
  !         c%bravais_type = "cP"
  !         c%cfam = "c"
  !         c%delsort = 3
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 1d0, 1d0 /)
  !      else
  !         ! row=4, K3, cP: 0 0 14 0 14 14 0 
  !         c%delaunay = "K3"
  !         c%bravais_type = "cP"
  !         c%cfam = "c"
  !         c%delsort = 4
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      end if
  !   elseif (nzero == 2 .and. ndiff == 2) then
  !      if (nzero123 == 2) then
  !         ! row=5, H , hP: 12 0 12 0 12 34
  !         c%delaunay = "H "
  !         c%bravais_type = "hP"
  !         c%cfam = "h"
  !         c%delsort = 5
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      elseif (nsame(1) == 3 .or. nsame(2) == 3) then
  !         ! row=7, R2, hR: 0 13 13 13 24 0
  !         c%delaunay = "R2"
  !         c%bravais_type = "hR"
  !         c%cfam = "h"
  !         c%delsort = 7
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 0d0, 3d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 1d0, 2d0 /)
  !      elseif (abs(scv(3) - scv(4)) < eps) then
  !         ! row=16, O4, oI: 0 13 14 14 13 0
  !         c%delaunay = "O4"
  !         c%bravais_type = "oI"
  !         c%cfam = "o"
  !         c%delsort = 16
  !         c%rmat_conventional(:,1) = (/ 0d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 1d0, 1d0, 0d0 /)
  !      else
  !         ! row=17, O4, oI: 0 13 13 23 23 0
  !         c%delaunay = "O4"
  !         c%bravais_type = "oI"
  !         c%cfam = "o"
  !         c%delsort = 17
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 2d0 /)
  !      end if
  !   elseif (nzero == 0 .and. ndiff == 2) then
  !      if (nsame(1) == 3 .or. nsame(2) == 3) then
  !         ! row=6, R1, hR: 12 12 14 12 14 14
  !         c%delaunay = "R1"
  !         c%bravais_type = "hR"
  !         c%cfam = "h"
  !         c%delsort = 6
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ -1d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, -1d0, 1d0 /)
  !      else
  !         ! row=8, Q1, tI: 12 13 13 13 13 12
  !         c%delaunay = "Q1"
  !         c%bravais_type = "tI"
  !         c%cfam = "t"
  !         c%delsort = 8
  !         c%rmat_conventional(:,1) = (/ 0d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 1d0, 1d0, 0d0 /)
  !      end if
  !   elseif (nzero == 1 .and. (ndiff == 2 .or. ndiff == 1)) then
  !      ! row=9, Q2, tI: 0 13 13 13 13 34
  !      c%delaunay = "Q2"
  !      c%bravais_type = "tI"
  !      c%cfam = "t"
  !      c%delsort = 9
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 1d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 1d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 2d0 /)
  !   elseif (nzero == 3 .and. ndiff == 2) then
  !      if (nzero123 == 3) then
  !         ! row=10, Q3, tP: 0 0 14 0 14 34
  !         c%delaunay = "Q3"
  !         c%bravais_type = "tP"
  !         c%cfam = "t"
  !         c%delsort = 10
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      elseif (abs(scv(6)) < eps) then
  !         ! row=11, Q3, tP: 0 0 14 14 24 0
  !         c%delaunay = "Q3"
  !         c%bravais_type = "tP"
  !         c%cfam = "t"
  !         c%delsort = 11
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 1d0, 1d0 /)
  !      else
  !         ! row=12, Q3, tP: 0 0 14 23 0 23
  !         c%delaunay = "Q3"
  !         c%bravais_type = "tP"
  !         c%cfam = "t"
  !         c%delsort = 12
  !         c%rmat_conventional(:,1) = (/ 0d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 1d0, 0d0 /)
  !      end if
  !   elseif (nzero == 0 .and. ndiff == 3) then
  !      if (any(nsame(1:3) == 4)) then
  !         ! row=13, O1, oF: 12 13 13 13 13 34
  !         c%delaunay = "O1"
  !         c%bravais_type = "oF"
  !         c%cfam = "o"
  !         c%delsort = 13
  !         c%rmat_conventional(:,1) = (/ 1d0, -1d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 2d0 /)
  !      else
  !         ! row=14, O2, oI: 12 13 14 14 13 12
  !         c%delaunay = "O2"
  !         c%bravais_type = "oI"
  !         c%cfam = "o"
  !         c%delsort = 14
  !         c%rmat_conventional(:,1) = (/ 0d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 1d0, 1d0, 0d0 /)
  !      end if
  !   elseif (nzero == 1 .and. ndiff == 3) then
  !      if (abs(scv(2) - scv(3)) < eps) then
  !         ! row=15, O3, oI: 0 13 13 23 23 34
  !         c%delaunay = "O3"
  !         c%bravais_type = "oI"
  !         c%cfam = "o"
  !         c%delsort = 15
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 1d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 2d0 /)
  !      elseif (abs(scv(2) - scv(5)) < eps) then
  !         ! row=25, M4, mI: 0 13 14 14 13 34
  !         c%delaunay = "M4"
  !         c%bravais_type = "mI"
  !         c%cfam = "o"
  !         c%delsort = 25
  !         c%rmat_conventional(:,1) = (/ 0d0, 1d0, -1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 1d0, 0d0, -1d0 /)
  !      else
  !         ! row=26, M4, mI: 0 13 14 13 14 34
  !         c%delaunay = "M4"
  !         c%bravais_type = "mI"
  !         c%cfam = "o"
  !         c%delsort = 26
  !         c%rmat_conventional(:,1) = (/ -1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ -1d0, -1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ -1d0, 0d0, 1d0 /)
  !      end if
  !   elseif (nzero == 2 .and. ndiff == 3) then
  !      if (abs(scv(1) - scv(5)) < eps) then
  !         c%delaunay = "O5"
  !         ! row=18, O5, oC: 12 0 14 0 12 34
  !         c%bravais_type = "oC"
  !         c%cfam = "o"
  !         c%delsort = 18
  !         c%rmat_conventional(:,1) = (/ 2d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      elseif (abs(scv(3) - scv(5)) < eps) then
  !         ! row=19, O5, oC: 12 0 14 0 14 34
  !         c%delaunay = "O5"
  !         c%bravais_type = "oC"
  !         c%cfam = "o"
  !         c%delsort = 19
  !         c%rmat_conventional(:,1) = (/ 1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ -1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      elseif (abs(scv(4) - scv(5)) < eps) then
  !         ! row=27, M5, mI: 0 13 14 23 23 0
  !         c%delaunay = "M5"
  !         c%bravais_type = "mI"
  !         c%cfam = "m"
  !         c%delsort = 27
  !         c%rmat_conventional(:,1) = (/ -1d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,2) = (/ -1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ -2d0, 0d0, 0d0 /)
  !      else
  !         ! row=28, M5, mI: 0 13 14 23 13 0
  !         c%delaunay = "M5"
  !         c%bravais_type = "mI"
  !         c%cfam = "m"
  !         c%delsort = 28
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, -1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, -1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, -1d0, -1d0 /)
  !      end if
  !   elseif (nzero == 3 .and. ndiff == 3) then
  !      if (abs(scv(4)) < eps) then
  !         ! row=20, O6, oP: 0 0 14 0 24 34
  !         c%delaunay = "O6"
  !         c%bravais_type = "oP"
  !         c%cfam = "o"
  !         c%delsort = 20
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      else
  !         ! row=21, O6, oP: 0 0 14 23 24 0
  !         c%delaunay = "O6"
  !         c%bravais_type = "oP"
  !         c%cfam = "o"
  !         c%delsort = 21
  !         c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ 0d0, 0d0, 1d0 /)
  !         c%rmat_conventional(:,3) = (/ 0d0, 1d0, 1d0 /)
  !      end if
  !   elseif (nzero == 0 .and. ndiff == 4) then
  !      if (abs(scv(2) - scv(4)) < eps) then
  !         ! row=22, M1, mI: 12 13 14 13 14 34
  !         c%delaunay = "M1"
  !         c%bravais_type = "mI"
  !         c%cfam = "m"
  !         c%delsort = 22
  !         c%rmat_conventional(:,1) = (/ -1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,2) = (/ -1d0, -1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ -1d0, 0d0, 1d0 /)
  !      else
  !         ! row=23, M2, mI: 12 13 14 14 13 34
  !         c%delaunay = "M2"
  !         c%bravais_type = "mI"
  !         c%cfam = "m"
  !         c%delsort = 23
  !         c%rmat_conventional(:,1) = (/ 0d0, 1d0, -1d0 /)
  !         c%rmat_conventional(:,2) = (/ 1d0, 1d0, 0d0 /)
  !         c%rmat_conventional(:,3) = (/ 1d0, 0d0, -1d0 /)
  !      end if
  !   elseif (nzero == 1 .and. ndiff == 4) then
  !      ! row=24, M3, mI: 0 13 14 23 23 34
  !      c%delaunay = "M3"
  !      c%bravais_type = "mI"
  !      c%cfam = "m"
  !      c%delsort = 24
  !      c%rmat_conventional(:,1) = (/ -1d0, 0d0, 1d0 /)
  !      c%rmat_conventional(:,2) = (/ -1d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ -2d0, 0d0, 0d0 /)
  !   elseif (nzero == 2 .and. ndiff == 4) then
  !      ! row=29, M6, mP: 0 13 14 0 24 34
  !      c%delaunay = "M6"
  !      c%bravais_type = "mP"
  !      c%cfam = "m"
  !      c%delsort = 29
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !   elseif (nzero == 0 .and. ndiff == 6) then
  !      ! row=30, T1, aP: 12 13 14 23 24 34
  !      c%delaunay = "T1"
  !      c%bravais_type = "aP"
  !      c%cfam = "a"
  !      c%delsort = 30
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !   elseif (nzero == 1 .and. ndiff == 5) then
  !      ! row=31, T2, aP: 0 13 14 23 24 34
  !      c%delaunay = "T2"
  !      c%bravais_type = "aP"
  !      c%cfam = "a"
  !      c%delsort = 31
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !   elseif (nzero == 2 .and. ndiff == 4) then
  !      ! row=32, T3, aP: 0 13 14 23 24 0
  !      c%delaunay = "T3"
  !      c%bravais_type = "aP"
  !      c%cfam = "a"
  !      c%delsort = 32
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !   else
  !      c%delaunay = "??"
  !      c%bravais_type = "??"
  !      c%cfam = "?"
  !      c%delsort = 0
  !      c%rmat_conventional(:,1) = (/ 1d0, 0d0, 0d0 /)
  !      c%rmat_conventional(:,2) = (/ 0d0, 1d0, 0d0 /)
  !      c%rmat_conventional(:,3) = (/ 0d0, 0d0, 1d0 /)
  !      call ferror('classify','could not classify lattice delsort',warning)
  !   end if
  ! 
  ! end subroutine classify

  ! !> Transform to the standard conventional cell. (???)
  ! subroutine conventional_standard(c,verbose,rmat)
  !   class(crystal), intent(inout) :: c
  !   logical, intent(in) :: verbose
  !   real*8, intent(out), optional :: rmat(3,3)
  ! 
  !   real*8 :: m(3,3)
  ! 
  !   real*8, parameter :: eps = 1d-6
  ! 
  !   if (present(rmat)) rmat = eye
  ! 
  !   ! ignore molecules
  !   if (c%ismolecule) return
  ! 
  !   ! Only available if havesym >= 1
  !   if (c%havesym < 1) return
  ! 
  !   ! transform to the primitive or output through rmat
  !   if (present(rmat)) then
  !      rmat = c%rmat_conventional
  !   else
  !      call c%primitive_delaunay(.true.)
  !      call c%newcell(c%rmat_conventional,verbose0=.true.)
  !   end if
  ! 
  ! end subroutine conventional_standard
