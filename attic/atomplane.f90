  !> Calculate the axis of the plane determined by 3 atoms. 
  subroutine atomplane(i,j,k,scalx,scaly,sizex,sizey,p0,p1,p2)
    use struct_basic
    use global
    use tools_io

    integer, intent(in) :: i, j, k
    real*8, intent(in) :: scalx, scaly, sizex, sizey
    real*8, dimension(3), intent(out) :: p0, p1, p2

    real*8, allocatable :: vec(:,:)
    real*8, dimension(3) :: v2, x
    real*8 :: dist, dmax, det
    integer :: m, l, mult
    integer :: lvec(3), ln

    p0 = cr%atcel(i)%x
    p1 = cr%atcel(j)%x
    call cr%symeqv(p1,mult,vec)
    dmax = 1d10
    do l = 1, mult
       do m = 0, 26
          ln = m
          lvec(1) = mod(ln,3) - 1
          ln = ln / 3
          lvec(2) = mod(ln,3) - 1
          ln = ln / 3
          lvec(3) = mod(ln,3) - 1

          x = cr%x2c(vec(:,l) - p0 + real(lvec,8))
          dist = dot_product(x,x)
          if (dist < 1d-8) cycle
          if (dist < dmax) then
             p1 = vec(:,l) + real(lvec,8)
             dmax = dist
          end if
       end do
    end do

    p2 = cr%atcel(k)%x
    call cr%symeqv(p2,mult,vec)
    dmax = 1d10
    do l = 1, mult
       do m = 0, 26
          ln = m
          lvec(1) = mod(ln,3) - 1
          ln = ln / 3
          lvec(2) = mod(ln,3) - 1
          ln = ln / 3
          lvec(3) = mod(ln,3) - 1
          
          x = cr%x2c(vec(:,l) - p0 + real(lvec,8))
          dist = dot_product(x,x)
          if (dist < 1d-8) cycle
          if (dist < dmax) then
             ! check determinant
             v2 = vec(:,l) + real(lvec,8)
             det = p0(1)*p1(2)*v2(3) + p0(2)*p1(3)*v2(1) + &
                  p0(3)*p1(1)*v2(2) - p0(3)*p1(2)*v2(1) - &
                  p0(2)*p1(1)*v2(3) - p0(1)*p1(3)*v2(2)
             if (abs(det) < 1d-5) cycle
             p2 = v2
             dmax = dist
          end if
       end do
    end do
    deallocate(vec)
    call buildplane(p0,p1,p2,scalx,scaly,sizex,sizey)

  end subroutine atomplane

