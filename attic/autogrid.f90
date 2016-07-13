  ! doesn't work -- missing points in promolecular urea.
  subroutine autogrid()
    use fields
    use struct_basic
    use global
    use tools_io

    type(crystal) :: caux

    integer :: ix, iy, iz, jx, jy, jz, kx, ky, kz, iv
    integer :: i, j, nvec, nvec0, vec(3,16), n(3)
    real*8 :: al(40), x0(3)
    logical :: ifound(16), ok, iscp

    ! is the reference field a grid? 
    if (f(refden)%type /= type_grid) &
       call ferror("autogrid","autogrid can only be used with fields on a grid",faterr)
    n = f(refden)%n

    ! calculate areas*lengths and grid vectors. Use a smaller crystal
    ! where the "lattice" corresponds to the cube grid points.
    caux%isinit = .true.
    caux%aa = cr%aa / real(n,8)
    caux%bb = cr%bb
    call caux%set_cryscar()
    call caux%wigner((/0d0,0d0,0d0/),.false.,nvec,vec,al)

    ! only vectors not related by inversion
    ifound = .false.
    do i = 1, nvec
       if (ifound(i)) cycle
       ok = .false.
       do j = i+1, nvec
          if (all(vec(:,i) + vec(:,j) == 0)) then
             ifound(j) = .true.
             ok = .true.
             exit
          end if
       end do
       if (.not.ok) &
          call ferror("autogrid","lattice vector pair not found",faterr)
    end do

    nvec0 = 0
    do i = 1, nvec
       if (.not.ifound(i)) then
          nvec0 = nvec0 + 1
          vec(:,nvec0) = vec(:,i)
       end if
    end do
    nvec = nvec0

    do ix = 1, n(1)
       do iy = 1, n(2)
          do iz = 1, n(3)
             iscp = .true.
             do iv = 1, nvec
                call translate(ix,iy,iz, vec(1,iv), vec(2,iv), vec(3,iv),jx,jy,jz)
                call translate(ix,iy,iz,-vec(1,iv),-vec(2,iv),-vec(3,iv),kx,ky,kz)
                if (f(refden)%f(jx,jy,jz) < f(refden)%f(ix,iy,iz) .and. f(refden)%f(ix,iy,iz) < f(refden)%f(kx,ky,kz) .or.&
                    f(refden)%f(jx,jy,jz) > f(refden)%f(ix,iy,iz) .and. f(refden)%f(ix,iy,iz) > f(refden)%f(kx,ky,kz)) then
                   iscp = .false.
                   exit
                end if
             end do
             if (iscp) then
                x0 = real((/ix, iy, iz/)-1,8) / n
                x0 = cr%x2c(x0)
                call addcp(x0)
             end if
          end do
       end do
    end do

  contains
    ! Translate point x0 by xd and give the result back in x1.
    ! Uses the n(:) from the main routine
    subroutine translate(i0,j0,k0,id,jd,kd,i1,j1,k1)
      integer :: i0, j0, k0, id, jd, kd, i1, j1, k1

      i1 = modulo(i0 + id - 1,n(1)) + 1
      j1 = modulo(j0 + jd - 1,n(2)) + 1
      k1 = modulo(k0 + kd - 1,n(3)) + 1

    end subroutine translate
  end subroutine autogrid

