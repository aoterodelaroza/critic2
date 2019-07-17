  !> Find asterisms. For every atom in the unit cell, find the atoms in the 
  !> main cell or adjacent cells that are connected to it. 
  module subroutine find_asterisms(c)
    use global, only: bondfactor
    use param, only: atmcov, vsmall
    use types, only: realloc

    class(crystal), intent(inout) :: c

    integer :: i, j, k
    real*8 :: rvws(3), x0(3), dist, dist2
    real*8 :: d0
    integer :: lvec0(3), lvec(3)

    if (allocated(c%nstar)) deallocate(c%nstar)
    if (.not.allocated(c%nstar)) allocate(c%nstar(c%ncel))

    ! allocate the neighbor star
    do i = 1, c%ncel
       allocate(c%nstar(i)%idcon(4))
       allocate(c%nstar(i)%lcon(3,4))
    end do

    if (c%ismolecule) then
       ! run over all pairs of atoms in the molecule
       lvec = 0
       do i = 1, c%ncel
          do j = i+1, c%ncel
             d0 = bondfactor * (atmcov(c%spc(c%atcel(i)%is)%z)+atmcov(c%spc(c%atcel(j)%is)%z))
             ! use the Cartesian directly
             x0 = c%atcel(j)%r - c%atcel(i)%r
             if (any(abs(x0) > d0)) cycle
             d0 = d0 * d0
             dist2 = x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3)
             if (dist2 < d0) then
                call addpair(i,j,lvec)
                call addpair(j,i,lvec)
             end if
          end do
       end do
    else
       ! run over all pairs of atoms in the unit cell
       do i = 1, c%ncel
          do j = i, c%ncel
             x0 = c%atcel(j)%x - c%atcel(i)%x
             lvec0 = nint(x0)
             x0 = x0 - lvec0
             d0 = bondfactor * (atmcov(c%spc(c%atcel(i)%is)%z)+atmcov(c%spc(c%atcel(j)%is)%z))

             do k = 0, c%ws_nf
                if (k == 0) then
                   rvws = x0
                   lvec = lvec0
                else
                   rvws = x0 - c%ws_ineighx(:,k)
                   lvec = lvec0 + c%ws_ineighx(:,k)
                endif
                rvws = matmul(c%m_x2c,rvws)
                dist = sqrt(rvws(1)*rvws(1)+rvws(2)*rvws(2)+rvws(3)*rvws(3))
                if (all(abs(rvws) < d0+1d-6)) then
                   dist = sqrt(rvws(1)*rvws(1)+rvws(2)*rvws(2)+rvws(3)*rvws(3))
                   if (dist > vsmall .and. dist < d0) then
                      call addpair(i,j,lvec)
                      call addpair(j,i,-lvec)
                   end if
                end if
             end do
          end do
       end do
    end if
    do i = 1, c%ncel
       call realloc(c%nstar(i)%idcon,c%nstar(i)%ncon)
       call realloc(c%nstar(i)%lcon,3,c%nstar(i)%ncon)
    end do

  contains
    subroutine addpair(i,j,lvec)
      integer :: i, j, lvec(3)

      c%nstar(i)%ncon = c%nstar(i)%ncon + 1
      if (c%nstar(i)%ncon > size(c%nstar(i)%idcon)) then
         call realloc(c%nstar(i)%idcon,2*c%nstar(i)%ncon)
         call realloc(c%nstar(i)%lcon,3,2*c%nstar(i)%ncon)
      end if
      c%nstar(i)%idcon(c%nstar(i)%ncon) = j
      c%nstar(i)%lcon(:,c%nstar(i)%ncon) = -lvec

    end subroutine addpair
  end subroutine find_asterisms

