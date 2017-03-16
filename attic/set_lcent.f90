  !> Calculate the lcent from the centering vectors (ncv and cen)
  subroutine set_lcent(c)
    class(crystal), intent(inout) :: c

    logical :: ok

    c%lcent = 9 ! unknown
    if (c%ncv == 0) then
       c%lcent = 0 ! unset
    elseif (c%ncv == 1) then
       c%lcent = 1 ! P
    elseif (c%ncv == 2) then
       if (c%are_lclose(c%cen(:,2),(/0.0d0,0.5d0,0.5d0/),1d-5)) then
          c%lcent = 2 ! A
       elseif (c%are_lclose(c%cen(:,2),(/0.5d0,0.0d0,0.5d0/),1d-5)) then
          c%lcent = 3 ! B
       elseif (c%are_lclose(c%cen(:,2),(/0.5d0,0.5d0,0.0d0/),1d-5)) then
          c%lcent = 4 ! C
       elseif (c%are_lclose(c%cen(:,2),(/0.5d0,0.5d0,0.5d0/),1d-5)) then
          c%lcent = 6 ! I
       end if
    elseif (c%ncv == 3) then
       if (c%are_lclose(c%cen(:,2),(/2d0/3d0,1d0/3d0,1d0/3d0/),1d-5) .and.&
           c%are_lclose(c%cen(:,3),(/-2d0/3d0,-1d0/3d0,-1d0/3d0/),1d-5) .or.&
           c%are_lclose(c%cen(:,2),(/-2d0/3d0,-1d0/3d0,-1d0/3d0/),1d-5) .and.&
           c%are_lclose(c%cen(:,3),(/2d0/3d0,1d0/3d0,1d0/3d0/),1d-5)) then
           c%lcent = 7 ! R (obverse)
        elseif (c%are_lclose(c%cen(:,2),(/1d0/3d0,2d0/3d0,1d0/3d0/),1d-5) .and.&
           c%are_lclose(c%cen(:,3),(/-1d0/3d0,-2d0/3d0,-1d0/3d0/),1d-5) .or.&
           c%are_lclose(c%cen(:,2),(/-1d0/3d0,-2d0/3d0,-1d0/3d0/),1d-5) .and.&
           c%are_lclose(c%cen(:,3),(/1d0/3d0,2d0/3d0,1d0/3d0/),1d-5)) then
           c%lcent = 8 ! R (reverse)
        end if
    elseif (c%ncv == 4) then
       ok = c%are_lclose(c%cen(:,2),(/0d0,0.5d0,0.5d0/),1d-5) .or.&
          c%are_lclose(c%cen(:,2),(/0.5d0,0d0,0.5d0/),1d-5) .or.&
          c%are_lclose(c%cen(:,2),(/0.5d0,0.5d0,0d0/),1d-5)
       ok = ok .and. &
          (c%are_lclose(c%cen(:,3),(/0d0,0.5d0,0.5d0/),1d-5) .or.&
          c%are_lclose(c%cen(:,3),(/0.5d0,0d0,0.5d0/),1d-5) .or.&
          c%are_lclose(c%cen(:,3),(/0.5d0,0.5d0,0d0/),1d-5))
       ok = ok .and. &
          (c%are_lclose(c%cen(:,4),(/0d0,0.5d0,0.5d0/),1d-5) .or.&
          c%are_lclose(c%cen(:,4),(/0.5d0,0d0,0.5d0/),1d-5) .or.&
          c%are_lclose(c%cen(:,4),(/0.5d0,0.5d0,0d0/),1d-5))
       if (ok) c%lcent = 6 ! F
    endif
    
  end subroutine set_lcent

