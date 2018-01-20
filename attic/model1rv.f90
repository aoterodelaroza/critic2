!In systemmod@proc.f90:

    elseif (trim(id) == "model1r") then
       ! xxxx
       iid = 0
       xp = syl%c%c2x(x0)
       call syl%c%nearest_atom(xp,iid,dist,lvec)
       if (dist < rc) then
          u = dist/rc
          field_feval%f = (1 - 10*u**3 + 15*u**4 - 6*u**5) * 21d0 / (5d0 * pi * rc**3)
       else
          field_feval%f = 0d0
       end if
       field_feval%fval = field_feval%f
    elseif (trim(id) == "model1v") then
       ! xxxx requires the +1 charges
       iid = 0
       xp = syl%c%c2x(x0)
       call syl%c%nearest_atom(xp,iid,dist,lvec)
       field_feval%f = syl%c%ewald_pot(xp,.false.) + syl%c%qsum * syl%c%eta**2 * pi / syl%c%omega 
       if (dist < rc .and. dist > 1d-6) then ! correlates with the value in the ewald routine
          u = dist/rc
          field_feval%f = field_feval%f + (12 - 14*u**2 + 28*u**5 - 30*u**6 + 9*u**7) / (5d0 * rc) - 1d0 / dist
       elseif (dist <= 1d-6) then
          u = dist/rc
          field_feval%f = field_feval%f + (12 - 14*u**2 + 28*u**5 - 30*u**6 + 9*u**7) / (5d0 * rc) - 2d0 / sqpi / syl%c%eta
       end if
       field_feval%fval = field_feval%f

!In arithmetic@proc.F90:

  !> Does this identifier correspond to a special field. This routine
  !> is thread-safe.
  module function isspecialfield(fid)
    use tools_io, only: lower
    character*(*), intent(in) :: fid
    logical :: isspecialfield

    isspecialfield = (trim(fid) == "ewald" .or. trim(fid) == "model1r" .or.&
       trim(fid) == "model1v") 
    
  end function isspecialfield

