! Copyright (c) 2016 Alberto Otero de la Roza <alberto@fluor.quimica.uniovi.es>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and Víctor Luaña
! <victor@fluor.quimica.uniovi.es>. 
!
! critic2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! critic2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Interface to DFTB+ wavefunctions.
module dftb_private
  implicit none
  
  private

  public :: dftb_read
  public :: dftb_register_struct
  public :: dftb_rho2

contains

  !> Read the information for a DFTB+ field from the detailed.xml, eigenvec.bin,
  !> and the basis set definition in HSD format.
  subroutine dftb_read(f,filexml,filebin,filehsd)
    use tools_io
    use types
    
    character*(*), intent(in) :: filexml !< The detailed.xml file
    character*(*), intent(in) :: filebin !< The eigenvec.bin file
    character*(*), intent(in) :: filehsd !< The definition of the basis set in hsd format
    type(field), intent(out) :: f !< Output field

    integer :: lu, i, j, k, idum, n
    character(len=:), allocatable :: line
    logical :: ok, iread(5)
    type(dftbatom) :: at

    ! detailed.xml, first pass
    lu = fopen_read(filexml)
    iread = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = adjustl(lower(line))
       if (index(line,"<real>") > 0) then
          f%isreal = next_logical(lu,line,"real")
          iread(1) = .true.
       elseif (index(line,"<nrofkpoints>") > 0) then
          f%nkpt = next_integer(lu,line,"nrofkpoints")
          iread(2) = .true.
       elseif (index(line,"<nrofspins>") > 0) then
          f%nspin = next_integer(lu,line,"nrofspins")
          iread(3) = .true.
       elseif (index(line,"<nrofstates>") > 0) then
          f%nstates = next_integer(lu,line,"nrofstates")
          iread(4) = .true.
       elseif (index(line,"<nroforbitals>") > 0) then
          f%norb = next_integer(lu,line,"nroforbitals")
          iread(5) = .true.
       end if
    end do
    if (.not.all(iread)) &
       call ferror('dftb_read','missing information in xml',faterr)

    ! detailed.xml, second pass
    allocate(f%dkpt(3,f%nkpt),f%dw(f%nkpt),f%docc(f%nstates,f%nkpt,f%nspin))
    rewind(lu)
    call read_kpointsandweights(lu,f%dkpt,f%dw)
    call read_occupations(lu,f%docc)
    call fclose(lu)

    ! use occupations times weights
    do i = 1, f%nkpt
       f%docc(:,i,:) = f%docc(:,i,:) * f%dw(i)
    end do

    ! eigenvec.bin
    lu = fopen_read(filebin,"unformatted")
    read (lu) idum ! this is the identity number

    ! read the eigenvectors
    if (f%isreal) then
       if (allocated(f%evecr)) deallocate(f%evecr)
       allocate(f%evecr(f%norb,f%norb,f%nspin))
       do i = 1, f%nspin
          do k = 1, f%norb
             read (lu) f%evecr(:,k,i)
          end do
       end do
    else
       if (allocated(f%evecc)) deallocate(f%evecc)
       allocate(f%evecc(f%norb,f%norb,f%nkpt,f%nspin))
       do i = 1, f%nspin
          do j = 1, f%nkpt
             do k = 1, f%norb
                read (lu) f%evecc(:,k,j,i)
             end do
          end do
       end do
    end if
    call fclose(lu)

    ! open the hsd file with the basis definition
    lu = fopen_read(filehsd)
    allocate(f%bas(10))
    n = 0
    do while(next_hsd_atom(lu,at))
       n = n + 1 
       if (n > size(f%bas)) call realloc(f%bas,2*n)
       f%bas = at
    end do
    call realloc(f%bas,n)
    call fclose(lu)

    ! finished
    f%init = .true.

  end subroutine dftb_read

  !> Calculate the density and derivatives of a DFTB+ field.
  subroutine dftb_rho2(f,xpos,nder,rho,grad,h)
    use tools_io
    use types
    use param

    type(field), intent(in) :: f !< Input field
    real*8, intent(in) :: xpos(3) !< Position in Cartesian
    integer, intent(in) :: nder  !< Number of derivatives
    real*8, intent(out) :: rho !< Density
    real*8, intent(out) :: grad(3) !< Gradient
    real*8, intent(out) :: h(3,3) !< Hessian 

    integer, parameter :: imax(0:2) = (/1,4,10/)
    
    rho = 0d0
    grad = 0d0
    h = 0d0

  end subroutine dftb_rho2

  !> Register structural information.
  subroutine dftb_register_struct()
    use types

    ! integer, intent(in) :: ncel
    ! type(celatom), intent(in) :: atcel(:)
    ! integer, intent(in) :: nneq
    ! type(atom), intent(in) :: at(:)

    ! integer :: i
    ! 
    ! nat = ncel
    ! if (allocated(xat)) deallocate(xat)
    ! allocate(xat(3,nat))
    ! do i = 1, nat
    !    xat(:,i) = atcel(i)%r
    ! end do
    
  end subroutine dftb_register_struct

  !xx! private !xx! 

  function next_logical(lu,line0,key0) result(next)
    use tools_io
    integer, intent(in) :: lu
    character*(*), intent(in) :: line0, key0
    logical :: next

    character(len=:), allocatable :: line, key
    integer :: idx, idxn, idxy
    logical :: found, ok, lastpass
    
    ! lowercase and build key
    key = "<" // trim(adjustl(key0)) // ">"
    line = lower(line0)

    ! delete up to the key in the current line, leave the rest
    idx = index(line,key)
    line = line(idx+len_trim(key):)
    line = trim(adjustl(line))

    ! parse all lines until </key> is found
    key = "</" // trim(adjustl(key0)) // ">"
    found = .false.
    lastpass = .false.
    next = .false.
    do while (.true.)
       ! are there any nos or yes? if so, which is first?
       idxn = index(line,"no")
       idxy = index(line,"yes")
       if (idxn > 0 .and. idxy == 0) then
          found = .true.
          next = .false.
       elseif (idxy > 0 .and. idxn == 0) then
          found = .true.
          next = .true.
       elseif (idxy > 0 .and. idxn > 0) then
          found = .true.
          next = (idxy < idxn)
       end if
       
       ! exit if found or if this was the last pass (after a </key> was read)
       if (found .or. lastpass) exit

       ! get a new line
       ok = getline_raw(lu,line)
       if (.not.ok) exit

       ! clean up and lowercase
       line = trim(adjustl(lower(line)))

       ! if the </key> is detected, read up to the </key> and activate the last pass flag
       idx = index(line,key)
       if (idx > 0) then
          line = line(1:idx-1)
          lastpass = .true.
       end if
    end do

    if (.not.found) &
       call ferror("next_logical","could not parse xml file",faterr)

  end function next_logical

  function next_integer(lu,line0,key0) result(next)
    use tools_io
    integer, intent(in) :: lu
    character*(*), intent(in) :: line0, key0
    integer :: next

    character(len=:), allocatable :: line, key
    integer :: idx, idxn, idxy, lp 
    logical :: found, ok, lastpass
    
    ! lowercase and build key
    key = "<" // trim(adjustl(key0)) // ">"
    line = lower(line0)

    ! delete up to the key in the current line, leave the rest
    idx = index(line,key)
    line = line(idx+len_trim(key):)
    line = trim(adjustl(line))
    lp = 1

    ! parse all lines until </key> is found
    key = "</" // trim(adjustl(key0)) // ">"
    found = .false.
    lastpass = .false.
    next = 0
    do while (.true.)
       ! are there any integers in this line? 
       found = isinteger(next,line,lp)
       
       ! exit if found or if this was the last pass (after a </key> was read)
       if (found .or. lastpass) exit

       ! get a new line
       ok = getline_raw(lu,line)
       if (.not.ok) exit

       ! clean up and lowercase
       line = trim(adjustl(lower(line)))

       ! if the </key> is detected, read up to the </key> and activate the last pass flag
       idx = index(line,key)
       if (idx > 0) then
          line = line(1:idx-1)
          lastpass = .true.
       end if
    end do

    if (.not.found) &
       call ferror("next_integer","could not parse xml file",faterr)

  end function next_integer

  subroutine read_kpointsandweights(lu,kpts,w)
    use tools_io
    integer, intent(in) :: lu
    real*8, intent(out) :: kpts(:,:)
    real*8, intent(out) :: w(:)

    logical :: ok, found
    character(len=:), allocatable :: line, key
    integer :: i

    key = "<kpointsandweights>"
    found = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = trim(adjustl(lower(line)))
       if (index(line,key) > 0) then
          found = .true.
          do i = 1, size(w)
             read (lu,*) kpts(:,i), w(i)
          end do
          exit
       end if
    end do

    if (.not.found) &
       call ferror("read_kpointsandweights","could not parse xml file",faterr)
    
  end subroutine read_kpointsandweights

  subroutine read_occupations(lu,occ)
    use tools_io
    integer, intent(in) :: lu
    real*8, intent(out) :: occ(:,:,:)

    logical :: ok, found
    character(len=:), allocatable :: line, key
    integer :: is, ik

    key = "<occupations>"
    found = .false.
    do while (.true.)
       ok = getline_raw(lu,line)
       if (.not.ok) exit
       line = trim(adjustl(lower(line)))
       if (index(line,key) > 0) then
          found = .true.

          do is = 1, size(occ,3)
             do ik = 1, size(occ,2)
                ! advance to the k-point with this name
                key = "<k" // string(ik) // ">"
                do while (.true.)
                   ok = getline_raw(lu,line,.true.)
                   if (index(line,key) > 0) exit
                end do
                ! read nstate real numbers
                occ(:,ik,is) = dftb_read_reals1(lu,size(occ,1))
             end do
          end do
          exit

       end if
    end do

    if (.not.found) &
       call ferror("read_occupations","could not parse xml file",faterr)
    
  end subroutine read_occupations

  !> Read a list of n reals from a logical unit
  function dftb_read_reals1(lu,n) result(x)
    use tools_io
    integer, intent(in) :: lu, n
    real*8 :: x(n)

    integer :: kk, lp
    real*8 :: rdum
    character(len=:), allocatable :: line
    logical :: ok

    kk = 0
    lp = 1
    ok = getline_raw(lu,line,.true.)
    line = trim(adjustl(line))
    do while(.true.)
       if (.not.isreal(rdum,line,lp)) then
          lp = 1
          ok = getline_raw(lu,line)
          line = trim(adjustl(line))
          if (.not.ok .or. line(1:1) == "<") exit
       else
          kk = kk + 1
          if (kk > n) call ferror("dftb_read_reals1","exceeded size of the array",2)
          x(kk) = rdum
       endif
    enddo

  endfunction dftb_read_reals1

  function next_hsd_atom(lu,at) result(ok)
    use tools_io
    use types
    integer, intent(in) :: lu
    type(dftbatom), intent(out) :: at
    logical :: ok

    character(len=:), allocatable :: line, word
    integer :: idx, nb, lp, i, n
    real*8 :: rdum

    ok = .false.
    at%norb = 0
    at%nexp = 0
    at%ncoef = 0
    nb = 0
    do while(getline(lu,line,.false.))
       if (len_trim(line) > 0) then
          ! must be an atom?
          idx = index(line,"{")
          at%name = line(1:idx-1)
          at%name = trim(adjustl(at%name))
          nb = nb + 1

          ! read the keywords at this level: atomicnumber and orbital
          do while(getline(lu,line,.true.))
             lp = 1
             word = lgetword(line,lp)
             idx = index(line,"{")
             if (idx > 0) word = word(1:idx-1)
             word = trim(adjustl(word))

             if (equal(word,"atomicnumber")) then
                idx = index(line,"=")
                word = line(idx+1:)
                read (word,*) at%z
             elseif (equal(word,"orbital")) then
                nb = nb + 1
                at%norb = at%norb + 1
                ! read the keywords at this level: 
                ! angularmomentum, occupation, cutoff, exponents, coefficients
                do while(getline(lu,line,.true.))
                   lp = 1
                   word = lgetword(line,lp)
                   idx = index(line,"{")
                   if (idx > 0) word = word(1:idx-1)
                   word = trim(adjustl(word))
                   if (equal(word,"angularmomentum")) then
                      idx = index(line,"=")
                      word = line(idx+1:)
                      read (word,*) at%l(at%norb)
                   elseif (equal(word,"occupation")) then
                      idx = index(line,"=")
                      word = line(idx+1:)
                      read (word,*) at%occ(at%norb)
                   elseif (equal(word,"cutoff")) then
                      idx = index(line,"=")
                      word = line(idx+1:)
                      read (word,*) at%cutoff(at%norb)
                   elseif (equal(word,"exponents")) then
                      nb = nb + 1
                      ok = getline(lu,line,.true.)
                      at%nexp(at%norb) = 0
                      lp = 1
                      do while (isreal(rdum,line,lp))
                         at%nexp(at%norb) = at%nexp(at%norb) + 1
                         at%eexp(at%nexp(at%norb),at%norb) = rdum
                      end do
                   elseif (equal(word,"coefficients")) then
                      nb = nb + 1
                      if (at%nexp(at%norb) == 0) &
                         call ferror('next_hsd_atom','coefficients must come after exponents',faterr)
                      do i = 1, at%nexp(at%norb)
                         ok = getline(lu,line,.true.)
                         n = 0
                         lp = 1
                         do while (isreal(rdum,line,lp))
                            n = n + 1
                            at%coef(n,i,at%norb) = rdum
                         end do
                         at%ncoef(i,at%norb) = n
                      end do
                   elseif (equal(word,"}")) then
                      nb = nb - 1
                      if (nb == 1) exit
                   else
                      call ferror("next_hsd_atom","unknown keyword (2): " // string(word),faterr)
                   end if
                end do
             elseif (equal(word,"}")) then
                nb = nb - 1
                if (nb == 0) exit
                exit
             else
                call ferror("next_hsd_atom","unknown keyword (1): " // string(word),faterr)
             end if
          end do

          ! succesfully read an atom
          ok = .true.
          exit
       end if
    end do

  end function next_hsd_atom

end module dftb_private
