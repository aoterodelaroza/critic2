! Copyright (c) 2007-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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

!> Scalar field class
submodule (fieldmod) proc
  implicit none

  !xx! private procedures
  ! function adaptive_stepper(fid,xpoint,h0,maxstep,eps,res)
  ! subroutine stepper_euler1(xpoint,grdt,h0,xout)
  ! subroutine stepper_heun(fid,xpoint,grdt,h0,xout,res)
  ! subroutine stepper_bs(fid,xpoint,grdt,h0,xout,xerr,res)
  ! subroutine stepper_rkck(fid,xpoint,grdt,h0,xout,xerr,res)
  ! subroutine stepper_dp(fid,xpoint,grdt,h0,xout,xerr,res)

  ! eps to move to the main cell
  real*8, parameter :: flooreps = 1d-4 ! border around unit cell

  ! numerical differentiation parameters
  real*8, parameter :: derw = 1.4d0, derw2 = derw*derw, big = 1d30, safe = 2d0

contains

  !> Adapt the size of an allocatable 1D type(field) array
  module subroutine realloc_field(a,nnew)
    use tools_io, only: ferror, faterr

    type(field), intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew

    type(field), allocatable :: temp(:)
    integer :: l1, u1

    if (.not.allocated(a)) &
       call ferror('realloc_field','array not allocated',faterr)
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))

    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc_field

  !> Deallocate and uninitialize
  module subroutine field_end(f)
    class(field), intent(inout) :: f

    nullify(f%c)
    f%isinit = .false.
    f%type = type_uninit
    f%usecore = .false.
    f%numerical = .false.
    f%exact = .false.
    f%typnuc = -3
    f%name = ""
    f%file = ""
    if (allocated(f%elk)) deallocate(f%elk)
    if (allocated(f%wien)) deallocate(f%wien)
    if (allocated(f%pi)) deallocate(f%pi)
    if (allocated(f%grid)) deallocate(f%grid)
    if (allocated(f%wfn)) deallocate(f%wfn)
    if (allocated(f%dftb)) deallocate(f%dftb)
    if (allocated(f%fr)) deallocate(f%fr)
    if (allocated(f%zpsp)) deallocate(f%zpsp)
    f%expr = ""
    f%fcp_deferred = .true.
    if (allocated(f%cp)) deallocate(f%cp)
    if (allocated(f%cpcel)) deallocate(f%cpcel)
    f%ncp = 0
    f%ncpcel = 0

  end subroutine field_end

  !> Sets the default options for the given field.
  module subroutine field_set_default_options(ff)
    class(field), intent(inout) :: ff

    if (allocated(ff%grid)) call ff%grid%setmode('default')
    if (allocated(ff%wien)) ff%wien%cnorm = .true.
    ff%exact = .false.
    ff%usecore = .false.
    ff%numerical = .false.
    ff%typnuc = -3
    if (.not.allocated(ff%zpsp)) allocate(ff%zpsp(ff%c%nspc))
    ff%zpsp = -1

  end subroutine field_set_default_options

  !> Set field flags. fid is the field slot. line is the input line to
  !> parse. oksyn is true in output if the syntax was OK. dormt is
  !> true in output if the MT test discontinuity test needs to be run
  !> for elk/wien2k fields.
  module subroutine field_set_options(ff,line,errmsg)
    use grid1mod, only: grid1_register_core
    use global, only: eval_next
    use tools_io, only: string, lgetword, equal, zatguess,&
       isinteger, getword, equali
    use param, only: sqfp
    class(field), intent(inout) :: ff
    character*(*), intent(in) :: line
    character(len=:), allocatable, intent(out) :: errmsg

    character(len=:), allocatable :: word, word2, aux
    integer :: lp, lp2, iz, iq, i, is
    logical :: ok
    real*8 :: norm

    errmsg = ""
    ! parse the rest of the line
    lp = 1
    do while (.true.)
       word = lgetword(line,lp)
       if (equal(word,'tricubic') .or. equal(word,'trispline') .or. &
          equal(word,'trilinear') .or. equal(word,'nearest') .or. &
          equal(word,'smoothrho')) then
          if (.not.allocated(ff%grid)) then
             errmsg = "tricubic/... incompatible with fields other than grids"
             return
          end if
          call ff%grid%setmode(word)

          do while (.true.)
             lp2 = lp
             word = lgetword(line,lp)
             if (equal(word,'nenv')) then
                ok = isinteger(ff%grid%smr_nenv,line,lp)
                if (.not.ok) then
                   errmsg = "wrong nenv option in smoothrho"
                   return
                end if
             elseif (equal(word,'fdmax')) then
                ok = eval_next(ff%grid%smr_fdmax,line,lp)
                if (.not.ok) then
                   errmsg = "wrong fdmax option in smoothrho"
                   return
                end if
             else
                lp = lp2
                exit
             end if
          end do

       else if (equal(word,'exact')) then
          ff%exact = .true.
       else if (equal(word,'approximate')) then
          ff%exact = .false.
       else if (equal(word,'rhonorm')) then
          if (.not.ff%type == type_wien.or..not.allocated(ff%wien)) then
             errmsg = "rhonorm incompatible with fields other than wien2k"
             return
          end if
          if (.not.ff%wien%cnorm) then
             ff%wien%cnorm = .true.
             if (allocated(ff%wien%slm)) &
                ff%wien%slm(:,1,:) = ff%wien%slm(:,1,:) / sqfp
          end if
       else if (equal(word,'vnorm')) then
          if (.not.ff%type == type_wien.or..not.allocated(ff%wien)) then
             errmsg = "vnorm incompatible with fields other than wien2k"
             return
          end if
          if (ff%wien%cnorm) then
             ff%wien%cnorm = .false.
             if (allocated(ff%wien%slm)) &
                ff%wien%slm(:,1,:) = ff%wien%slm(:,1,:) * sqfp
          end if
       else if (equal(word,'numerical')) then
          ff%numerical = .true.
       else if (equal(word,'analytical')) then
          ff%numerical = .false.
       else if (equal(word,'nocore')) then
          ff%usecore = .false.
          if (.not.allocated(ff%zpsp)) allocate(ff%zpsp(ff%c%nspc))
          ff%zpsp = -1
       else if (equal(word,'typnuc')) then
          ok = eval_next(ff%typnuc,line,lp)
          if (.not.ok) then
             errmsg = "wrong typnuc"
             return
          end if
          if (ff%typnuc /= -3 .and. ff%typnuc /= -1 .and. &
             ff%typnuc /= +1 .and. ff%typnuc /= +3) then
             errmsg = "wrong typnuc"
             return
          end if
       else if (equal(word,'normalize')) then
          if (.not.ff%type == type_grid) then
             errmsg = "normalize incompatible with fields other than grids"
             return
          end if
          ok = eval_next(norm,line,lp)
          if (.not. ok) then
             errmsg = "value for normalize keyword missing"
             return
          end if
          call ff%grid%normalize(norm,ff%c%omega)
       else if (equal(word,'zpsp')) then
          ff%usecore = .true.
          if (.not.allocated(ff%zpsp)) allocate(ff%zpsp(ff%c%nspc))
          do while (.true.)
             lp2 = lp

             ! try matching to species names
             word2 = getword(line,lp)
             is = -1
             do i = 1, ff%c%nspc
                if (equali(ff%c%spc(i)%name,word2)) then
                   is = i
                   exit
                end if
             end do
             if (is > 0) then
                aux = getword(line,lp)
                if (.not.isinteger(iq,aux)) then
                   errmsg = "wrong syntax in ZPSP"
                   return
                end if
                ff%zpsp(is) = iq
             else
                ! try matching the atomic number
                if (len_trim(word2) > 0 .and. len_trim(word2) <= 2) then
                   iz = zatguess(word2)
                   if (iz > 0) then
                      aux = getword(line,lp)
                      if (.not.isinteger(iq,aux)) then
                         errmsg = "wrong syntax in ZPSP"
                         return
                      end if
                      do i = 1, ff%c%nspc
                         if (ff%c%spc(i)%z == iz) then
                            ff%zpsp(i) = iq
                         end if
                      end do
                   else
                      ! must be some other keyword
                      lp = lp2
                      exit
                   end if
                else
                   ! must be some other keyword
                   lp = lp2
                   exit
                end if
             end if
          end do

          ! register the species
          ff%usecore = any(ff%zpsp > 0)
          if (ff%usecore) then
             do i = 1, ff%c%nspc
                if (ff%zpsp(i) > 0) &
                   call grid1_register_core(ff%c%spc(i)%z,ff%zpsp(i))
             end do
          end if
       else if (len_trim(word) > 0) then
          errmsg = "unknown extra keyword"
          return
       else
          exit
       end if
    end do

  end subroutine field_set_options

  !> Load a new field using the given field seed, the crystal
  !> structure pointer, the ID of the new field in the system (id)
  !> and the parent system's C pointer (sptr). If an error was
  !> found, returns a non-zero-length error message (errmsg).
  !> Some new fields are handled at system level (see load_field_string).
  module subroutine field_new(f,seed,c,id,sptr,errmsg,ti)
    use iso_c_binding, only: c_loc
    use types, only: realloc
    use fieldseedmod, only: fieldseed
    use arithmetic, only: eval, pretokenize, token
    use tools_io, only: equal, isinteger
    use param, only: ifformat_unknown, ifformat_wien, ifformat_elk, ifformat_pi,&
       ifformat_cube, ifformat_bincube, ifformat_abinit, ifformat_vasp,&
       ifformat_vaspnov, ifformat_qub,&
       ifformat_xsf, ifformat_elkgrid, ifformat_siestagrid, ifformat_fplogrid,&
       ifformat_dftb, ifformat_pwc, ifformat_fmt, ifformat_txt,&
       ifformat_wfn, ifformat_wfx, ifformat_fchk, ifformat_molden, ifformat_as,&
       ifformat_as_promolecular, ifformat_as_core, ifformat_as_resample,&
       ifformat_as_ft_x,ifformat_as_ft_y,ifformat_as_ft_z,ifformat_as_ft_xx,&
       ifformat_as_ft_xy,ifformat_as_ft_xz,ifformat_as_ft_yy,ifformat_as_ft_yz,&
       ifformat_as_ft_zz,ifformat_as_ft_grad,ifformat_as_ft_lap,ifformat_as_ft_pot,&
       ifformat_as_clm, ifformat_as_clm_sub, ifformat_as_ghost, &
       ifformat_copy, ifformat_promolecular, ifformat_promolecular_fragment
    use hashmod, only: hash
    use iso_c_binding, only: c_ptr
    use types, only: vstring
    class(field), intent(inout) :: f !< Input field
    type(fieldseed), intent(in) :: seed
    type(crystal), intent(inout), target :: c
    integer, intent(in) :: id
    type(c_ptr), intent(in) :: sptr
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: i, j, k, n(3)
    type(fragment) :: fr
    real*8 :: xdelta(3,3), x(3), rho
    real*8, allocatable :: faux(:,:,:)
    logical :: ifail
    type(vstring) :: lerrmsg
    type(token), allocatable :: toklist(:)

    errmsg = ""
    if (.not.c%isinit) then
       errmsg = "crystal not initialized"
       return
    end if
    call f%end()
    f%c => c
    f%sptr = sptr
    f%id = id
    f%name = adjustl(trim(seed%fid))

    ! set the initial pseudopotential charges
    if (.not.allocated(f%zpsp)) allocate(f%zpsp(f%c%nspc))
    f%zpsp = -1

    ! set the default field flags
    call f%set_default_options()

    ! interpret the seed and load the field
    if (seed%iff == ifformat_unknown) then
       errmsg = "unknown seed format"

    elseif (seed%iff == ifformat_wien) then
       if (.not.allocated(f%wien)) allocate(f%wien)
       call f%wien%end()
       call f%wien%read_clmsum(seed%file(1),seed%file(2),errmsg,ti=ti)
       f%type = type_wien
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_elk) then
       if (seed%nfile == 1) then
          if (.not.allocated(f%grid)) allocate(f%grid)
          call f%grid%end()
          call f%grid%read_elk(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
          f%type = type_grid
          f%file = seed%file(1)
       elseif (seed%nfile == 2) then
          if (.not.allocated(f%elk)) allocate(f%elk)
          call f%elk%end()
          call f%elk%read_out(c_loc(f%c),seed%file(1),seed%file(2),errmsg=errmsg,ti=ti)
          f%type = type_elk
          f%file = seed%file(1)
       else
          if (.not.allocated(f%elk)) allocate(f%elk)
          call f%elk%end()
          call f%elk%read_out(c_loc(f%c),seed%file(1),seed%file(2),seed%file(3),&
             errmsg=errmsg,ti=ti)
          f%type = type_elk
          f%file = seed%file(3)
       endif

    elseif (seed%iff == ifformat_pi) then
       if (.not.allocated(f%pi)) allocate(f%pi)
       call f%pi%read(c_loc(f%c),seed%nfile,seed%piat,seed%file,errmsg,ti=ti)
       f%type = type_pi
       f%file = "<pi ion files>"

    elseif (seed%iff == ifformat_cube) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_cube(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_bincube) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_bincube(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_abinit) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_abinit(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_vasp .or. seed%iff == ifformat_vaspnov) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_vasp(c_loc(c),seed%file(1),c%m_x2c,(seed%iff == ifformat_vasp),seed%vaspblk,&
          errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_qub) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_qub(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_xsf) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_xsf(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_fmt) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_fmt(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_txt) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_txt(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_elkgrid) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_elk(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_siestagrid) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_siesta(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_fplogrid) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       call f%grid%read_fplo(c_loc(c),seed%file(1),c%m_x2c,errmsg,ti=ti)
       f%type = type_grid
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_dftb) then
       if (.not.allocated(f%dftb)) allocate(f%dftb)
       call f%dftb%end()
       call f%dftb%read(c_loc(f%c),seed%file(1),seed%file(2),seed%file(3),errmsg,ti=ti)
       f%type = type_dftb
       f%file = seed%file(1)

    elseif (seed%iff == ifformat_pwc) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()

       f%type = type_grid
       f%file = seed%file(1)
       call f%grid%read_pwc(c_loc(c),seed%file(1),seed%pwcspin,seed%pwcikpt,seed%pwcibnd,&
          seed%pwcemin,seed%pwcemax,c%m_x2c,errmsg,ti=ti)
       if (seed%nfile == 2) then
          call f%grid%read_wannier_chk(seed%file(2),errmsg=errmsg,ti=ti)
       elseif (seed%nfile == 3) then
          call f%grid%read_wannier_chk(seed%file(2),seed%file(3),errmsg=errmsg,ti=ti)
       end if

    elseif (seed%iff == ifformat_wfn) then
       if (.not.allocated(f%wfn)) allocate(f%wfn)
       call f%wfn%end()
       call f%wfn%read_wfn(c_loc(f%c),seed%file(1),errmsg,ti=ti)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_wfx) then
       if (.not.allocated(f%wfn)) allocate(f%wfn)
       call f%wfn%end()
       call f%wfn%read_wfx(c_loc(f%c),seed%file(1),errmsg,ti=ti)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_fchk) then
       if (.not.allocated(f%wfn)) allocate(f%wfn)
       call f%wfn%end()
       call f%wfn%read_fchk(c_loc(f%c),seed%file(1),seed%readvirtual,errmsg,ti=ti)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_molden) then
       if (.not.allocated(f%wfn)) allocate(f%wfn)
       call f%wfn%end()
       call f%wfn%read_molden(c_loc(f%c),seed%file(1),seed%molden_type,seed%readvirtual,errmsg,ti=ti)
       f%type = type_wfn
       f%file = trim(seed%file(1))

    elseif (seed%iff == ifformat_promolecular) then
       call f%load_promolecular(f%c,id,"<promolecular>")
       f%name = "<promolecular>"

    elseif (seed%iff == ifformat_promolecular_fragment) then
       fr = f%c%identify_fragment_from_xyz(seed%file(1),errmsg,ti=ti)
       if (len_trim(errmsg) == 0) then
          if (fr%nat > 0) then
             call f%load_promolecular(f%c,id,trim(seed%file(1)),fr)
          else
             errmsg = "fragment contains unknown atoms"
          end if
       end if
       f%name = "<promolecular using fragment>"

    elseif (seed%iff == ifformat_as_promolecular.or.seed%iff == ifformat_as_core) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       call f%grid%end()
       if (seed%iff == ifformat_as_promolecular) then
          if (seed%nfile > 0) then
             fr = c%identify_fragment_from_xyz(seed%file(1),errmsg,ti=ti)
             if (len_trim(errmsg) == 0) then
                if (fr%nat > 0) then
                   call c%promolecular_array3(faux,seed%n,fr=fr)
                   call f%grid%from_array3(faux,c%m_x2c,c_loc(c))
                   deallocate(faux)
                else
                   errmsg = "zero atoms in the fragment"
                end if
             end if
             f%name = "<generated>, promolecular grid using fragment"
          else
             call c%promolecular_array3(faux,seed%n)
             call f%grid%from_array3(faux,c%m_x2c,c_loc(c))
             deallocate(faux)
             f%name = "<generated>, promolecular grid"
          end if
       else
          call f%set_options(seed%elseopt,errmsg) ! I need the zpsp now
          if (all(f%zpsp < 0)) then
             errmsg = "need zpsp values for load as core"
             return
          end if
          call c%promolecular_array3(faux,seed%n,zpsp=f%zpsp)
          call f%grid%from_array3(faux,c%m_x2c,c_loc(c))
          deallocate(faux)
          f%name = "<generated>, core grid"
       end if
       f%type = type_grid
       f%file = ""

    elseif (seed%iff == ifformat_as_ghost) then
       call f%load_ghost(c,id,"<ghost>",seed%expr,sptr)
       f%name = "<ghost>, " // trim(seed%expr)

    elseif (seed%iff == ifformat_as) then
       if (.not.allocated(f%grid)) allocate(f%grid)
       f%type = type_grid
       f%file = ""
       n = seed%n
       call f%grid%new_eval(sptr,c_loc(c),n,seed%expr,c%m_x2c)
       if (.not.f%grid%isinit) then
          ! This handles the cases with xc(.), structural variables,
          ! chemfunctions or fields.
          call f%grid%end()
          f%grid%n = n
          allocate(f%grid%f(n(1),n(2),n(3)))

          do i = 1, 3
             xdelta(:,i) = 0d0
             xdelta(i,i) = 1d0 / real(n(i),8)
          end do

          ! tokenizing takes a while so let's do it here and use the token list later
          call pretokenize(seed%expr,toklist,errmsg,sptr)

          ifail = .false.
          errmsg = ""
          !$omp parallel do private(x,rho,lerrmsg)
          do k = 1, n(3)
             do j = 1, n(2)
                do i = 1, n(1)
                   if (ifail) cycle
                   x = (i-1) * xdelta(:,1) + (j-1) * xdelta(:,2) + (k-1) * xdelta(:,3)
                   x = c%x2c(x)
                   rho = eval(seed%expr,lerrmsg%s,x,sptr,.true.,toklist)
                   !$omp critical(write)
                   if (len_trim(lerrmsg%s) > 0) then
                      errmsg = lerrmsg%s
                      ifail = .true.
                   end if
                   f%grid%f(i,j,k) = rho
                   !$omp end critical(write)
                end do
             end do
          end do
          !$omp end parallel do
          if (len_trim(errmsg) > 0) goto 999
          f%grid%isinit = .true.
       end if
       f%name = "<generated>, grid: " // trim(seed%expr)

    elseif (seed%iff == ifformat_copy .or. seed%iff == ifformat_as_ft_x .or.&
       seed%iff == ifformat_as_ft_y .or. seed%iff == ifformat_as_ft_z .or. &
       seed%iff == ifformat_as_ft_xx .or. seed%iff == ifformat_as_ft_xy .or. &
       seed%iff == ifformat_as_ft_xz .or. seed%iff == ifformat_as_ft_yy .or. &
       seed%iff == ifformat_as_ft_yz .or. seed%iff == ifformat_as_ft_zz .or. &
       seed%iff == ifformat_as_ft_grad .or. seed%iff == ifformat_as_ft_lap .or. &
       seed%iff == ifformat_as_ft_pot .or.&
       seed%iff == ifformat_as_clm .or. seed%iff == ifformat_as_clm_sub .or.&
       seed%iff == ifformat_as_resample) then
       errmsg = "error in file format for field_new"
    else
       errmsg = "unknown seed format"
    end if
    if (len_trim(errmsg) > 0) goto 999

    ! set the rest of the variables passed with the field
    call f%set_options(seed%elseopt,errmsg)
    if (len_trim(errmsg) > 0) goto 999

    ! if this is an "as core" field, deactivate the use of an additonal core (double count)
    if (seed%iff == ifformat_as_core) then ! same as nocore
       f%usecore = .false.
       if (.not.allocated(f%zpsp)) allocate(f%zpsp(f%c%nspc))
       f%zpsp = -1
    end if

    f%isinit = .true.
    call f%init_cplist()

    return
999 continue ! error termination with non-zero errmsg
    call f%end()

  end subroutine field_new

  !> Load a ghost field using expression expr. c = pointer to the
  !> crystal structure, id = numerical ID of the new field. name =
  !> name of the new field. sptr = C pointer to the parent system.
  module subroutine load_ghost(f,c,id,name,expr,sptr)
    use grid3mod, only: grid3
    use fragmentmod, only: fragment
    use hashmod, only: hash
    use iso_c_binding, only: c_ptr
    class(field), intent(inout) :: f !< Input/output field
    type(crystal), intent(in), target :: c
    integer, intent(in) :: id
    character*(*), intent(in) :: name
    character*(*), intent(in) :: expr
    type(c_ptr), intent(in) :: sptr

    if (.not.c%isinit) return
    call f%end()
    f%c => c
    f%id = id
    f%isinit = .true.
    f%type = type_ghost
    f%numerical = .true.
    f%exact = .false.
    f%name = adjustl(name)
    f%file = ""
    if (.not.allocated(f%zpsp)) allocate(f%zpsp(f%c%nspc))
    f%zpsp = -1
    f%expr = expr
    f%sptr = sptr
    call f%init_cplist()

  end subroutine load_ghost

  !> Load a promolecular density field using the given crystal
  !> structure. The ID and name of the field are also set using
  !> the provided arguments.
  module subroutine load_promolecular(f,c,id,name,fr)
    use fragmentmod, only: fragment
    class(field), intent(inout) :: f !< Input field
    type(crystal), intent(in), target :: c
    integer, intent(in) :: id
    character*(*), intent(in) :: name
    type(fragment), intent(in), optional :: fr

    if (.not.c%isinit) return
    call f%end()
    f%c => c
    f%id = id
    f%isinit = .true.
    if (present(fr)) then
       f%type = type_promol_frag
       f%fr = fr
    else
       f%type = type_promol
    end if
    f%numerical = .false.
    f%exact = .false.
    f%name = adjustl(name)
    f%file = ""
    f%typnuc = -3
    if (.not.allocated(f%zpsp)) allocate(f%zpsp(f%c%nspc))
    f%zpsp = -1
    call f%init_cplist()

  end subroutine load_promolecular

  !> Load field as a transformation of the 3d grid given in g. The
  !> transformation is given by ityp (ifformat_as_* flags in param).
  !> The resulting field is returned in f. The crystal structure (c)
  !> is used to inform the transformation. id and name are the
  !> numerical id and name of the new field. If ifformat_as_pot and
  !> isry_ is present and true, the output potential is in Rydberg. If
  !> ifformat_as_resample, the optional argument n contains the size
  !> of the resample grid, and it is mandatory.
  module subroutine load_as_fftgrid(f,c,id,name,g,ityp,isry_,n)
    use grid3mod, only: grid3
    use fragmentmod, only: fragment
    use tools_io, only: ferror, faterr
    use param, only: ifformat_as_resample
    class(field), intent(inout) :: f !< Input/output field
    type(crystal), intent(in), target :: c
    integer, intent(in) :: id
    character*(*), intent(in) :: name
    type(grid3), intent(in) :: g
    integer, intent(in) :: ityp
    logical, intent(in), optional :: isry_
    integer, intent(in), optional :: n(3)

    logical :: isry

    isry = .false.
    if (present(isry_)) isry = isry_

    if (ityp == ifformat_as_resample .and..not.present(n)) &
       call ferror("load_as_fftgrid","need size for fft resample",faterr)

    if (.not.c%isinit) return
    call f%end()
    f%c => c
    f%id = id
    f%isinit = .true.
    f%type = type_grid
    if (.not.allocated(f%grid)) allocate(f%grid)
    if (ityp == ifformat_as_resample) then
       call f%grid%resample(g,n)
    else
       call f%grid%fft(g,ityp)
    end if
    f%numerical = .false.
    f%exact = .false.
    f%name = adjustl(name)
    f%file = ""
    f%typnuc = -3
    if (.not.allocated(f%zpsp)) allocate(f%zpsp(f%c%nspc))
    f%zpsp = -1
    call f%init_cplist()

  end subroutine load_as_fftgrid

  !> Calculate the scalar field f at point v (Cartesian) and its
  !> derivatives up to nder. Return the results in res0. If periodic
  !> is present and false, consider the field is defined in a
  !> non-periodic system. This routine is thread-safe.
  recursive module subroutine grd(f,v,nder,res,fder,periodic)
    use wfn_private, only: wfn_uhf
    use arithmetic, only: eval
    use types, only: scalar_value
    use global, only: CP_hdegen
    use tools_math, only: rsindex
    use tools_io, only: ferror, faterr
    use param, only: icrd_cart
    class(field), intent(inout) :: f !< Input field
    real*8, intent(in) :: v(3) !< Target point in Cartesian coordinates
    integer, intent(in) :: nder !< Number of derivatives to calculate (or -1 for special field)
    type(scalar_value), intent(out) :: res !< Output density and related scalar properties
    character*(*), intent(in), optional :: fder !< modifier for the special field
    logical, intent(in), optional :: periodic !< Whether the system is to be considered periodic (molecules only)

    real*8 :: wx(3), wxr(3), wc(3), wcr(3), x(3)
    integer :: i, nid, idx(3)
    real*8 :: rho, grad(3), h(3,3), rhox(3), rhov(3), gradx(3,3), gradv(3,3)
    real*8 :: hx(3,3,3), hv(3,3,3), gkinx(3)
    real*8 :: fval(-ndif_jmax:ndif_jmax,3), fzero
    logical :: isgrid, per, skipvalassign
    character(len=:), allocatable :: errmsg

    real*8, parameter :: hini = 1d-3, errcnv = 1d-8
    real*8, parameter :: neargrideps = 1d-12

    if (.not.f%isinit) call ferror("grd","field not initialized",faterr)

    ! initialize output quantities
    call res%clear()
    res%avail_der1 = (nder > 0)
    res%avail_der2 = (nder > 1)

    ! check consistency
    if (nder < 0 .and. f%type /= type_wfn) &
       call ferror("grd","MO values can only be calculated for molecular wavefunctions",faterr)

    ! initialize flags
    if (present(periodic)) then
       per = periodic
    else
       per = .true.
    end if

    ! numerical derivatives
    res%isnuc = .false.
    if (f%numerical .and. nder >= 0) then
       wc = v
       fzero = grd0(f,v,periodic)
       res%f = fzero
       res%gf = 0d0
       res%hf = 0d0
       if (nder > 0) then
          ! x
          fval(:,1) = 0d0
          fval(0,1) = fzero
          res%gf(1) = f%der1i((/1d0,0d0,0d0/),v,hini,errcnv,fval(:,1),periodic)
          ! y
          fval(:,2) = 0d0
          fval(0,2) = fzero
          res%gf(2) = f%der1i((/0d0,1d0,0d0/),v,hini,errcnv,fval(:,2),periodic)
          ! z
          fval(:,3) = 0d0
          fval(0,3) = fzero
          res%gf(3) = f%der1i((/0d0,0d0,1d0/),v,hini,errcnv,fval(:,3),periodic)
          if (nder > 1) then
             ! xx, yy, zz
             res%hf(1,1) = f%der2ii((/1d0,0d0,0d0/),v,0.5d0*hini,errcnv,fval(:,1),periodic)
             res%hf(2,2) = f%der2ii((/0d0,1d0,0d0/),v,0.5d0*hini,errcnv,fval(:,2),periodic)
             res%hf(3,3) = f%der2ii((/0d0,0d0,1d0/),v,0.5d0*hini,errcnv,fval(:,3),periodic)
             ! xy, xz, yz
             res%hf(1,2) = f%der2ij((/1d0,0d0,0d0/),(/0d0,1d0,0d0/),v,hini,hini,errcnv,periodic)
             res%hf(1,3) = f%der2ij((/1d0,0d0,0d0/),(/0d0,0d0,1d0/),v,hini,hini,errcnv,periodic)
             res%hf(2,3) = f%der2ij((/0d0,1d0,0d0/),(/0d0,0d0,1d0/),v,hini,hini,errcnv,periodic)
             ! final
             res%hf(2,1) = res%hf(1,2)
             res%hf(3,1) = res%hf(1,3)
             res%hf(3,2) = res%hf(2,3)
          end if
       end if
       ! valence quantities
       res%fval = res%f
       res%gfmodval = res%gfmod
       res%del2fval = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)
       goto 999
    end if

    ! To the main cell. Add a small safe zone around the limits of the unit cell
    ! to prevent precision problems.
    wx = f%c%c2x(v)
    wxr = f%c%c2xr(v)
    if (per) then
       do i = 1, 3
          if (wx(i) < -flooreps .or. wx(i) > 1d0+flooreps) &
             wx(i) = wx(i) - real(floor(wx(i)),8)
          if (wxr(i) < -flooreps .or. wxr(i) > 1d0+flooreps) &
             wxr(i) = wxr(i) - real(floor(wxr(i)),8)
       end do
    else
       if (.not.f%c%ismolecule) &
          call ferror("grd","non-periodic calculation in a crystal",faterr)
       ! if outside the main cell and the field is limited to a certain region of
       ! space, nullify the result and exit
       if ((any(wx < -flooreps) .or. any(wx > 1d0+flooreps)) .and. &
          f%type == type_grid .or. f%type == type_wien .or. f%type == type_elk .or.&
          f%type == type_pi) then
          return
       end if
    end if
    wc = f%c%x2c(wx)
    wcr = f%c%xr2c(wxr)

    ! type selector
    skipvalassign = .false.
    select case(f%type)
    case(type_grid)
       isgrid = .false.
       if (nder == 0) then
          ! maybe we can get the grid point directly
          x = modulo(wx,1d0) * f%grid%n
          idx = nint(x)
          isgrid = all(abs(x-idx) < neargrideps)
       end if
       if (isgrid) then
          idx = modulo(idx,f%grid%n)+1
          res%f = f%grid%f(idx(1),idx(2),idx(3))
          res%gf = 0d0
          res%hf = 0d0
       else
          call f%grid%interp(wx,res%f,res%gf,res%hf)
       endif

    case(type_wien)
       call f%wien%rho2(wx,nder,res%f,res%gf,res%hf)
       res%gf = matmul(transpose(f%c%m_c2x),res%gf)
       res%hf = matmul(matmul(transpose(f%c%m_c2x),res%hf),f%c%m_c2x)

    case(type_elk)
       call f%elk%rho2(wx,nder,res%f,res%gf,res%hf)
       res%gf = matmul(transpose(f%c%m_c2x),res%gf)
       res%hf = matmul(matmul(transpose(f%c%m_c2x),res%hf),f%c%m_c2x)

    case(type_pi)
       call f%pi%rho2(wcr,f%exact,res%f,res%gf,res%hf)
       ! all work done in Cartesians in a finite environment.

    case(type_wfn)
       if (nder >= 0) then
          call f%wfn%rho2(wcr,nder,rhox,rhov,gradx,gradv,hx,hv,gkinx,res%vir,res%stress)
          skipvalassign = .true.
          ! arrange output values
          res%f = rhox(1)
          res%fval = rhov(1)
          res%fspin = rhox(2:3)
          res%fspinval = rhov(2:3)

          res%gf = gradx(:,1)
          res%gfmodval = norm2(gradv(:,1))
          res%gfmodspin(1) = norm2(gradx(:,2))
          res%gfmodspin(2) = norm2(gradx(:,3))
          res%gfmodspinval(1) = norm2(gradv(:,2))
          res%gfmodspinval(2) = norm2(gradv(:,3))

          res%gkin = gkinx(1)
          res%gkinspin = gkinx(2:3)

          res%hf = hx(:,:,1)
          res%del2fval = hv(1,1,1)+hv(2,2,1)+hv(3,3,1)
          res%lapspin(1) = hx(1,1,2) + hx(2,2,2) + hx(3,3,2)
          res%lapspin(2) = hx(1,1,3) + hx(2,2,3) + hx(3,3,3)
          res%lapspinval(1) = hv(1,1,2) + hv(2,2,2) + hv(3,3,2)
          res%lapspinval(2) = hv(1,1,3) + hv(2,2,3) + hv(3,3,3)

          res%avail_gkin = .true.
          res%avail_stress = .true.
          res%avail_vir = .true.
          res%avail_spin = .true.
          res%spinpol = (f%wfn%wfntyp == wfn_uhf)
       else
          call f%wfn%calculate_mo(wcr,res%fspc,fder)
          return
       end if
       ! all work done in Cartesians in a finite environment.

    case(type_dftb)
       call f%dftb%rho2(wcr,f%exact,nder,res%f,res%gf,res%hf,res%gkin)
       res%avail_gkin = .true.
       ! all work done in Cartesians in a finite environment.

    case(type_promol)
       call f%c%promolecular_atom(wcr,icrd_cart,res%f,res%gf,res%hf,nder)
       ! not needed because grd_atomic uses struct.

    case(type_promol_frag)
       call f%c%promolecular_atom(wcr,icrd_cart,res%f,res%gf,res%hf,nder,fr=f%fr)
       ! not needed because grd_atomic uses struct.

    case(type_ghost)
       res%f = eval(f%expr,errmsg,wc,f%sptr,periodic)
       res%gf = 0d0
       res%hf = 0d0
       if (len_trim(errmsg) > 0) &
          call ferror("grd",errmsg,faterr)

    case default
       call ferror("grd","unknown scalar field type",faterr)
    end select

    if (.not.skipvalassign) then
       ! save the valence-only value
       res%fval = res%f
       res%gfmodval = res%gfmod
       res%del2fval = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)
    end if

    ! augment with the core if applicable
    if (f%usecore .and. any(f%zpsp /= -1)) then
       call f%c%promolecular_atom(wc,icrd_cart,rho,grad,h,nder,zpsp=f%zpsp)
       res%f = res%f + rho
       res%gf  = res%gf + grad
       res%hf = res%hf + h
    end if

    ! if numerical, skip to here
999 continue

    ! If it's on a nucleus, nullify the gradient (may not be zero in
    ! grid fields, for instance)
    nid = f%c%identify_atom(wc,icrd_cart,distmax=1d-5)
    res%isnuc = (nid > 0)

    ! gradient
    if (nder > 0) then
       if (res%isnuc) res%gf = 0d0
       res%gfmod = norm2(res%gf)
    end if

    ! hessian
    if (nder > 1) then
       res%del2f = res%hf(1,1) + res%hf(2,2) + res%hf(3,3)
       res%hfevec = res%hf
       call rsindex(res%hfevec,res%hfeval,res%r,res%s,CP_hdegen)
    end if

  end subroutine grd

  !> Calculate only the value of the scalar field at the given point
  !> (v in Cartesian). If periodic is present and false, consider the
  !> field is defined in a non-periodic system. This routine is
  !> thread-safe.
  recursive module function grd0(f,v,periodic)
    use arithmetic, only: eval
    use tools_io, only: ferror, faterr
    use param, only: icrd_cart
    class(field), intent(inout) :: f
    real*8, dimension(3), intent(in) :: v !< Target point in cartesian or spherical coordinates.
    real*8 :: grd0
    logical, intent(in), optional :: periodic !< Whether the system is to be considered periodic (molecules only)

    real*8 :: wx(3), wxr(3), wc(3), wcr(3)
    integer :: i
    real*8 :: h(3,3), grad(3), rho, rhoaux, gkin, vir, stress(3,3)
    real*8 :: hx(3,3,3), hv(3,3,3), gradx(3,3), gradv(3,3)
    real*8 :: rhox(3), rhov(3), gkinx(3)
    logical :: per
    character(len=:), allocatable :: errmsg

    ! initialize
    if (present(periodic)) then
       per = periodic
    else
       per = .true.
    end if
    grd0 = 0d0

    ! To the main cell. Add a small safe zone around the limits of the unit cell
    ! to prevent precision problems.
    wx = f%c%c2x(v)
    wxr = f%c%c2xr(v)
    if (per) then
       do i = 1, 3
          if (wx(i) < -flooreps .or. wx(i) > 1d0+flooreps) &
             wx(i) = wx(i) - real(floor(wx(i)),8)
          if (wxr(i) < -flooreps .or. wxr(i) > 1d0+flooreps) &
             wxr(i) = wxr(i) - real(floor(wxr(i)),8)
       end do
    else
       if (.not.f%c%ismolecule) &
          call ferror("grd","non-periodic calculation in a crystal",faterr)
       ! if outside the main cell and the field is limited to a certain region of
       ! space, nullify the result and exit
       if ((any(wx < -flooreps) .or. any(wx > 1d0+flooreps)) .and. &
          f%type == type_grid .or. f%type == type_wien .or. f%type == type_elk .or.&
          f%type == type_pi) then
          return
       end if
    end if
    wc = f%c%x2c(wx)
    wcr = f%c%xr2c(wxr)

    ! type selector
    select case(f%type)
    case(type_grid)
       call f%grid%interp(wx,rho,grad,h)
    case(type_wien)
       call f%wien%rho2(wx,0,rho,grad,h)
    case(type_elk)
       call f%elk%rho2(wx,0,rho,grad,h)
    case(type_pi)
       call f%pi%rho2(wcr,f%exact,rho,grad,h)
    case(type_wfn)
       call f%wfn%rho2(wcr,0,rhox,rhov,gradx,gradv,hx,hv,gkinx,vir,stress)
       rho = rhox(1)
    case(type_dftb)
       call f%dftb%rho2(wcr,f%exact,0,rho,grad,h,gkin)
    case(type_promol)
       call f%c%promolecular_atom(wcr,icrd_cart,rho,grad,h,0)
    case(type_promol_frag)
       call f%c%promolecular_atom(wcr,icrd_cart,rho,grad,h,0,fr=f%fr)
    case(type_ghost)
       rho = eval(f%expr,errmsg,wc,f%sptr,periodic)
       if (len_trim(errmsg) > 0) &
          call ferror("grd",errmsg,faterr)
    case default
       call ferror("grd","unknown scalar field type",faterr)
    end select

    if (f%usecore .and. any(f%zpsp /= -1)) then
       call f%c%promolecular_atom(wc,icrd_cart,rhoaux,grad,h,0,zpsp=f%zpsp)
       rho = rho + rhoaux
    end if
    grd0 = rho

  end function grd0

  !> Function derivative using finite differences and Richardson's
  !> extrapolation formula. This routine is thread-safe.
  recursive module function der1i(f,dir,x,h,errcnv,pool,periodic)
    class(field), intent(inout) :: f
    real*8 :: der1i
    real*8, intent(in) :: dir(3)
    real*8, intent(in) :: x(3), h, errcnv
    real*8, intent(inout) :: pool(-ndif_jmax:ndif_jmax)
    logical, intent(in), optional :: periodic

    real*8 :: err
    real*8 :: ww, hh, erract, n(ndif_jmax,ndif_jmax)
    real*8 :: f0, fp, fm
    integer :: i, j
    integer :: nh

    der1i = 0d0
    nh = 0
    hh = h
    if (pool(nh+1) == 0d0) then
       f0 = f%grd0(x+dir*hh,periodic)
       pool(nh+1) = f0
    else
       f0 = pool(nh+1)
    end if
    fp = f0
    if (pool(nh+1) == 0d0) then
       f0 = f%grd0(x-dir*hh,periodic)
       pool(-nh-1) = f0
    else
       f0 = pool(-nh-1)
    end if
    fm = f0
    n(1,1) = (fp - fm) / (hh+hh)

    err = big
    do j = 2, ndif_jmax
       hh = hh / derw
       nh = nh + 1
       if (pool(nh+1) == 0d0) then
          f0 = f%grd0(x+dir*hh,periodic)
          pool(nh+1) = f0
       else
          f0 = pool(nh+1)
       end if
       fp = f0
       if (pool(-nh-1) == 0d0) then
          f0 = f%grd0(x-dir*hh,periodic)
          pool(-nh-1) = f0
       else
          f0 = pool(-nh-1)
       end if
       fm = f0
       n(1,j) = (fp - fm) / (hh+hh)
       ww = derw2
       do i = 2, j
          n(i,j) = (ww*n(i-1,j)-n(i-1,j-1))/(ww-1d0)
          ww = ww * derw2
          erract = max (abs(n(i,j)-n(i-1,j)), abs(n(i,j)-n(i-1,j-1)))
          if (erract.le.err) then
             err = erract
             der1i = n(i,j)
          endif
       enddo
       if ( abs(n(j,j)-n(j-1,j-1)) .gt. safe*err .or. err .le. errcnv ) return
    enddo

  end function der1i

  !> Function second derivative using finite differences and
  !> Richardson's extrapolation formula. This routine is thread-safe.
  recursive module function der2ii(f,dir,x,h,errcnv,pool,periodic)
    class(field), intent(inout) :: f
    real*8 :: der2ii
    real*8, intent(in) :: dir(3)
    real*8, intent(in) :: x(3), h, errcnv
    real*8, intent(inout) :: pool(-ndif_jmax:ndif_jmax)
    logical, intent(in), optional :: periodic

    real*8 :: err
    real*8 :: ww, hh, erract, n(ndif_jmax, ndif_jmax)
    real*8 :: fx, fp, fm, f0
    integer :: i, j
    integer :: nh

    der2ii = 0d0
    nh = 0
    hh = h
    if (pool(nh) == 0d0) then
       f0 = f%grd0(x,periodic)
       pool(nh) = f0
    else
       f0 = pool(nh)
    end if
    fx = 2 * f0
    if (pool(nh-1) == 0d0) then
       f0 = f%grd0(x-dir*(hh+hh),periodic)
       pool(nh-1) = f0
    else
       f0 = pool(nh-1)
    end if
    fm = f0
    if (pool(nh+1) == 0d0) then
       f0 = f%grd0(x+dir*(hh+hh),periodic)
       pool(nh+1) = f0
    else
       f0 = pool(nh+1)
    end if
    fp = f0
    n(1,1) = (fp - fx + fm) / (4d0*hh*hh)
    err = big
    do j = 2, ndif_jmax
       hh = hh / derw
       nh = nh + 1
       if (pool(nh+1) == 0d0) then
          f0 = f%grd0(x-dir*(hh+hh),periodic)
          pool(nh+1) = f0
       else
          f0 = pool(nh+1)
       end if
       fm = f0
       if (pool(-nh-1) == 0d0) then
          f0 = f%grd0(x+dir*(hh+hh),periodic)
          pool(-nh-1) = f0
       else
          f0 = pool(-nh-1)
       end if
       fp = f0
       n(1,j) = (fp - fx + fm) / (4d0*hh*hh)
       ww = derw2
       erract = 0d0
       do i = 2, j
          n(i,j) = (ww*n(i-1,j)-n(i-1,j-1))/(ww-1d0)
          ww = ww * derw2
          erract = max (abs(n(i,j)-n(i-1,j)), abs(n(i,j)-n(i-1,j-1)))
          if (erract.le.err) then
             err = erract
             der2ii = n(i,j)
          endif
       enddo
       if (erract .gt. safe*err .or. err .le. errcnv) return
    enddo
    return
  end function der2ii

  !> Function mixed second derivative using finite differences and
  !> Richardson's extrapolation formula. This routine is thread-safe.
  recursive module function der2ij(f,dir1,dir2,x,h1,h2,errcnv,periodic)
    class(field), intent(inout) :: f
    real*8 :: der2ij
    real*8, intent(in) :: dir1(3), dir2(3)
    real*8, intent(in) :: x(3), h1, h2, errcnv
    logical, intent(in), optional :: periodic

    real*8 :: err, fpp, fmp, fpm, fmm, f0
    real*8 :: hh1, hh2, erract, ww
    real*8 :: n(ndif_jmax,ndif_jmax)
    integer :: i, j

    der2ij = 0d0
    hh1 = h1
    hh2 = h2
    f0 = f%grd0(x+dir1*hh1+dir2*hh2,periodic)
    fpp = f0
    f0 = f%grd0(x+dir1*hh1-dir2*hh2,periodic)
    fpm = f0
    f0 = f%grd0(x-dir1*hh1+dir2*hh2,periodic)
    fmp = f0
    f0 = f%grd0(x-dir1*hh1-dir2*hh2,periodic)
    fmm = f0
    n(1,1) = (fpp - fmp - fpm + fmm ) / (4d0*hh1*hh2)
    err = big
    do j = 2, ndif_jmax
       hh1 = hh1 / derw
       hh2 = hh2 / derw
       f0 = f%grd0(x+dir1*hh1+dir2*hh2,periodic)
       fpp = f0
       f0 = f%grd0(x+dir1*hh1-dir2*hh2,periodic)
       fpm = f0
       f0 = f%grd0(x-dir1*hh1+dir2*hh2,periodic)
       fmp = f0
       f0 = f%grd0(x-dir1*hh1-dir2*hh2,periodic)
       fmm = f0
       n(1,j) = (fpp - fmp - fpm + fmm) / (4d0*hh1*hh2)
       ww = derw2
       do i = 2, j
          n(i,j) = (ww*n(i-1,j)-n(i-1,j-1))/(ww-1d0)
          ww = ww * derw2
          erract = max (abs(n(i,j)-n(i-1,j)), abs(n(i,j)-n(i-1,j-1)))
          if (erract.le.err) then
             err = erract
             der2ij = n(i,j)
          endif
       enddo
       if ( abs(n(j,j)-n(j-1,j-1)) .gt. safe*err .or. err .le. errcnv ) return
    enddo
  end function der2ij

  !> Return a string description of the field type.
  module function typestring(f,short) result(s)
    use tools_io, only: faterr, ferror
    class(field), intent(in) :: f
    character(len=:), allocatable :: s
    logical, intent(in) :: short

    if (.not.short) then
       select case (f%type)
       case (type_uninit)
          s = "not used"
       case (type_promol)
          s = "promolecular"
       case (type_grid)
          s = "grid"
       case (type_wien)
          s = "wien2k"
       case (type_elk)
          s = "elk"
       case (type_pi)
          s = "pi"
       case (type_wfn)
          s = "molecular wavefunction"
       case (type_dftb)
          s = "dftb+"
       case (type_promol_frag)
          s = "promolecular fragment"
       case (type_ghost)
          s = "ghost field"
       case default
          call ferror('typestring','unknown field type',faterr)
       end select
    else
       select case (f%type)
       case (type_uninit)
          s = "??"
       case (type_promol)
          s = "promol"
       case (type_grid)
          s = "grid"
       case (type_wien)
          s = "wien2k"
       case (type_elk)
          s = "elk"
       case (type_pi)
          s = "pi"
       case (type_wfn)
          s = "molwfn"
       case (type_dftb)
          s = "dftb+"
       case (type_promol_frag)
          s = "profrg"
       case (type_ghost)
          s = "ghost"
       case default
          call ferror('typestring','unknown field type',faterr)
       end select
    end if

  endfunction typestring

  !> Write information about the field to the standard output. If
  !> isload is true, show load-time information. If isset is true,
  !> show flags for this field.
  module subroutine printinfo(f,isload,isset)
    use grid3mod, only: mode_smr
    use wfn_private, only: molwfn, wfn_rhf, wfn_uhf, wfn_frac, &
       molden_type_orca, molden_type_psi4, molden_type_adf_sto
    use global, only: dunit0, iunit, iunitname0
    use tools_io, only: uout, string, ferror, faterr, ioj_center, nameguess, &
       ioj_right
    use param, only: hartoev
    class(field), intent(in) :: f
    logical, intent(in) :: isload
    logical, intent(in) :: isset
    real*8 :: fspin

    character(len=:), allocatable :: str, aux
    integer :: i, is, j, k, n(3)

    ! header
    if (.not.f%isinit) then
       write (uout,'("  Not initialized ")')
       return
    end if

    ! general information about the field
    if (len_trim(f%name) > 0) &
       write (uout,'("  Name: ",A)') string(f%name)
    if (len_trim(f%file) > 0) then
       write (uout,'("  Source: ",A)') string(f%file)
    else
       write (uout,'("  Source: <generated>")')
    end if
    write (uout,'("  Type: ",A)') f%typestring(.false.)

    ! type-specific
    if (f%type == type_promol) then
       !
    elseif (f%type == type_grid) then
       ! grids
       n = f%grid%n
       if (isload) then
          write (uout,'("  Grid dimensions : ",3(A,"  "))') (string(n(j)),j=1,3)
          write (uout,'("  Max. length Voronoi-relevant vector (",A,"):  ",2(A,"  "))') iunitname0(iunit), &
             string(f%grid%dmax*dunit0(iunit),'f',decimal=5)
          write (uout,'("  First elements... ",3(A,"  "))') (string(f%grid%f(1,1,j),'e',decimal=12),j=1,min(3,f%grid%n(3)))
          write (uout,'("  Last elements... ",3(A,"  "))') &
             (string(f%grid%f(n(1),n(2),n(3)-2+j),'e',decimal=12),j=3-min(3,f%grid%n(3)),2)
          write (uout,'("  Sum of elements... ",A)') string(sum(f%grid%f(:,:,:)),'e',decimal=12)
          write (uout,'("  Sum of squares of elements... ",A)') string(sum(f%grid%f(:,:,:)**2),'e',decimal=12)
          write (uout,'("  Cell integral (",A,") = ",A)') string(f%id),&
             string(sum(f%grid%f) * f%c%omega / real(product(n),8),'f',decimal=8)
          write (uout,'("  Min: ",A)') string(minval(f%grid%f),'e',decimal=8)
          write (uout,'("  Average: ",A)') string(sum(f%grid%f) / real(product(n),8),'e',decimal=8)
          write (uout,'("  Max: ",A)') string(maxval(f%grid%f),'e',decimal=8)
       end if
       if (isset) then
          write (uout,'("  Interpolation mode (1=nearest,2=linear,3=spline,4=tricubic,5=smoothrho): ",A)') string(f%grid%mode)
          if (f%grid%mode == mode_smr) then
             write (uout,'("    smoothrho with ",A," stencil nodes and smoothing dmax factor = ",A)') &
                string(f%grid%smr_nenv), string(f%grid%smr_fdmax,'f',decimal=4)
             write (uout,'("    Please cite: A. Otero-de-la-Roza, J. Chem. Phys. 156 (2022) 224116 (10.1063/5.0090232)")')
          end if
       end if
    elseif (f%type == type_wien) then
       if (isload) then
          write (uout,'("  Complex?: ",L)') f%wien%cmpl
          write (uout,'("  Spherical harmonics expansion LMmax: ",A)') string(size(f%wien%lm,2))
          write (uout,'("  Max. points in radial grid: ",A)') string(size(f%wien%slm,1))
          write (uout,'("  Total number of plane waves (new/orig): ",A,"/",A)') string(f%wien%nwav), string(f%wien%lastind)
       end if
       if (isset) then
          write (uout,'("  Density-style normalization? ",L)') f%wien%cnorm
       end if
    elseif (f%type == type_elk) then
       if (isload) then
          write (uout,'("  Number of LM pairs: ",A)') string(size(f%elk%rhomt,2))
          write (uout,'("  Max. points in radial grid: ",A)') string(size(f%elk%rhomt,1))
          write (uout,'("  Total number of plane waves: ",A)') string(size(f%elk%rhok))
       end if
    elseif (f%type == type_pi) then
       if (isset) then
          write (uout,'("  Exact calculation? ",L)') f%exact
       end if
       write (uout,'("  List of species and associated ion files")')
       write (uout,'("# spc  name    Z  ion file")')
       do i = 1, f%c%nspc
          if (f%pi%bas(i)%pi_used) then
             str = trim(f%pi%bas(i)%piname)
          else
             str = "<not used>"
          end if
          write (uout,'(99("  ",A))') &
             string(i,length=3,justify=ioj_right), &
             string(f%c%spc(f%c%at(i)%is)%name,length=5,justify=ioj_center), &
             string(f%c%spc(f%c%at(i)%is)%z,length=2,justify=ioj_right), &
             trim(str)
       end do

    elseif (f%type == type_wfn) then
       if (isload) then
          if (f%wfn%issto) then
             str = "(STO)"
          else
             str = "(GTO)"
          end if
          if (f%wfn%molden_type == molden_type_orca) then
             write (uout,'("  Molden file dialect: orca")')
          else if (f%wfn%molden_type == molden_type_psi4) then
             write (uout,'("  Molden file dialect: psi4")')
          else if (f%wfn%molden_type == molden_type_adf_sto) then
             write (uout,'("  Molden file dialect: ADF (STOs)")')
          end if
          if (f%wfn%wfntyp == wfn_rhf) then
             write (uout,'("  Wavefunction type: restricted")')
             write (uout,'("  Number of MOs (total): ",A)') string(f%wfn%nmoall)
             write (uout,'("  Number of MOs (occupied): ",A)') string(f%wfn%nmoocc)
          elseif (f%wfn%wfntyp == wfn_uhf) then
             write (uout,'("  Wavefunction type: unrestricted")')
             write (uout,'("  Number of MOs (total): ",A," (alpha=",A,",beta=",A,")")') &
                string(f%wfn%nmoall), string(f%wfn%nalpha+f%wfn%nalpha_virt), string(f%wfn%nmoall-(f%wfn%nalpha+f%wfn%nalpha_virt))
             write (uout,'("  Number of MOs (occupied): ",A," (alpha=",A,",beta=",A,")")') &
                string(f%wfn%nmoocc), string(f%wfn%nalpha), string(f%wfn%nmoocc-f%wfn%nalpha)
          elseif (f%wfn%wfntyp == wfn_frac) then
             write (uout,'("  Wavefunction type: fractional occupation")')
             write (uout,'("  Number of MOs: ",A)') string(f%wfn%nmoocc)
             write (uout,'("  Number of electrons: ",A)') string(nint(sum(f%wfn%occ(1:f%wfn%nmoocc))))
          end if
          write (uout,'("  Number of primitives ",A,": ",A)') str, string(f%wfn%npri)
          write (uout,'("  Number of EDFs: ",A)') string(f%wfn%nedf)
       end if
    elseif (f%type == type_dftb) then
       if (isload) then
          write (uout,'("  Number of states: ",A)') string(f%dftb%nstates)
          write (uout,'("  Number of spin channels: ",A)') string(f%dftb%nspin)
          write (uout,'("  Number of orbitals: ",A)') string(f%dftb%norb)
          write (uout,'("  Number of kpoints: ",A)') string(f%dftb%nkpt)
          write (uout,'("  Real wavefunction? ",L)') f%dftb%isreal
       end if
       if (isset) then
          write (uout,'("  Exact calculation? ",L)') f%exact
       end if
    elseif (f%type == type_promol_frag) then
       if (isload) then
          write (uout,'("  Number of atoms in fragment: ",A)') string(f%fr%nat)
       end if
    elseif (f%type == type_ghost) then
       write (uout,'("  Expression: ",A)') string(f%expr)
    else
       call ferror("printinfo","unknown field type",faterr)
    end if

    ! flags for any field
    if (isset) then
       write (uout,'("  Use core densities? ",L)') f%usecore
       if (any(f%zpsp > 0)) then
          str = ""
          do i = 1, f%c%nspc
             if (f%zpsp(i) > 0) then
                aux = str // string(f%c%spc(i)%name) // "(" // string(f%zpsp(i)) // "), "
                str = aux
             end if
          end do
          aux = str
          str = aux(1:len_trim(aux)-1)
          write (uout,'("  Core charges (ZPSP): ",A)') str
       end if
       write (uout,'("  Numerical derivatives? ",L)') f%numerical
       write (uout,'("  Nuclear CP signature: ",A)') string(f%typnuc)
    end if

    write (uout,'("  Number of non-equivalent critical points: ",A)') string(f%ncp)
    write (uout,'("  Number of critical points in the unit cell: ",A)') string(f%ncpcel)

    if (f%type == type_grid) then
       if (f%grid%isqe) then
          if (f%grid%qe%nspin == 1) then
             fspin = 2d0
          else
             fspin = 1d0
          end if

          write (uout,*)
          write (uout,'("+ Plane-wave Kohn-Sham states for this field")')
          write (uout,'("  Spin: ",A," K-points: ",A," Bands: ",A)') string(f%grid%qe%nspin), string(f%grid%qe%nks), &
             string(f%grid%qe%nbnd)
          do i = 1, f%grid%qe%nks
             write (uout,'("# kpt ",A," (",A," ",A," ",A,") w=",A)') string(i,2,justify=ioj_right), &
                (string(f%grid%qe%kpt(j,i),'f',decimal=4,justify=ioj_right),j=1,3),&
                string(f%grid%qe%wk(i),'e',decimal=8)
             do is = 1, f%grid%qe%nspin
                write (uout,'("  Ek(",A,"):",99(" ",A))') string(is),&
                   (string(f%grid%qe%ek(j,i,1) * hartoev,'f',6,2,justify=ioj_right),j=1,f%grid%qe%nbnd)
                write (uout,'(" Occ(",A,"):",99(" ",A))') string(is),&
                   (string(f%grid%qe%occ(j,i,1)/f%grid%qe%wk(i)*fspin,'f',6,2,justify=ioj_right),j=1,f%grid%qe%nbnd)
             end do
          end do
       end if

       if (f%grid%iswan) then
          write (uout,*)
          write (uout,'("+ Wannier functions information for this field")')
          write (uout,'("  Real-space lattice vectors: ",3(A," "))') (string(f%grid%qe%nk(i)),i=1,3)
          write (uout,'("  Spin: ",A," Bands: ",A)') string(f%grid%qe%nspin), string(f%grid%qe%nbnd)
          write (uout,'("  Wannier function centers (cryst. coords.) and spreads: ")')
          write (uout,'("# bnd spin        ----  center  ----        spread(",A,")")') iunitname0(iunit)
          do i = 1, f%grid%qe%nspin
             do j = 1, f%grid%qe%nbndw(i)
                write (uout,'("  ",99(A," "))') string(j,4,ioj_center), string(i,2,ioj_center), &
                   (string(f%grid%qe%center(k,j,i),'f',10,6,4),k=1,3),&
                   string(f%grid%qe%spread(j,i) * dunit0(iunit),'f',14,8,4)
             end do
          end do
       end if
    end if

  end subroutine printinfo

  !> Write the field info to a JSON object. The field object to json
  !> with root p.
  module subroutine write_json(f,json,p)
    use json_module, only: json_value, json_core
    class(field), intent(in) :: f
    type(json_core), intent(inout) :: json
    type(json_value), pointer, intent(inout) :: p

    type(json_value), pointer :: s

    if (.not.f%isinit) return
    call json%create_object(s,'field')
    call json%add(p,s)

    call json%add(s,'id',f%id)
    call json%add(s,'name',trim(f%name))
    if (len_trim(f%file) > 0) then
       call json%add(s,'source',trim(f%file))
    else
       call json%add(s,'source','generated')
    end if
    call json%add(s,'type',trim(f%typestring(.false.)))

  end subroutine write_json

  !> Initialize the critical point list with the atoms in the crystal
  !> structure.
  module subroutine init_cplist(f)
    use global, only: rbetadef, atomeps
    class(field), intent(inout) :: f

    integer :: i, j

    if (.not.f%c%isinit) return

    f%fcp_deferred = .true.
    if (allocated(f%cp)) deallocate(f%cp)
    if (allocated(f%cpcel)) deallocate(f%cpcel)

    !.initialize
    f%ncp = f%c%nneq
    allocate(f%cp(f%ncp))
    f%ncpcel = f%c%ncel
    allocate(f%cpcel(f%ncpcel))

    ! insert the nuclei in the list
    do i = 1, f%ncp
       f%cp(i)%x = f%c%at(i)%x - floor(f%c%at(i)%x)
       f%cp(i)%r = f%c%at(i)%r
       f%cp(i)%typ = f%typnuc
       f%cp(i)%typind = (f%cp(i)%typ+3)/2
       f%cp(i)%mult = f%c%at(i)%mult
       f%cp(i)%isdeg = .false.
       f%cp(i)%isnnm = .false.
       f%cp(i)%isnuc = .true.
       f%cp(i)%idx = i
       f%cp(i)%ir = 1
       f%cp(i)%ic = 1
       f%cp(i)%lvec = 0
       f%cp(i)%name = f%c%spc(f%c%at(i)%is)%name
       f%cp(i)%rbeta = rbetadef

       ! properties at the nuclei
       f%cp(i)%brdist = 0d0
       f%cp(i)%brpathlen = 0d0
       f%cp(i)%brang = 0d0

       ! calculate the point group symbol
       f%cp(i)%pg = f%c%sitesymm(f%cp(i)%x,atomeps)
    enddo

    ! add positions to the complete cp list
    do j = 1, f%c%ncel
       f%cpcel(j) = f%cp(f%c%atcel(j)%idx)
       f%cpcel(j)%idx = f%c%atcel(j)%idx
       f%cpcel(j)%ir = f%c%atcel(j)%ir
       f%cpcel(j)%ic = f%c%atcel(j)%ic
       f%cpcel(j)%lvec = f%c%atcel(j)%lvec

       f%cpcel(j)%x = matmul(f%c%rotm(1:3,1:3,f%cpcel(j)%ir),f%cp(f%cpcel(j)%idx)%x) + &
          f%c%rotm(:,4,f%cpcel(j)%ir) + f%c%cen(:,f%cpcel(j)%ic) + f%cpcel(j)%lvec
       f%cpcel(j)%r = f%c%x2c(f%cpcel(j)%x)
    end do

  end subroutine init_cplist

  !> Calculate the scalar field at the nuclei and fill the
  !> corresponding fields in f%cp and f%cpcel. This operation is
  !> deferred from init_cplist to avoid the slow down of field
  !> initalization it causes.
  module subroutine init_cplist_deferred(f)
    class(field), intent(inout) :: f

    integer :: i

    if (.not.f%fcp_deferred) return
    f%fcp_deferred = .false.

    do i = 1, f%c%nneq
       call f%grd(f%c%at(i)%r,2,f%cp(i)%s)
       f%cp(i)%s%r = 3
       f%cp(i)%s%s = f%typnuc
    end do

  end subroutine init_cplist_deferred

  !> Given the point xp in crystallographic coordinates, calculates
  !> the nearest CP. In output, nid is the id from the complete CP
  !> list and dist is the distance. If type is given, only search
  !> among the CPs of that type. If nid0, consider only CPs with index
  !> nid0 (complete list) If id0, consider only CPs with index id0
  !> (non-equivalent list). If nozero, skip zero-distance CPs.
  module subroutine nearest_cp(f,xp,nid,dist,lvec,type,nid0,id0,nozero)
    class(field), intent(in) :: f
    real*8, intent(in) :: xp(:)
    integer, intent(out) :: nid
    real*8, intent(out) :: dist
    integer, intent(out), optional :: lvec(3)
    integer, intent(in), optional :: type
    integer, intent(in), optional :: nid0
    integer, intent(in), optional :: id0
    logical, intent(in), optional :: nozero

    real*8, parameter :: eps = 1d-10

    real*8 :: temp(3), dd, dmin
    integer :: j

    nid = 0
    dmin = 1d30
    do j = 1, f%ncpcel
       if (present(type)) then
          if (f%cpcel(j)%typ /= type) cycle
       end if
       if (present(nid0)) then
          if (f%cpcel(j)%idx /= nid0) cycle
       end if
       if (present(id0)) then
          if (j /= id0) cycle
       end if

       if (f%c%ismolecule) then
          temp = f%cpcel(j)%x - (xp - floor(xp))
          temp = f%c%x2c(temp)
          dd = norm2(temp)
       else
          temp = f%cpcel(j)%x - xp
          call f%c%shortest(temp,dd)
       end if
       if (present(nozero)) then
          if (dd < eps) cycle
       end if
       if (dd < dmin) then
          nid = j
          dmin = dd
          if (present(lvec)) then
             lvec = nint(f%cpcel(j)%x - xp - temp)
          end if
       end if
    end do
    dist = dmin

  end subroutine nearest_cp

  !> Identify a CP in the unit cell. Input: position in cryst
  !> coords. Output: the ncpcel CP index, or 0 if none found. eps is
  !> the distance threshold.
  module function identify_cp(f,x0,eps)
    class(field), intent(in) :: f
    integer :: identify_cp
    real*8, intent(in) :: x0(3)
    real*8, intent(in) :: eps

    real*8 :: x(3), dist2
    integer :: i

    identify_cp = 0
    do i = 1, f%ncpcel
       x = x0 - f%cpcel(i)%x
       call f%c%shortest(x,dist2)
       if (dist2 < eps) then
          identify_cp = i
          return
       end if
    end do

  end function identify_cp

  !> Test the muffin tin discontinuity. ilvl = 0: quiet. ilvl = 1:
  !> normal output.  ilvl = 2: verbose output.
  module subroutine testrmt(f,ilvl,errmsg,ti)
    use tools_io, only: uout, ferror, warning, string, fopen_write, fclose
    use types, only: scalar_value
    use param, only: pi
    class(field), intent(inout) :: f
    integer, intent(in) :: ilvl
    character(len=:), allocatable, intent(out) :: errmsg
    type(thread_info), intent(in), optional :: ti

    integer :: n, i, j
    integer :: ntheta, nphi
    integer :: nt, np
    real*8 :: phi, theta, dir(3), xnuc(3), xp(3)
    real*8 :: fin, fout, gfin, gfout
    real*8 :: r, rmt
    character*4 :: label
    character(len=:), allocatable :: linefile
    integer :: luline, luplane
    integer  :: npass(f%c%nneq), nfail(f%c%nneq)
    logical :: ok
    real*8 :: epsm, epsp, mepsm, mepsp, dif, dosum, mindif, maxdif
    type(scalar_value) :: res

    real*8, parameter :: eps = 1d-3

    errmsg = ""
    if (f%type /= type_wien .and. f%type /= type_elk) then
       errmsg = "field must be of wien or elk type"
       return
    end if

    write (uout,'("* Muffin-tin discontinuity test")')

    ntheta = 10
    nphi = 10
    do n = 1, f%c%nneq
       if (f%type == type_wien) then
          rmt = f%wien%rmt_atom(f%c%at(n)%x)
       else
          rmt = f%elk%rmt(f%c%at(n)%is)
       end if
       mepsm = 0d0
       mepsp = 0d0
       if (ilvl > 1) then
          write (uout,'("+ Analysis of the muffin tin discontinuity for atom ",A)') string(n)
          write (uout,'("  ntheta = ",A)') string(ntheta)
          write (uout,'("  nphi = ",A)') string(nphi)
       end if

       xnuc = f%c%at(n)%x
       if (ilvl > 1) write (uout,'("  coords = ",3(A," "))') (string(xnuc(j),'f',decimal=9),j=1,3)
       xnuc = f%c%x2c(xnuc)

       if (ilvl > 1) then
          write (uout,'("  rmt = ",A)') string(rmt,'f',decimal=7)
          write (uout,'(2(A8," "),6(A12," "),A4)') "Azim.", "Polar", "f_in",&
             "f_out", "f_in-f_out", "gf_in", "gf_out", "gf_in-gf_out", "ok?"
          write (uout,'(100("-"))')
       end if
       npass(n) = 0
       nfail(n) = 0
       dosum = 0d0
       mindif = 1d10
       maxdif = -1d10
       if (ilvl > 1) then
          ! write line
          linefile = "plane_" // string(n,2,pad0=.true.) // ".dbg"
          luplane = fopen_write(linefile,ti=ti)
          write (luplane,'("#",A,I3)') " atom: ", n
          write (luplane,'("#",A,1p,3(E20.13," "))') " at: ", f%c%at(n)%x
          write (luplane,'("#",A,1p,3(E20.13," "))') " atc: ", xnuc
          write (luplane,'("#",A,1p,E20.13)') " rmt: ", rmt
          write (luplane,'("#  theta phi in out")')
       end if
       do nt = 1, ntheta
          do np = 0, nphi
             if ((np == 0 .or. np == nphi) .and. nt /= 1) cycle
             phi = real(np,8) * pi / nphi
             theta = real(nt,8) * 2d0 * pi / ntheta
             dir(1) = 1d0 * cos(theta) * sin(phi)
             dir(2) = 1d0 * sin(theta) * sin(phi)
             dir(3) = 1d0 * cos(phi)

             xp = xnuc + (rmt+eps) * dir
             call f%grd(xp,1,res)
             fout = res%f
             gfout = dot_product(res%gf,xp-xnuc) / (rmt+eps)
             xp = xnuc + (rmt-eps) * dir
             call f%grd(xp,1,res)
             fin = res%f
             gfin = dot_product(res%gf,xp-xnuc) / (rmt-eps)

             dif = fout - fin
             dosum = dosum + dif * dif
             if (dif < mindif) mindif = dif
             if (dif > maxdif) maxdif = dif

             if (ilvl > 1) then
                ! write line
                write (luplane,'(1p,4(E20.13," "))') theta, phi, fin, fout
             end if

             if (gfin*gfout > 0d0) then
                label = "pass"
                npass(n) = npass(n) + 1
             else
                label = "fail"
                if (ilvl > 2) then
                   ! write line
                   linefile = "line_" // string(n,2,pad0=.true.) // "_" // string(nt,3,pad0=.true.) //&
                      "_" // string(np,3,pad0=.true.) // ".dbg"
                   luline = fopen_write(linefile,ti=ti)
                   write (luline,'("#",A,I3)') " atom: ", n
                   write (luline,'("#",A,1p,3(E20.13," "))') " at: ", f%c%at(n)%x
                   write (luline,'("#",A,1p,3(E20.13," "))') " atc: ", xnuc
                   write (luline,'("#",A,1p,E20.13)') " rmt: ", rmt
                   write (luline,'("#",A,1p,3(E20.13," "))') " dir: ", dir
                   write (luline,'("#",A,1p,E20.13)') " r_ini: ", 0.50d0 * rmt
                   write (luline,'("#",A,1p,E20.13)') " r_end: ", 4.50d0 * rmt
                   do i = 0, 1000
                      r = 0.50d0 * rmt + (real(i,8) / 1000) * 4d0 * rmt
                      xp = xnuc + r * dir
                      call f%grd(xp,1,res)
                      write (luline,'(1p,3(E20.13," "))') r, res%f, dot_product(res%gf,xp-xnuc) / r
                   end do
                   call fclose(luline)
                end if
                epsm = 0d0
                epsp = 0d0
                mepsm = min(epsm,mepsm)
                mepsp = max(epsp,mepsp)
                nfail(n) = nfail(n) + 1
             end if
             if (ilvl > 1) write (uout,'(2(F8.4," "),1p,6(E12.4," "),0p,A4)') &
                theta, phi, fin, fout, fin-fout, gfin, gfout, gfin-gfout, label
          end do
       end do
       dosum = sqrt(dosum / (npass(n)+nfail(n)))
       if (ilvl > 1) then
          write (uout,'(100("-"))')
          write (uout,*)
       end if
       if (ilvl > 1) then
          call fclose(luplane)
       end if
       if (nfail(n) > 0) then
          write (uout,'("  Atom: ",A," delta_m = ",A," delta_p = ",A)') &
             string(n), string(mepsm+1d-3,'f',decimal=6), string(mepsp+1d-3,'f',decimal=6)
          write (uout,*)
       end if
       write (uout,'("  Atom: ",A," rmt= ",A," RMS/max/min(fout-fin) = ",3(A,"  "))') &
          string(n), string(rmt,'f',decimal=7), string(dosum,'f',decimal=6), &
          string(maxdif,'f',decimal=6), string(mindif,'f',decimal=6)
    end do

    ok = .true.
    if (ilvl > 0) then
       write (uout,'("+ Summary ")')
       write (uout,'(A4,3(" ",A7))') "Atom", "Pass", "Fail", "Total"
    end if
    do n = 1, f%c%nneq
       if (ilvl > 0) write (uout,'(I4,3(" ",I7))') n, npass(n), nfail(n), npass(n)+nfail(n)
       ok = ok .and. (nfail(n) == 0)
    end do
    if (ilvl > 0) write (uout,*)
    write (uout,'("+ Assert - no spurious CPs on the muffin tin surface: ",L1/)') ok
    if (.not.ok) call ferror('testrmt','Spurious CPs on the muffin tin surface!',warning)

  end subroutine testrmt

  !> Test the speed of the grd call by calculating npts random point
  !> in the unit cell. Write the results to the standard output.
  module subroutine benchmark(f,npts)
    use tools_io, only: uout, string
    use types, only: scalar_value
    implicit none
    class(field), intent(inout) :: f
    integer, intent(in) :: npts

    integer :: i, j
    real*8 :: x(3), aux(3), dist
    integer :: c1, c2, rate
    real*8, allocatable :: randn(:,:,:)
    type(scalar_value) :: res
    logical :: inrmt
    integer, allocatable :: wpts(:)

    allocate(wpts(f%c%nneq+1))
    write (uout,'("* Benchmark of the field ")')
    write (uout,'("* Field : ",A)') string(f%id)

    if (f%type == type_wien .or. f%type == type_elk) then
       wpts = 0
       allocate(randn(3,npts,f%c%nneq+1))

       out: do while (any(wpts(1:f%c%nneq+1) < npts))
          call random_number(x)
          do i = 1, f%c%ncel
             if (f%c%atcel(i)%idx > f%c%nneq) cycle
             aux = x - f%c%atcel(i)%x
             aux = f%c%x2c(aux - nint(aux))
             dist = norm2(aux)
             if (f%type == type_wien) then
                inrmt = (dist < f%wien%rmt_atom(f%c%at(f%c%atcel(i)%idx)%x))
             else
                inrmt = (dist < f%elk%rmt(f%c%atcel(i)%is))
             end if
             if (inrmt) then
                if (wpts(f%c%atcel(i)%idx) >= npts) cycle out
                wpts(f%c%atcel(i)%idx) = wpts(f%c%atcel(i)%idx) + 1
                randn(:,wpts(f%c%atcel(i)%idx),f%c%atcel(i)%idx) = f%c%x2c(x)
                cycle out
             end if
          end do
          if (wpts(f%c%nneq+1) >= npts) cycle out
          wpts(f%c%nneq+1) = wpts(f%c%nneq+1) + 1
          randn(:,wpts(f%c%nneq+1),f%c%nneq+1) = f%c%x2c(x)
       end do out

       write (uout,'("* Benchmark of muffin / interstitial grd ")')
       write (uout,'("* Number of points per zone : ",A)') string(npts)
       do i = 1, f%c%nneq+1
          call system_clock(count=c1,count_rate=rate)
          do j = 1, npts
             call f%grd(randn(:,j,i),0,res)
          end do
          call system_clock(count=c2)
          if (i <= f%c%nneq) then
             write (uout,'("* Atom : ",A)') string(i)
          else
             write (uout,'("* Interstitial ")')
          end if
          write (uout,'("* Total wall time : ",A)') trim(string(real(c2-c1,8) / rate,'f',12,6)) // " s"
          write (uout,'("* Avg. wall time per call : ",A)') &
             trim(string(real(c2-c1,8) / rate / real(npts,8) * 1d6,'f',12,4)) // " us"
       end do
       write (uout,*)
       deallocate(randn)

    else
       allocate(randn(3,npts,1))

       call random_number(randn)
       do i = 1, npts
          randn(1:3,i,1) = f%c%x2c(randn(1:3,i,1))
       end do

       ! grd
       write (uout,'("  Benchmark of the grd call ")')
       write (uout,'("  Number of points : ",A)') string(npts)
       call system_clock(count=c1,count_rate=rate)
       do i = 1, npts
          call f%grd(randn(:,i,1),0,res)
       end do
       call system_clock(count=c2)
       write (uout,'("  Total wall time : ",A)') trim(string(real(c2-c1,8) / rate,'f',12,6)) // " s"
       write (uout,'("  Avg. wall time per call : ",A)') &
          trim(string(real(c2-c1,8) / rate / real(npts,8) * 1d6,'f',12,4)) // " us"
       write (uout,*)

       deallocate(randn)

    end if
    deallocate(wpts)

  end subroutine benchmark

  !> Do a Newton-Raphson search at point r (Cartesian). A CP is found
  !> when the gradient is less than gfnormeps. ier is the exit code: 0
  !> (success), 1 (singular Hessian or inverse fail), and 2 (too many iterations).
  module subroutine newton(f,r,gfnormeps,ier)
    use tools_math, only: matinvsym
    use types, only: scalar_value
    class(field), intent(inout) :: f
    real*8, dimension(3), intent(inout) :: r
    integer, intent(out) :: ier
    real*8, intent(in) :: gfnormeps

    integer :: it
    real*8 :: wx(3)
    type(scalar_value) :: res

    integer, parameter :: maxit = 200

    do it = 1, maxit
       ! evaluate and stop criterion
       call f%grd(r,2,res)
       wx = f%c%c2x(r)
       ! write (*,'("xx ",I4," ",3(F14.6," "),1p,3(E14.6," "))') it, wx, res%gfmod

       if (res%gfmod < gfnormeps) then
          ier = 0
          return
       end if

       ! invert the Hessian
       call matinvsym(res%hf,3,ier)
       if (ier /= 0) return

       ! advance
       r = r - matmul(res%hf,res%gf)
    end do

    ! too many iterations
    ier = 2

  end subroutine newton

  !> Try to add the candidate CP x0 to the CP list. This routine
  !> checks if it is already known, and calculates all the relevant
  !> information: multiplicity, type, shell, assoc. nucleus, etc.
  !> Also, updates the complete CP list. Input x0 in Cartesian. If
  !> itype is present, assign itype as the type of critical bond
  !> (-3,-1,1,3). If typeok is present, only accept a CP type if the
  !> corresponding element in typeok is true.
  module subroutine addcp(f,x0,cpeps,nuceps,nucepsh,gfnormeps,itype,typeok)
    use arithmetic, only: eval
    use tools_io, only: ferror, string
    use types, only: scalar_value, realloc
    use global, only: rbetadef
    use types, only: realloc
    use param, only: icrd_crys
    class(field), intent(inout) :: f
    real*8, intent(in) :: x0(3) !< Position of the CP, in Cartesian coordinates
    real*8, intent(in) :: cpeps !< Discard CPs closer than cpeps from other CPs
    real*8, intent(in) :: nuceps !< Discard CPs closer than nuceps from nuclei
    real*8, intent(in) :: nucepsh !< Discard CPs closer than nucepsh from hydrogen
    real*8, intent(in) :: gfnormeps !< Discard CPs with gradient norm higher than this value
    integer, intent(in), optional :: itype !< Force a CP type (useful in grids)
    logical, intent(in), optional :: typeok(4) !< Filter out by CP type

    real*8 :: xc(3), x1(3)
    integer :: nid
    real*8 :: dist
    integer :: n, i, num, idx
    real*8, allocatable  :: sympos(:,:)
    integer, allocatable :: symrotm(:), symcenv(:)
    integer :: lnuc, lshell
    character*(1), parameter :: smallnamecrit(0:3) =(/'n','b','r','c'/)
    type(scalar_value) :: res

    ! Transform to cryst. and main cell
    xc = f%c%c2x(x0)
    xc = xc - floor(xc)

    ! If the scalar field has shells, then determine associated nucleus (lnuc) and shell (lshell).
    ! valence / core associated shell
    lnuc = 0
    lshell = 0

    ! If it's a molecule and the CP is in the border, return
    if (f%c%ismolecule) then
       if (xc(1) < f%c%molborder(1) .or. xc(1) > (1d0-f%c%molborder(1)) .or.&
          xc(2) < f%c%molborder(2) .or. xc(2) > (1d0-f%c%molborder(2)) .or.&
          xc(3) < f%c%molborder(3) .or. xc(3) > (1d0-f%c%molborder(3))) goto 999
    endif

    ! distance to other CPs
    call f%nearest_cp(xc,nid,dist)
    if (nid > f%c%ncel .and. dist < cpeps) then
       goto 999
    end if

    ! distance to atoms
    nid = f%c%identify_atom(xc,icrd_crys,dist=dist,distmax=max(nuceps,nucepsh))
    if (nid > 0) then
       if (f%c%spc(f%c%atcel(nid)%is)%z == 1) then
          if (dist < nucepsh) goto 999
       else
          if (dist < nuceps) goto 999
       end if
    end if

    ! density info
    ! verify that the gradient norm is higher than the eps
    ! apply type filter
    x1 = f%c%x2c(xc)
    call f%grd(x1,2,res)
    if (res%gfmod >= gfnormeps) goto 999
    if (present(typeok)) then
       idx = (res%s+5)/2
       if (idx >= 1 .and. idx <= 4) then
          if (.not.typeok(idx)) goto 999
       end if
    end if

    ! reallocate if more slots are needed for the new cp
    if (f%ncp >= size(f%cp)) then
       call realloc(f%cp,2*f%ncp)
    end if

    ! Write critical point in the list
    f%ncp = f%ncp + 1
    n = f%ncp
    f%cp(n)%x = xc
    f%cp(n)%r = x1

    ! Density info
    f%cp(n)%s = res

    ! Symmetry info
    f%cp(n)%pg = f%c%sitesymm(f%cp(n)%x,cpeps)
    f%cp(n)%idx = n
    f%cp(n)%ir = 1
    f%cp(n)%ic = 1
    f%cp(n)%lvec = 0

    ! Type of critical point
    if (.not.present(itype)) then
       f%cp(n)%isdeg = (f%cp(n)%s%r /= 3)
       f%cp(n)%typ = f%cp(n)%s%s
       f%cp(n)%typind = (f%cp(n)%s%s+3)/2
       f%cp(n)%isnuc = .false.
       f%cp(n)%isnnm = (f%cp(n)%s%s == f%typnuc)
    else
       f%cp(n)%s%r = 3
       f%cp(n)%s%s = itype
       f%cp(n)%isdeg = .false.
       f%cp(n)%typ = itype
       f%cp(n)%typind = (itype+3)/2
       f%cp(n)%isnuc = .false.
       f%cp(n)%isnnm = (itype == f%typnuc)
    end if

    ! discard degenerate critical points
    if (f%cp(n)%isdeg) then
       f%ncp = f%ncp - 1
       goto 999
    end if

    ! Wait until reconstruction to calculate graph properties
    f%cp(n)%brdist = 0d0
    f%cp(n)%brpathlen = 0d0
    f%cp(n)%brang = 0d0
    f%cp(n)%ipath = 0
    f%cp(n)%rbeta = Rbetadef
    f%cp(n)%brvec = 0d0

    ! Name
    num = count(f%cp(1:n)%typ == f%cp(n)%s%s)
    f%cp(n)%name = smallnamecrit(f%cp(n)%typind) // string(num)

    ! Add positions to the complete CP list
    call f%c%symeqv(f%cp(n)%x,f%cp(n)%mult,sympos,symrotm,symcenv,cpeps)
    do i = 1, f%cp(n)%mult
       f%ncpcel = f%ncpcel + 1
       if (f%ncpcel >= size(f%cpcel)) then
          call realloc(f%cpcel,2*f%ncpcel)
       end if

       f%cpcel(f%ncpcel) = f%cp(n)
       f%cpcel(f%ncpcel)%x = sympos(:,i)
       f%cpcel(f%ncpcel)%r = f%c%x2c(f%cpcel(f%ncpcel)%x)
       f%cpcel(f%ncpcel)%idx = n
       f%cpcel(f%ncpcel)%ir = symrotm(i)
       f%cpcel(f%ncpcel)%ic = symcenv(i)
       f%cpcel(f%ncpcel)%lvec = nint(f%cpcel(f%ncpcel)%x - &
          (matmul(f%c%rotm(1:3,1:3,symrotm(i)),f%cp(n)%x) + &
          f%c%rotm(1:3,4,symrotm(i)) + f%c%cen(:,symcenv(i))))
    end do
    deallocate(sympos,symrotm,symcenv)
999 continue

  end subroutine addcp

  ! Sort the non-equivalent cp list
  module subroutine sortcps(f,cpeps)
    use tools, only: mergesort
    use types, only: cp_type
    class(field), intent(inout) :: f
    real*8, intent(in) :: cpeps !< Discard CPs closer than cpeps from other CPs

    integer :: i, j, mi
    real*8, allocatable :: sympos(:,:)
    integer, allocatable :: symrotm(:), symcenv(:), iaux(:)
    real*8, allocatable :: raux(:)
    integer, allocatable :: iperm(:)

    if (f%fcp_deferred) then
       call f%init_cplist_deferred()
       f%fcp_deferred = .false.
    end if

    ! Sort the nneq CP list: first by type then by density, the atoms at the top
    allocate(iperm(f%ncp),raux(f%ncp),iaux(f%ncp))
    do i = 1, f%ncp
       raux(i) = -f%cp(i)%s%f
       iaux(i) = f%cp(i)%typind
       iperm(i) = i
    end do
    call mergesort(raux,iperm,f%c%nneq+1,f%ncp)
    call mergesort(iaux,iperm,f%c%nneq+1,f%ncp)
    f%cp = f%cp(iperm)

    ! Rewrite the complete CP list
    f%ncpcel = 0
    do i = 1, f%ncp
       call f%c%symeqv(f%cp(i)%x,mi,sympos,symrotm,symcenv,cpeps)
       do j = 1, mi
          f%ncpcel = f%ncpcel + 1
          f%cpcel(f%ncpcel) = f%cp(i)
          f%cpcel(f%ncpcel)%x = sympos(:,j)
          f%cpcel(f%ncpcel)%r = f%c%x2c(f%cpcel(f%ncpcel)%x)
          f%cpcel(f%ncpcel)%idx = i
          f%cpcel(f%ncpcel)%ir = symrotm(j)
          f%cpcel(f%ncpcel)%ic = symcenv(j)
          f%cpcel(f%ncpcel)%lvec = nint(f%cpcel(f%ncpcel)%x - &
             (matmul(f%c%rotm(1:3,1:3,symrotm(j)),f%cp(i)%x) + &
             f%c%rotm(1:3,4,symrotm(j)) + f%c%cen(:,symcenv(j))))
       end do
    end do

  end subroutine sortcps

  !> Generalized gradient tracing routine. The gp integration starts
  !> at xpoint (Cartesian). iup = 1 if the gp is traced up the
  !> density, -1 if down. nstep = number of steps (output) in the
  !> gradient path. ier = 0 (correct), 1 (short step), 2 (too many
  !> iterations), 3 (outside molcell in molecules). up2beta = .true.,
  !> stop when reaching a beta-sphere. plen, length of the gradient
  !> path. If path is present, return the gradient path. If prune is
  !> present, prune the gradient path to have point-to-point distances
  !> equal to prune. prune also sets the maximum step size. If pathini is
  !> present, initialize the output path with this point (Cartesian).
  module subroutine gradient(fid,xpoint,iup,nstep,ier,up2beta,plen,path,prune,pathini)
    use global, only: nav_step, nav_gradeps, rbetadef
    use types, only: scalar_value, gpathp, realloc
    use param, only: icrd_crys
    class(field), intent(inout) :: fid
    real*8, dimension(3), intent(inout) :: xpoint
    integer, intent(in) :: iup
    integer, intent(out) :: nstep
    integer, intent(out) :: ier
    logical, intent(in) :: up2beta
    real*8, intent(out) :: plen
    type(gpathp), intent(inout), allocatable, optional :: path(:)
    real*8, intent(in), optional :: prune
    real*8, intent(in), optional :: pathini(3)

    real*8, parameter :: minstep = 1d-12
    integer, parameter :: mhist = 5
    integer, parameter :: mstep0 = 10000
    real*8, parameter :: hfirst = 0.01d0

    integer :: mstep
    integer :: i, npath
    real*8 :: h0, hini
    real*8 :: dx(3), scalhist(mhist)
    real*8 :: xlast(3), xlast2(3), dd, prune0, xcart(3), xcaux(3), xxaux(3)
    real*8 :: sphrad, xnuc(3), xnucr(3), cprad, xcp(3), xcpr(3), dist2
    integer :: idnuc, idcp, nhist
    integer :: lvec(3)
    logical :: ok, incstep
    type(scalar_value) :: res, resaux

    ! initialization
    ier = 0
    hini = abs(NAV_step) * iup
    h0 = hfirst * iup
    scalhist = 1d0
    nhist = 0
    xlast2 = xpoint
    xlast = xpoint
    npath = 0
    plen = 0d0
    prune0 = -1d0
    if (present(prune)) then
       prune0 = prune
       hini = min(abs(hini),1.1d0*prune0) * iup
    end if
    incstep = .false.
    mstep = nint(mstep0 * 0.1d0 / abs(hini))

    if (present(path)) then
       if (allocated(path)) deallocate(path)
       allocate(path(100))
    end if
    if (present(pathini)) then
       xcaux = pathini
       xxaux = fid%c%c2x(xcaux)
       if (present(path)) then
          call fid%grd(xcaux,2,resaux)
          call addtopath(0,xcaux,xxaux,resaux)
       end if
       plen = norm2(xpoint - pathini)
    end if

    ! properties at point
    call fid%grd(xpoint,2,res)

    do nstep = 1, mstep
       ! point in cartesian and crystallographic. res contains evaluation at
       ! this point
       xcart = xpoint
       xpoint = fid%c%c2x(xpoint)
       plen = plen + norm2(xcart-xlast)

       ! save to the gradient path
       if (present(path)) then
          if (nstep == 1) then
             call addtopath(nstep,xcart,xpoint,res)
          else
             dd = norm2(xcart - path(npath)%r)
             if (dd > prune0) then
                if (prune0 < 0d0) then
                   call addtopath(nstep,xcart,xpoint,res)
                else
                   do while (dd > prune0)
                      xcaux = path(npath)%r + (xcart - path(npath)%r) * prune0 / dd
                      xxaux = fid%c%c2x(xcaux)
                      call fid%grd(xcaux,2,resaux)
                      call addtopath(nstep,xcaux,xxaux,resaux)
                      dd = norm2(xcart - path(npath)%r)
                   end do
                end if
             end if
          end if
       end if

       ! nearest nucleus
       call fid%c%nearest_atom(xpoint,icrd_crys,idnuc,sphrad,lvec=lvec)
       xnuc = fid%c%x2c(fid%c%atcel(idnuc)%x + lvec)
       xnucr = fid%c%atcel(idnuc)%x + lvec
       idnuc = fid%c%atcel(idnuc)%idx

       ! get nearest -3 CP (idncp) and +3 CP (idccp), skip hydrogens
       if ((fid%typnuc==-3 .and. iup==1 .or. fid%typnuc==3 .and. iup==-1) .and.&
          fid%c%spc(fid%c%at(idnuc)%is)%z /= 1) then
          idcp = idnuc
          cprad = sphrad
          xcp = xnuc
          xcpr = xnucr
       else
          idcp = 0
          cprad = 1d15
          xcp = 0d0
          xcpr = 0d0
       end if
       if (fid%ncpcel > 0) then
          do i = 1, fid%ncpcel
             if (.not.(fid%cpcel(i)%typ==-3 .and. iup==1 .or. fid%cpcel(i)%typ==3 .and. iup==-1)) cycle  ! only cages if down, only ncps if up
             dx = xpoint - fid%cpcel(i)%x
             call fid%c%shortest(dx,dist2)
             if (dist2 < cprad) then
                dx = fid%c%c2x(dx)
                cprad = dist2
                idcp = fid%cpcel(i)%idx
                xcp = fid%c%x2c(fid%cpcel(i)%x + nint(xpoint-fid%cpcel(i)%x-dx))
                xcpr = fid%cpcel(i)%x + nint(xpoint-fid%cpcel(i)%x-dx)
             end if
          end do
       end if

       ! is it a nuclear position?
       ok = .false.
       ! beta-sphere if up2beta is activated
       if (up2beta .and. idcp/=0) then
          ok = ok .or. (cprad <= fid%cp(idcp)%rbeta)
       else
          ok = ok .or. (cprad <= min(Rbetadef,0.5d0*abs(hini)))
       end if
       if (ok) then
          ! Found a CP as the last point in the path
          incstep = .true.
          plen = plen + norm2(xcart - xcp)
          xpoint = xcp
          ier = 0
          if (present(path)) then
             xcaux = xcp
             xxaux = xcpr
             call fid%grd(xcaux,2,resaux)
             call addtopath(nstep+1,xcaux,xxaux,resaux)
          end if
          goto 999
       end if

       ! is it outside the molcell?
       if (fid%c%ismolecule .and. iup==-1) then
          if (xpoint(1) < fid%c%molborder(1) .or. xpoint(1) > (1d0-fid%c%molborder(1)) .or.&
             xpoint(2) < fid%c%molborder(2) .or. xpoint(2) > (1d0-fid%c%molborder(2)) .or.&
             xpoint(3) < fid%c%molborder(3) .or. xpoint(3) > (1d0-fid%c%molborder(3))) then
             ! The gradient path exited the molecular cell
             ier = 3
             goto 999
          end if
       endif

       ! tasks in cartesian
       xpoint = fid%c%x2c(xpoint)

       ! Is it a new CP?
       ok = (res%gfmod < NAV_gradeps)
       if (ok) then
          ! Found a CP that was not on the list
          ier = 0
          if (present(path)) then
             if (path(npath)%i /= nstep) then
                xcaux = xpoint
                xxaux = fid%c%c2x(xcaux)
                call fid%grd(xcaux,2,resaux)
                call addtopath(nstep,xcaux,xxaux,resaux)
             end if
          end if
          goto 999
       end if

       ! take step
       xlast2 = xlast
       xlast = xpoint
       ok = adaptive_stepper(fid,xpoint,h0,hini,NAV_gradeps,res)

       ! add to the trajectory angle history, terminate the gradient if the
       ! trajectory bounces around mhist times
       nhist = mod(nhist,mhist) + 1
       scalhist(nhist) = dot_product(xlast-xlast2,xpoint-xlast)
       ok = ok .and. .not.(all(scalhist < 0d0))

       if (.not.ok .or. abs(h0) < minstep) then
          ! The step is too small and the gradient path is bouncing around
          xpoint = xlast
          ier = 1
          goto 999
       end if
    end do

    nstep = mstep
    ier = 2

999 continue
    if (incstep) nstep = nstep + 1
    if (present(path)) call realloc(path,npath)

  contains
    subroutine addtopath(istep,cart,crys,res)
      integer, intent(in) :: istep
      real*8, intent(in) :: cart(3),crys(3)
      type(scalar_value), intent(in) :: res

      npath = npath + 1
      if (npath > size(path,1)) call realloc(path,2*npath)
      path(npath)%i = istep
      path(npath)%x = crys
      path(npath)%r = cart
      path(npath)%f = res%f
      path(npath)%gf = res%gf
      path(npath)%hf = res%hf
    end subroutine addtopath
  end subroutine gradient

  !xx! private procedures

  !> Integration using adaptive_stepper step.
  function adaptive_stepper(fid,xpoint,h0,maxstep,eps,res)
    use global, only: nav_stepper, nav_stepper_heun, nav_stepper_rkck, nav_stepper_dp,&
       nav_stepper_bs, nav_stepper_euler, nav_maxerr
    use types, only: scalar_value
    use param, only: vsmall
    logical :: adaptive_stepper
    type(field), intent(inout) :: fid
    real*8, intent(inout) :: xpoint(3)
    real*8, intent(inout) :: h0
    real*8, intent(in) :: maxstep, eps
    type(scalar_value), intent(inout) :: res

    integer :: ier
    real*8 :: grdt(3), ogrdt(3)
    real*8 :: xtemp(3), escalar, xerrv(3)
    real*8 :: nerr
    logical :: ok, first

    real*8, parameter :: h0break = 1.d-10

    adaptive_stepper = .true.
    ier = 1
    grdt = res%gf / (res%gfmod + VSMALL)
    first = .true.
    do while (ier /= 0)
       ! calculate the new grdt
       ogrdt = grdt

       ! new point
       if (NAV_stepper == NAV_stepper_euler) then
          call stepper_euler1(xpoint,grdt,h0,xtemp)
       else if (NAV_stepper == NAV_stepper_heun) then
          call stepper_heun(fid,xpoint,grdt,h0,xtemp,res)
       else if (NAV_stepper == NAV_stepper_rkck) then
          call stepper_rkck(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       else if (NAV_stepper == NAV_stepper_dp) then
          call stepper_dp(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       else if (NAV_stepper == NAV_stepper_bs) then
          call stepper_bs(fid,xpoint,grdt,h0,xtemp,xerrv,res)
       end if

       ! FSAL for BS stepper
       if (NAV_stepper /= NAV_stepper_bs) then
          call fid%grd(xtemp,2,res)
       end if
       grdt = res%gf / (res%gfmod + VSMALL)

       ! poor man's adaptive step size in Euler
       if (NAV_stepper == NAV_stepper_euler .or. NAV_stepper == NAV_stepper_heun) then
          ! angle with next step
          escalar = dot_product(ogrdt,grdt)

          ! gradient eps in cartesian
          ok = (res%gfmod < 0.99d0*eps)

          ! Check if they differ in > 90 deg.
          if (escalar < 0.d0.and..not.ok) then
             if (abs(h0) >= h0break) then
                h0 = 0.5d0 * h0
                ier = 1
             else
                adaptive_stepper = .false.
                return
             end if
          else
             ! Accept point. If angle is favorable, take longer steps
             if (escalar > 0.9 .and. first) &
                h0 = dsign(min(abs(maxstep), abs(1.6d0*h0)),maxstep)
             ier = 0
             xpoint = xtemp
          end if
       else
          ! use the error estimate
          nerr = norm2(xerrv)
          if (nerr < NAV_maxerr) then
             ! accept point
             ier = 0
             xpoint = xtemp
             ! if this is the first time through, and the norm is very small, propose a longer step
             if (first .and. nerr < NAV_maxerr/10d0) &
                h0 = dsign(min(abs(maxstep), abs(1.6d0*h0)),maxstep)
          else
             ! propose a new shorter step using the error estimate
             h0 = 0.9d0 * h0 * NAV_maxerr / nerr
             if (abs(h0) < VSMALL) then
                adaptive_stepper = .false.
                return
             end if
          endif
       end if
       first = .false.
    enddo

  end function adaptive_stepper

  !> Euler stepper.
  subroutine stepper_euler1(xpoint,grdt,h0,xout)

    real*8, intent(in) :: xpoint(3), h0, grdt(3)
    real*8, intent(out) :: xout(3)

    xout = xpoint + h0 * grdt

  end subroutine stepper_euler1

  !> Heun stepper.
  subroutine stepper_heun(fid,xpoint,grdt,h0,xout,res)
    use types, only: scalar_value
    use param, only: vsmall

    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), h0, grdt(3)
    real*8, intent(out) :: xout(3)
    type(scalar_value), intent(inout) :: res

    real*8 :: ak2(3)

    xout = xpoint + h0 * grdt

    call fid%grd(xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + 0.5d0 * h0 * (ak2 + grdt)

  end subroutine stepper_heun

  !> Bogacki-Shampine embedded 2(3) method, fsal
  subroutine stepper_bs(fid,xpoint,grdt,h0,xout,xerr,res)
    use types, only: scalar_value
    use param, only: vsmall

    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), h0, grdt(3)
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res

    real*8, dimension(3) :: ak1, ak2, ak3, ak4

    ak1 = grdt

    xout = xpoint + h0 * (0.5d0*ak1)
    call fid%grd(xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)

    xout = xpoint + h0 * (0.75d0*ak2)
    call fid%grd(xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)

    xout = xpoint + h0 * (2d0/9d0*ak1 + 1d0/3d0*ak2 + 4d0/9d0*ak3)
    call fid%grd(xout,2,res)
    ak4 = res%gf / (res%gfmod+VSMALL)

    xerr = xpoint + h0 * (7d0/24d0*ak1 + 1d0/4d0*ak2 + 1d0/3d0*ak3 + 1d0/8d0*ak4) - xout

  end subroutine stepper_bs

  !> Runge-Kutta-Cash-Karp embedded 4(5)-order, local extrapolation.
  subroutine stepper_rkck(fid,xpoint,grdt,h0,xout,xerr,res)
    use types, only: scalar_value
    use param, only: vsmall

    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), grdt(3), h0
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res

    real*8, parameter :: B21=.2d0, &
         B31=3.d0/40.d0, B32=9.d0/40.d0,&
         B41=.3d0, B42=-.9d0, B43=1.2d0,&
         B51=-11.d0/54.d0, B52=2.5d0, B53=-70.d0/27.d0, B54=35.d0/27.d0,&
         B61=1631.d0/55296.d0,B62=175.d0/512.d0, B63=575.d0/13824.d0, B64=44275.d0/110592.d0, B65=253.d0/4096.d0,&
         C1=37.d0/378.d0, C3=250.d0/621.d0, C4=125.d0/594.d0, C6=512.d0/1771.d0,&
         DC1=C1-2825.d0/27648.d0, DC3=C3-18575.d0/48384.d0, DC4=C4-13525.d0/55296.d0, DC5=-277.d0/14336.d0, DC6=C6-.25d0
    real*8, dimension(3) :: ak2, ak3, ak4, ak5, ak6

    xout = xpoint + h0*B21*grdt

    call fid%grd(xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B31*grdt+B32*ak2)

    call fid%grd(xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B41*grdt+B42*ak2+B43*ak3)

    call fid%grd(xout,2,res)
    ak4 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B51*grdt+B52*ak2+B53*ak3+B54*ak4)

    call fid%grd(xout,2,res)
    ak5 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(B61*grdt+B62*ak2+B63*ak3+B64*ak4+B65*ak5)

    call fid%grd(xout,2,res)
    ak6 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(C1*grdt+C3*ak3+C4*ak4+C6*ak6)
    xerr = h0*(DC1*grdt+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

  end subroutine stepper_rkck

  !> Doermand-Prince embedded 4(5)-order, local extrapolation.
  subroutine stepper_dp(fid,xpoint,grdt,h0,xout,xerr,res)
    use types, only: scalar_value
    use param, only: vsmall

    type(field), intent(inout) :: fid
    real*8, intent(in) :: xpoint(3), grdt(3), h0
    real*8, intent(out) :: xout(3), xerr(3)
    type(scalar_value), intent(inout) :: res

    real*8, parameter :: dp_a(7,7) = reshape((/&
       0.0d0,  0d0,0d0,0d0,0d0,0d0,0d0,&
       1d0/5d0,         0.0d0,0d0,0d0,0d0,0d0,0d0,&
       3d0/40d0,        9d0/40d0,       0.0d0,0d0,0d0,0d0,0d0,&
       44d0/45d0,      -56d0/15d0,      32d0/9d0,        0d0,0d0,0d0,0d0,&
       19372d0/6561d0, -25360d0/2187d0, 64448d0/6561d0, -212d0/729d0,  0d0,0d0,0d0,&
       9017d0/3168d0,  -355d0/33d0,     46732d0/5247d0,  49d0/176d0,  -5103d0/18656d0, 0d0,0d0,&
       35d0/384d0,      0d0,            500d0/1113d0,    125d0/192d0, -2187d0/6784d0,  11d0/84d0,      0d0&
       /),shape(dp_a))
    real*8, parameter :: dp_b2(7) = (/5179d0/57600d0, 0d0, 7571d0/16695d0, 393d0/640d0,&
       -92097d0/339200d0, 187d0/2100d0, 1d0/40d0/)
    real*8, parameter :: dp_b(7) = (/ 35d0/384d0, 0d0, 500d0/1113d0, 125d0/192d0, &
       -2187d0/6784d0, 11d0/84d0, 0d0 /)
    real*8, parameter :: dp_c(7) = dp_b2 - dp_b
    real*8, dimension(3) :: ak2, ak3, ak4, ak5, ak6, ak7

    xout = xpoint + h0*dp_a(2,1)*grdt

    call fid%grd(xout,2,res)
    ak2 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(3,1)*grdt+dp_a(3,2)*ak2)

    call fid%grd(xout,2,res)
    ak3 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(4,1)*grdt+dp_a(4,2)*ak2+dp_a(4,3)*ak3)

    call fid%grd(xout,2,res)
    ak4 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(5,1)*grdt+dp_a(5,2)*ak2+dp_a(5,3)*ak3+dp_a(5,4)*ak4)

    call fid%grd(xout,2,res)
    ak5 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_a(6,1)*grdt+dp_a(6,2)*ak2+dp_a(6,3)*ak3+dp_a(6,4)*ak4+dp_a(6,5)*ak5)

    call fid%grd(xout,2,res)
    ak6 = res%gf / (res%gfmod+VSMALL)
    xout = xpoint + h0*(dp_b(1)*grdt+dp_b(2)*ak2+dp_b(3)*ak3+dp_b(4)*ak4+dp_b(5)*ak5+dp_b(6)*ak6)

    call fid%grd(xout,2,res)
    ak7 = res%gf / (res%gfmod+VSMALL)
    xerr = h0*(dp_c(1)*grdt+dp_c(2)*ak2+dp_c(3)*ak3+dp_c(4)*ak4+dp_c(5)*ak5+dp_c(6)*ak6+dp_c(7)*ak7)
    xout = xout + xerr

  end subroutine stepper_dp

end submodule proc
