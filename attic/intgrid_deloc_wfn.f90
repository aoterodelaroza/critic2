    ! ! localization and delocalization indices - molecular wfn case
    ! call intgrid_deloc_wfn(res,nattr,xgatt,idg,imtype,luw)

    ! ! localization and delocalization indices, molecule wfn
    ! call int_output_deloc_wfn(nattr,icp,sij)

  !> Calculate localization and delocalization indices using the
  !> basin assignment found by YT or BADER and a wfn scalar field.
  subroutine intgrid_deloc_wfn(res,natt,xgatt,idg,imtype,luw)
    use yt, only: yt_weights
    use systemmod, only: sy, itype_deloc
    use fieldmod, only: type_wfn
    use tools_io, only: fopen_scratch, fclose
    type(int_result), intent(inout) :: res(:)
    integer, intent(in) :: natt
    real*8, intent(in) :: xgatt(3,natt)
    integer, intent(in) :: idg(:,:,:)
    integer, intent(in) :: imtype
    integer, intent(in) :: luw

    integer :: i, j, k, l, m, i1, i2, ndeloc
    integer :: fid, n(3), ntot, lumo, ix
    real*8, allocatable :: w(:,:,:)
    real*8, allocatable :: xmo(:)
    real*8 :: x(3), rho, auxg(3), auxh(3,3), gkin, vir
    real*8 :: stress(3,3)

    ! only wfn fields
    ndeloc = 0
    do l = 1, sy%npropi
       if (.not.sy%propi(l)%used) cycle
       if (.not.res(l)%int == int_delocwfn) cycle
       fid = sy%propi(l)%fid
       if (sy%f(fid)%type /= type_wfn) cycle
       ndeloc = ndeloc + 1

       if (allocated(res(l)%sij)) deallocate(res(l)%sij)
       allocate(res(l)%sij(sy%f(fid)%wfn%nmoocc,sy%f(fid)%wfn%nmoocc,natt,1))
    end do
    if (ndeloc == 0) return

    ! size of the grid
    n = sy%f(sy%iref)%grid%n
    ntot = n(1)*n(2)*n(3)

    if (imtype == imtype_yt) then
       allocate(w(n(1),n(2),n(3)))
    endif

    ! run over all properties for which multipole calculation is active
    do l = 1, sy%npropi
       if (.not.sy%propi(l)%used) cycle
       if (.not.sy%propi(l)%itype == itype_deloc) cycle
       fid = sy%propi(l)%fid
       if (sy%f(fid)%type /= type_wfn) cycle

       ! calculate the overlap matrix
       allocate(xmo(sy%f(fid)%wfn%nmoocc))
       res(l)%sij = 0d0
       xmo = 0d0
       if (imtype == imtype_bader) then
          !$omp parallel do private (i,j,k,x,rho,auxg,auxh,ix,i1,i2,gkin,vir) firstprivate(xmo) schedule(dynamic)
          do i = 1, n(1)
             do j = 1, n(2)
                do k = 1, n(3)
                   x = sy%c%x2c(real((/i-1,j-1,k-1/),8) / n)
                   call sy%f(fid)%wfn%rho2(x,0,rho,auxg,auxh,gkin,vir,stress,xmo)
                   ix = idg(i,j,k)
                   do i1 = 1, sy%f(fid)%wfn%nmoocc
                      do i2 = i1, sy%f(fid)%wfn%nmoocc
                         !$omp critical (sijwrite)
                         res(l)%sij(i1,i2,ix,1) = res(l)%sij(i1,i2,ix,1) + xmo(i1) * xmo(i2)
                         !$omp end critical (sijwrite)
                      end do
                   end do
                end do
             end do
          end do
          !$omp end parallel do
       else
          ! write the MO vectors to a file
          lumo = fopen_scratch()
          do i = 1, n(1)
             do j = 1, n(2)
                do k = 1, n(3)
                   x = sy%c%x2c(real((/i-1,j-1,k-1/),8) / n)
                   call sy%f(fid)%wfn%rho2(x,0,rho,auxg,auxh,gkin,vir,stress,xmo)
                   ix = idg(i,j,k)
                   write (lumo) xmo
                end do
             end do
          end do

          ! calculate the atomic overlap matrix and kill the file
          do m = 1, natt
             call yt_weights(luw=luw,idb=m,w=w)
             rewind(lumo)
             do i = 1, n(1)
                do j = 1, n(2)
                   do k = 1, n(3)
                      read(lumo) xmo
                      do i1 = 1, sy%f(fid)%wfn%nmoocc
                         do i2 = i1, sy%f(fid)%wfn%nmoocc
                            res(l)%sij(i1,i2,m,1) = res(l)%sij(i1,i2,m,1) + xmo(i1) * xmo(i2) * w(i,j,k)
                         end do
                      end do
                   end do
                end do
             end do
          end do
          call fclose(lumo)
       end if
       deallocate(xmo)

       ! symmetrize and scale
       res(l)%sij(:,:,:,:) = res(l)%sij(:,:,:,:) * sy%c%omega / ntot
       do ix = 1, natt
          do i1 = 1, sy%f(fid)%wfn%nmoocc
             do i2 = i1, sy%f(fid)%wfn%nmoocc
                res(l)%sij(i2,i1,ix,1) = res(l)%sij(i1,i2,ix,1)
             end do
          end do
       end do
    end do
    if (allocated(w)) deallocate(w)

  end subroutine intgrid_deloc_wfn

  !> Output delocalization indices in a molecular/wfn calculation.
  !> Inputs: number of attractors (nattr), CP identifiers for the
  !> attractors (icp), and the atomic overlap matrix (sij).
  subroutine int_output_deloc_wfn(nattr,icp,sij)
    use systemmod, only: sy, itype_deloc
    use fieldmod, only: type_wfn
    use wfn_private, only: wfn_rhf
    use tools_io, only: uout, string, ioj_left
    integer, intent(in) :: nattr
    integer, intent(in) :: icp(nattr)
    real*8, intent(in), allocatable, optional :: sij(:,:,:,:,:)

    integer :: i, j, l, fid
    integer :: ndeloc
    character(len=:), allocatable :: sncp, scp, sname, sz, smult
    character(len=:), allocatable :: sout1, sout2
    real*8, allocatable :: di(:,:)
    real*8 :: fspin

    if (.not.present(sij)) return
    if (.not.allocated(sij)) return
    if (.not.sy%c%ismolecule) return

    ! space for the DIs
    allocate(di(nattr,nattr))

    ndeloc = 0
    do l = 1, sy%npropi
       if (.not.sy%c%ismolecule) cycle
       if (.not.sy%propi(l)%used) cycle
       if (.not.sy%propi(l)%itype == itype_deloc) cycle
       fid = sy%propi(l)%fid
       if (sy%f(fid)%type /= type_wfn) cycle
       ndeloc = ndeloc + 1

       ! calculate localization and delocalization indices using
       ! the atomic overlap matrix
       di = 0d0
       if (sy%f(fid)%wfn%wfntyp == wfn_rhf) then
          fspin = 4d0
       else
          fspin = 1d0
       end if
       do i = 1, nattr
          do j = i, nattr
             di(i,j) = fspin * sum(sij(:,:,i,1,ndeloc)*sij(:,:,j,1,ndeloc))
             di(j,i) = di(i,j)
          end do
       end do

       ! Header
       write (uout,'("* Localization and delocalization indices")')
       write (uout,'("+ Integrated property (number ",A,"): ",A)') string(l), string(sy%propi(l)%prop_name)

       ! output the lambdas
       write (uout,'("+ Localization indices (lambda(A))")')
       write (uout,'("# Id   cp   ncp   Name  Z      lambda(A)  ")')
       do j = 1, nattr
          call assign_strings(j,icp(j),.false.,scp,sncp,sname,smult,sz)
          write (uout,'(2X,99(A,X))') &
             string(j,4,ioj_left), scp, sncp, sname, sz, &
             string(0.5d0*di(j,j),'e',15,8,4)
       enddo
       write (uout,*)

       write (uout,'("+ Delocalization indices (delta(A,B))")')
       write (uout,'("#     ----- atom A -----          ----- atom B -----                       ")')
       write (uout,'("#   Id   cp   ncp   Name  Z     Id   cp   ncp   Name  Z        delta(A,B)  ")')
       do i = 1, nattr
          call assign_strings(i,icp(i),.false.,scp,sncp,sname,smult,sz)
          sout1 = "| " // string(i,4,ioj_left) // " " // scp // " " // sncp // " " // sname // " " // sz // " " // " |"
          do j = i+1, nattr
             call assign_strings(j,icp(j),.false.,scp,sncp,sname,smult,sz)
             sout2 = string(j,4,ioj_left) // " " // scp // " " // sncp // " " // sname // " " // sz // " " // " |"
             write (uout,'(2X,99(A,X))') string(sout1), string(sout2), &
                string(di(i,j),'e',15,8,4)
          end do
       end do
       write (uout,*)
    end do
          
    ! clean up
    deallocate(di)

  end subroutine int_output_deloc_wfn

