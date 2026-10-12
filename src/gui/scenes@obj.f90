! Copyright (c) 2019-2022 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
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


! Export of a scene to a 3D model file (Wavefront obj + mtl). The
! spheres and cylinders are drawn as impostors, so their geometry is
! rebuilt here from the draw lists, in the world frame of the scene
! (critic2 Cartesian coordinates), converted to angstrom.
submodule (scenes) obj
  implicit none

  !xx! private procedures
  ! subroutine export_cone(g,x1,x2,r,rgb)
  ! function rgbint(rgb)

  ! number of segments of the exported cones
  integer, parameter :: ncone_seg = 24

contains

  !> Write the scene s, as currently drawn, to file (obj format, plus
  !> a mtl file with the same root). Text, the window-anchored objects,
  !> the 2D drawings and the legends/bars are not written. Returns a
  !> non-empty errmsg on error.
  module subroutine scene_export_obj(s,file,errmsg)
    use graphics, only: grhandle
    use param, only: bohrtoa, dirsep, mlen
    use tools_io, only: uout, string
    use tools_math, only: cross
    class(scene), intent(inout), target :: s
    character*(*), intent(in) :: file
    character(len=:), allocatable, intent(out) :: errmsg

    type(grhandle) :: g
    character(len=:), allocatable :: fobj, fmtl
    integer :: i, j, k, nseg, ndash, lu, ios, irep, nrep, nv, niso
    real(c_float) :: dl(3), r
    logical :: dsh(3)
    real*8 :: x1(3), x2(3), dir(3), vview(3), vperp(3), vout(3), blen, p, tc
    real*8 :: alpha, xq(3,4), xoff(3)
    real*8, allocatable :: xm(:,:), nm(:,:), cm(:,:)
    integer, allocatable :: im(:,:)
    logical :: exist

    errmsg = ""
    if (s%isinit == 0) then
       errmsg = "the scene is not initialized"
       return
    end if

    ! file names: add the extension if there is none (the dialog did
    ! not ask about overwriting that name)
    fobj = trim(file)
    if (index(fobj(index(fobj,dirsep,back=.true.)+1:),'.',back=.true.) == 0) then
       fobj = fobj // ".obj"
       inquire(file=fobj,exist=exist)
       if (exist) then
          errmsg = "file exists (give the name with its extension to replace it): " // fobj
          return
       end if
    end if
    fmtl = fobj(:index(fobj,'.',back=.true.)-1) // ".mtl"
    if (len(fobj) > mlen .or. len(fmtl) > mlen) then
       errmsg = "file name too long: " // fobj
       return
    end if

    ! check that both files can be written, without truncating them yet
    ! (the graphics writer stops the program if it cannot open them)
    open(newunit=lu,file=fobj,status='unknown',action='write',position='append',iostat=ios)
    if (ios == 0) then
       close(lu)
       open(newunit=lu,file=fmtl,status='unknown',action='write',position='append',iostat=ios)
       if (ios == 0) then
          close(lu)
       else
          errmsg = "could not open file for writing: " // fmtl
          return
       end if
    else
       errmsg = "could not open file for writing: " // fobj
       return
    end if

    ! build the draw lists if they are not current
    if (s%isinit == 1 .or. .not.allocated(s%obj%sph) .or. s%forcebuildlists) &
       call s%build_lists()

    ! the view direction in the world frame, for the lateral offset of the
    ! double and triple bonds (as on the screen; the world matrix is a
    ! rotation and a translation)
    vview = matmul(transpose(real(s%world(1:3,1:3),8)),real(s%camfront,8))
    if (norm2(vview) > 1d-10) then
       vview = vview / norm2(vview)
    else
       vview = (/0d0,0d0,-1d0/)
    end if

    call g%open("obj",fobj)

    ! spheres (atoms, critical points...)
    do i = 1, s%obj%nsph
       if (s%obj%sph(i)%ghost) cycle
       call g%ball(real(s%obj%sph(i)%x,8)*bohrtoa,rgbint(s%obj%sph(i)%rgb),&
          real(s%obj%sph(i)%r,8)*bohrtoa,real(s%obj%sph(i)%alpha,8))
    end do

    ! bonds: the bond order expanded into parallel (possibly dashed) runs,
    ! as in the renderer
    do i = 1, s%obj%ncyl
       x1 = s%obj%cyl(i)%x1
       x2 = s%obj%cyl(i)%x2
       blen = norm2(x2-x1)
       if (blen < 1d-7) cycle
       dir = (x2 - x1) / blen
       r = s%obj%cyl(i)%r
       alpha = s%obj%cyl(i)%alpha
       call cyl_bond_runs(s%obj%cyl(i)%order,r,nseg,dl,dsh)

       ! lateral direction: the ring exterior for aromatic bonds,
       ! perpendicular to the bond and the view direction otherwise
       vperp = 0d0
       if (nseg > 1) then
          vout = s%obj%cyl(i)%arvec
          if (s%obj%cyl(i)%order == -1 .and. norm2(vout) > 1d-10) then
             vperp = vout - dot_product(vout,dir) * dir
          else
             vperp = cross(vview,dir)
          end if
          if (norm2(vperp) < 1d-10) then
             ! the bond is along the view: any perpendicular direction
             if (abs(dir(1)) < 0.9d0) then
                vperp = cross((/1d0,0d0,0d0/),dir)
             else
                vperp = cross((/0d0,1d0,0d0/),dir)
             end if
          end if
          vperp = vperp / norm2(vperp)
       end if

       do k = 1, nseg
          xoff = dl(k) * vperp
          if (dsh(k)) then
             ndash = cyl_dash_count(s%obj%cyl(i))
             p = blen / real(ndash,8)
             do j = 0, ndash-1
                tc = (real(j,8) + 0.5d0) * p
                call g%stick((x1+xoff+dir*(tc-0.25d0*p))*bohrtoa,(x1+xoff+dir*(tc+0.25d0*p))*bohrtoa,&
                   rgbint(s%obj%cyl(i)%rgb),0.5d0*r*bohrtoa,alpha)
             end do
          else
             call g%stick((x1+xoff)*bohrtoa,(x2+xoff)*bohrtoa,rgbint(s%obj%cyl(i)%rgb),&
                0.5d0*r*bohrtoa,alpha)
          end if
       end do
    end do

    ! flat cylinders (unit cell and polyhedra edges)
    do i = 1, s%obj%ncylflat
       call g%stick(real(s%obj%cylflat(i)%x1,8)*bohrtoa,real(s%obj%cylflat(i)%x2,8)*bohrtoa,&
          rgbint(s%obj%cylflat(i)%rgb),0.5d0*real(s%obj%cylflat(i)%r,8)*bohrtoa)
    end do

    ! cones (arrowheads)
    do i = 1, s%obj%ncone
       call export_cone(g,real(s%obj%cone(i)%x1,8)*bohrtoa,real(s%obj%cone(i)%x2,8)*bohrtoa,&
          0.5d0*real(s%obj%cone(i)%r,8)*bohrtoa,rgbint(s%obj%cone(i)%rgb))
    end do

    ! flat rectangles
    do i = 1, s%obj%nplane
       xq(:,1) = s%obj%plane(i)%x - s%obj%plane(i)%e1 - s%obj%plane(i)%e2
       xq(:,2) = s%obj%plane(i)%x + s%obj%plane(i)%e1 - s%obj%plane(i)%e2
       xq(:,3) = s%obj%plane(i)%x + s%obj%plane(i)%e1 + s%obj%plane(i)%e2
       xq(:,4) = s%obj%plane(i)%x - s%obj%plane(i)%e1 + s%obj%plane(i)%e2
       call g%polygon(xq*bohrtoa,rgbint(s%obj%plane(i)%rgb),real(s%obj%plane(i)%alpha,8))
    end do

    ! flat triangles (polyhedra faces)
    do i = 1, s%obj%ntriangle
       xq(:,1) = s%obj%triangle(i)%x1
       xq(:,2) = s%obj%triangle(i)%x2
       xq(:,3) = s%obj%triangle(i)%x3
       call g%polygon(xq(:,1:3)*bohrtoa,rgbint(s%obj%triangle(i)%rgb),&
          real(s%obj%triangle(i)%alpha,8))
    end do

    ! triangle meshes (isosurfaces), one per drawn copy
    niso = 0
    do i = 1, s%obj%nmsh
       nv = s%obj%msh(i)%nv
       if (nv <= 0 .or. s%obj%msh(i)%nf <= 0) cycle
       nm = real(s%obj%msh(i)%nrm(:,1:nv),8)
       im = s%obj%msh(i)%idx(:,1:s%obj%msh(i)%nf) + 1
       if (allocated(cm)) deallocate(cm)
       if (allocated(s%obj%msh(i)%rgbv)) cm = real(s%obj%msh(i)%rgbv(:,1:nv),8)
       nrep = 1
       if (allocated(s%obj%msh(i)%xrep)) nrep = size(s%obj%msh(i)%xrep,2)
       do irep = 1, nrep
          xoff = 0d0
          if (allocated(s%obj%msh(i)%xrep)) xoff = s%obj%msh(i)%xrep(:,irep)
          xm = (real(s%obj%msh(i)%x(:,1:nv),8) + spread(xoff,2,nv)) * bohrtoa
          ! an unallocated cm is an absent optional argument
          call g%mesh(nv,xm,nm,s%obj%msh(i)%nf,im,rgbint(s%obj%msh(i)%rgb),&
             real(s%obj%msh(i)%alpha,8),cm)
          niso = niso + 1
       end do
    end do

    write (uout,'("* Scene written to: ",A)') fobj
    write (uout,'("  Spheres: ",A,", sticks: ",A,", cones: ",A,", surfaces: ",A,", materials: ",A/)') &
       string(g%nball), string(g%nstick), string(s%obj%ncone), string(niso), string(g%nmtl)
    call g%close()

  end subroutine scene_export_obj

  !xx! private procedures

  !> Write a cone with base center x1, apex x2, and base radius r, as a
  !> mesh with normals (side and base).
  subroutine export_cone(g,x1,x2,r,rgb)
    use graphics, only: grhandle
    use param, only: pi
    use tools_math, only: cross
    type(grhandle), intent(inout) :: g
    real*8, intent(in) :: x1(3), x2(3), r
    integer, intent(in) :: rgb(3)

    real*8 :: ax(3), u(3), v(3), h, ang, cs(3), slope
    real*8 :: x(3,3*ncone_seg+1), nrm(3,3*ncone_seg+1)
    integer :: idx(3,2*ncone_seg)
    integer :: i, i1

    h = norm2(x2 - x1)
    if (h < 1d-7 .or. r < 1d-7) return
    ax = (x2 - x1) / h

    ! two unit vectors perpendicular to the axis
    if (abs(ax(1)) < 0.9d0) then
       u = cross((/1d0,0d0,0d0/),ax)
    else
       u = cross((/0d0,1d0,0d0/),ax)
    end if
    u = u / norm2(u)
    v = cross(ax,u)
    slope = r / h

    ! vertices: the side ring (base points with side normals), the apex
    ! copies (one per segment), the base ring (with the base normal), and
    ! the base center
    do i = 1, ncone_seg
       ang = 2d0 * pi * real(i-1,8) / real(ncone_seg,8)
       cs = cos(ang) * u + sin(ang) * v
       x(:,i) = x1 + r * cs
       nrm(:,i) = (cs + slope * ax) / norm2(cs + slope * ax)
       ang = 2d0 * pi * (real(i,8) - 0.5d0) / real(ncone_seg,8)
       cs = cos(ang) * u + sin(ang) * v
       x(:,ncone_seg+i) = x2
       nrm(:,ncone_seg+i) = (cs + slope * ax) / norm2(cs + slope * ax)
       x(:,2*ncone_seg+i) = x(:,i)
       nrm(:,2*ncone_seg+i) = -ax
    end do
    x(:,3*ncone_seg+1) = x1
    nrm(:,3*ncone_seg+1) = -ax

    ! faces: side, then base
    do i = 1, ncone_seg
       i1 = mod(i,ncone_seg) + 1
       idx(:,i) = (/i, i1, ncone_seg+i/)
       idx(:,ncone_seg+i) = (/3*ncone_seg+1, 2*ncone_seg+i1, 2*ncone_seg+i/)
    end do
    call g%mesh(3*ncone_seg+1,x,nrm,2*ncone_seg,idx,rgb)

  end subroutine export_cone

  !> Color 0-1 to 0-255.
  function rgbint(rgb)
    real(c_float), intent(in) :: rgb(3)
    integer :: rgbint(3)

    rgbint = max(min(nint(rgb * 255._c_float),255),0)

  end function rgbint

end submodule obj
