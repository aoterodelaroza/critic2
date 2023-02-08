! -*-F90-*-

!!! stb procedures !!!

!x! int stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes);
function stbi_write_png(filename,w,h,comp,data,stride_in_bytes) bind(c,name="stbi_write_png")
  import c_int, c_ptr
  integer(c_int) :: stbi_write_png
  type(c_ptr), value :: filename
  integer(c_int), value :: w
  integer(c_int), value :: h
  integer(c_int), value :: comp
  type(c_ptr), value :: data
  integer(c_int), value :: stride_in_bytes
end function stbi_write_png
!x! int stbi_write_bmp(char const *filename, int w, int h, int comp, const void *data);
function stbi_write_bmp(filename,w,h,comp,data) bind(c,name="stbi_write_bmp")
  import c_int, c_ptr
  integer(c_int) :: stbi_write_bmp
  type(c_ptr), value :: filename
  integer(c_int), value :: w
  integer(c_int), value :: h
  integer(c_int), value :: comp
  type(c_ptr), value :: data
end function stbi_write_bmp
!x! int stbi_write_tga(char const *filename, int w, int h, int comp, const void *data);
function stbi_write_tga(filename,w,h,comp,data) bind(c,name="stbi_write_tga")
  import c_int, c_ptr
  integer(c_int) :: stbi_write_tga
  type(c_ptr), value :: filename
  integer(c_int), value :: w
  integer(c_int), value :: h
  integer(c_int), value :: comp
  type(c_ptr), value :: data
end function stbi_write_tga
!x! int stbi_write_jpg(char const *filename, int w, int h, int comp, const void *data, int quality);
function stbi_write_jpg(filename,w,h,comp,data,quality) bind(c,name="stbi_write_jpg")
  import c_int, c_ptr
  integer(c_int) :: stbi_write_jpg
  type(c_ptr), value :: filename
  integer(c_int), value :: w
  integer(c_int), value :: h
  integer(c_int), value :: comp
  type(c_ptr), value :: data
  integer(c_int), value :: quality
end function stbi_write_jpg
!x! int stbi_write_hdr(char const *filename, int w, int h, int comp, const float *data);
function stbi_write_hdr(filename,w,h,comp,data) bind(c,name="stbi_write_hdr")
  import c_int, c_ptr
  integer(c_int) :: stbi_write_hdr
  type(c_ptr), value :: filename
  integer(c_int), value :: w
  integer(c_int), value :: h
  integer(c_int), value :: comp
  type(c_ptr), value :: data
end function stbi_write_hdr
!x! STBIDEF void stbi_image_free(void *retval_from_stbi_load)
subroutine stbi_image_free(retval_from_stbi_load) bind(c,name="stbi_image_free")
  import c_ptr
  type(c_ptr), value :: retval_from_stbi_load
end subroutine stbi_image_free
!x! STBIDEF stbi_uc *stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp)
function stbi_load(filename,x,y,comp,req_comp) bind(c,name="stbi_load")
  import c_ptr, c_int
  type(c_ptr), value :: filename
  integer(c_int) :: x, y, comp
  integer(c_int), value :: req_comp
  type(c_ptr) :: stbi_load
end function stbi_load
