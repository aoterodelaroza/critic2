submodule (config) proc

  character(len=*), parameter :: package = "critic2"
#ifdef COMMIT
  character(len=*), parameter :: version = COMMIT
#else
  character(len=*), parameter :: version = "5aae766"
#endif
  character(len=*), parameter :: f77 = "gfortran"
  character(len=*), parameter :: fc = "gfortran"
  character(len=*), parameter :: cc = "gcc"
  character(len=*), parameter :: fflags = "-g -O2 "
  character(len=*), parameter :: fcflags = "-g -O2 -ffree-form -ffree-line-length-none -fopenmp"
  character(len=*), parameter :: cflags = "-g -O2 -Iqhull/"
  character(len=*), parameter :: ldflags = ""
  character(len=*), parameter :: atarget = "Linux oberon 4.13.0-1-amd64 #1 SMP Debian 4.13.4-2 (2017-10-15) x86_64 GNU/Linux"
#ifdef DATE
  character(len=*), parameter :: adate = DATE
#else
  character(len=*), parameter :: adate = "Sat Jan 13 12:14:23 PST 2018"
#endif
  character(len=*), parameter :: enable_debug = "no"
  character(len=*), parameter :: revision = "@AC_REVISION@"
#ifdef DATADIR
  character(len=*), parameter :: datadir = trim(adjustl(DATADIR)) // dirsep // "critic2"
#else
  character(len=*), parameter :: datadir = "."
#endif

contains

  !> Return the config string corresponding to the input id.
  function getstring(id) result(str)
    integer, intent(in) :: id
    character(len=:), allocatable :: str

    select case(id)
    case(istring_package)
       str = package
    case(istring_version)
       str = version
    case(istring_f77)
       str = f77
    case(istring_fc)
       str = fc
    case(istring_cc)
       str = cc
    case(istring_fflags)
       str = fflags
    case(istring_fcflags)
       str = fcflags
    case(istring_cflags)
       str = cflags
    case(istring_ldflags)
       str = ldflags
    case(istring_atarget)
       str = atarget
    case(istring_adate)
       str = adate
    case(istring_enabledebug)
       str = enable_debug
    case(istring_revision)
       str = revision
    case(istring_datadir)
       str = datadir
    end select

  end function getstring
  
end submodule proc
