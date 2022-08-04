  ! -*-F90-*-

  !xx! cimgui constants

  ! enum ImGuiConfigFlags_
  integer(c_int), parameter :: ImGuiConfigFlags_None = 0
  integer(c_int), parameter :: ImGuiConfigFlags_NavEnableKeyboard = int(b'0000000000000000000001') ! 1 << 0
  integer(c_int), parameter :: ImGuiConfigFlags_NavEnableGamepad = int(b'0000000000000000000010') ! 1 << 1
  integer(c_int), parameter :: ImGuiConfigFlags_NavEnableSetMousePos = int(b'0000000000000000000100') ! 1 << 2
  integer(c_int), parameter :: ImGuiConfigFlags_NavNoCaptureKeyboard = int(b'0000000000000000001000') ! 1 << 3,
  integer(c_int), parameter :: ImGuiConfigFlags_NoMouse = int(b'0000000000000000010000') ! 1 << 4,
  integer(c_int), parameter :: ImGuiConfigFlags_NoMouseCursorChange = int(b'0000000000000000100000') ! 1 << 5,
  integer(c_int), parameter :: ImGuiConfigFlags_DockingEnable = int(b'0000000000000001000000') ! 1 << 6,
  integer(c_int), parameter :: ImGuiConfigFlags_ViewportsEnable = int(b'0000000000010000000000') ! 1 << 10,
  integer(c_int), parameter :: ImGuiConfigFlags_DpiEnableScaleViewports = int(b'0000000100000000000000') ! 1 << 14,
  integer(c_int), parameter :: ImGuiConfigFlags_DpiEnableScaleFonts = int(b'0000001000000000000000') ! 1 << 15,
  integer(c_int), parameter :: ImGuiConfigFlags_IsSRGB = int(b'0100000000000000000000') ! 1 << 20,
  integer(c_int), parameter :: ImGuiConfigFlags_IsTouchScreen = int(b'1000000000000000000000') ! 1 << 21

  ! enum ImGuiKey_
  integer(c_int), parameter :: ImGuiKey_COUNT = 645
  integer(c_int), parameter :: ImGuiKey_KeysData_SIZE = 645

  ! enum ImGuiNavInput_
  integer(c_int), parameter :: ImGuiNavInput_COUNT = 20
