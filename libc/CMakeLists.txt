set(lib cerf)
set(intf "c_interface")
set(CERF_TARGET_CPP OFF)
set(cerf_LIBRARY ${lib} PARENT_SCOPE)

include("../cmake/buildLib.cmake")
