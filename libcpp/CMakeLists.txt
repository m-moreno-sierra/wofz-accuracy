set(lib cerfcpp)
set(intf "cpp_interface")
set(CERF_TARGET_CPP ON)
add_compile_definitions("CERF_AS_CPP=ON")
set(cerfcpp_LIBRARY ${lib} PARENT_SCOPE)

include("../cmake/buildLib.cmake")
