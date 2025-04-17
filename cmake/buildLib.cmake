message(STATUS "build library ${lib}, shared=${BUILD_SHARED_LIBS}")

set(src_files ../lib/erfcx.c ../lib/err_fcts.c ../lib/im_w_of_x.c ../lib/w_of_z.c ../lib/width.c)
set(inc_files ../lib/cerf.h)

if (CERF_TARGET_CPP)
    set_property(SOURCE ${src_files} PROPERTY LANGUAGE CXX)
else()
    set_property(SOURCE ${src_files} PROPERTY LANGUAGE C)
endif()

add_library(${lib} ${src_files})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${CERF_VERSION}
    SOVERSION ${CERF_SOVERSION})
if (CERF_TARGET_CPP)
    set_target_properties(${lib} PROPERTIES LINKER_LANGUAGE CXX)
endif()

if(NOT MSVC)
    target_link_libraries(${lib} m)
endif()

target_include_directories(${lib} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lib>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(
    TARGETS ${lib}
    EXPORT ${intf}
    RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT Libraries)
install(
    FILES ${inc_files}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Headers)

include(CMakePackageConfigHelpers)
set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/cerf")
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/cerfConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cerfConfig.cmake"
  INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/cerfConfigVersion.cmake
  VERSION ${CERF_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cerfConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cerfConfigVersion.cmake
  DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(EXPORT ${intf}
  NAMESPACE "cerf::"
  FILE "cerfTargets.cmake"
  DESTINATION ${CMAKECONFIG_INSTALL_DIR})

set(CPACK_PACKAGE_VENDOR "Forschungszentrum Jülich GmbH")
set(CPACK_PACKAGE_VERSION ${CERF_VERSION})
set(CPACK_GENERATOR "ZIP")
include(CPack)
