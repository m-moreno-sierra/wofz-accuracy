## Install header

install(
    FILES ${CMAKE_SOURCE_DIR}/lib/cerf.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT Headers)

## Make and install pkg_config ressource file

configure_file(var/libcerf.pc.in libcerf.pc @ONLY)
install(FILES "${CMAKE_BINARY_DIR}/libcerf.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

## Make and install CMake *Config* files

include(CMakePackageConfigHelpers)
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/cerfConfig.cmake.in
    "${CMAKE_BINARY_DIR}/cerfConfig.cmake"
    INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_BINARY_DIR}/cerfConfigVersion.cmake
    VERSION ${CERF_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_BINARY_DIR}/cerfConfig.cmake
    ${CMAKE_BINARY_DIR}/cerfConfigVersion.cmake
    DESTINATION ${CMAKECONFIG_INSTALL_DIR})

install(EXPORT interface
    NAMESPACE "cerf::"
    FILE "cerfTargets.cmake"
    DESTINATION ${CMAKECONFIG_INSTALL_DIR})

## Make zip package

set(CPACK_PACKAGE_VENDOR "Forschungszentrum Jülich GmbH")
set(CPACK_PACKAGE_VERSION ${CERF_VERSION})
set(CPACK_GENERATOR "ZIP")
include(CPack)
