message(STATUS "build library ${lib}, shared=${BUILD_SHARED_LIBS}")

set(src_files ../lib/erfcx.c ../lib/err_fcts.c ../lib/im_w_of_x.c ../lib/w_of_z.c ../lib/width.c)

if(CERF_TARGET_CPP)
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
if(CERF_TARGET_CPP)
    set_target_properties(${lib} PROPERTIES LINKER_LANGUAGE CXX)
endif()

if(NOT MSVC)
    target_link_libraries(${lib} m)
endif()

target_include_directories(${lib} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lib>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# will be used by 'install(EXPORT ...)'
install(
    TARGETS ${lib}
    EXPORT interface
    RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT Libraries)
