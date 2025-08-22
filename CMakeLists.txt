cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
if(NOT DEFINED LIB_MAN)
    option(LIB_MAN "Build libcerf man pages" ON)
endif()
if(NOT DEFINED LIB_RUN)
    option(LIB_RUN "Build executables for command-line computation" ON)
endif()
if(NOT DEFINED CERF_C)
    option(CERF_C "Build C library libcerf" ON)
endif()
if(NOT DEFINED CERF_CPP)
    option(CERF_CPP "Build C++ libcerfcpp" ON)
endif()
if(NOT DEFINED CERF_IEEE754)
    option(CERF_IEEE754 "Assume IEEE754 floating-point layout" ON)
    # Use our frexp2 instead of frexp from libc.
    # It is unclear though whether the remainder of the code will work on non-IEEE754 machines.
endif()

project(cerf LANGUAGES NONE)
if (CERF_C)
    if(WIN32)
        message(FATAL_ERROR "Compilation as C not supported by MSVC, use -DCERF_C=OFF")
    endif()
    enable_language(C)
    set(CMAKE_C_STANDARD 11)
endif()
if (CERF_CPP)
    enable_language(CXX)
    set(CMAKE_CXX_STANDARD 17)
endif()

set(CERF_SOVERSION                 3) # API version
set(CERF_VERSION ${CERF_SOVERSION}.1) # minor version

include(GNUInstallDirs)

SET(CERF_COMPILE_OPTIONS "" CACHE
    STRING "User-supplied compiler options, will be appended to hard-coded ones. Not for Windows.")

if (CERF_IEEE754)
    add_compile_definitions(CERF_IEEE754)
    if(DEFINED ENDIAN_IS_LITTLE)
        # Normalize value to 0 or 1
        if(ENDIAN_IS_LITTLE)
            add_compile_definitions(ENDIAN_IS_LITTLE=1)
        else()
            add_compile_definitions(ENDIAN_IS_LITTLE=0)
        endif()
    endif()
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel compilation
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS}")
    set(MS_NOWARN "/wd4018 /wd4068 /wd4101 /wd4244 /wd4267 /wd4305 /wd4715 /wd4996")
    set(MS_D "-D_CRT_SECURE_NO_WARNINGS -D_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${MS_NOWARN} ${MS_D}")
    set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
    if(BUILD_SHARED_LIBS)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    else()
        # required for post-build validation under vcpkg:
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd") # multithreaded, debug
    endif()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

elseif(CERF_COMPILE_OPTIONS)
    add_compile_options(${CERF_COMPILE_OPTIONS})

else()
    option(PEDANTIC "Compile with pedantic warnings" ON)
    option(WERROR "Treat warnings as errors" OFF)
    if(PEDANTIC)
        add_compile_options(-pedantic -Wall)
    endif()
    if(WERROR)
        add_compile_options(-Werror)
    endif()

    add_compile_options(-Wno-sign-compare -fno-omit-frame-pointer)
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        add_compile_options(-g)
    else()
        add_compile_options(-O3)
    endif()
endif()

include(CTest)

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

if(CERF_CPP)
    add_subdirectory(libcpp)
    add_subdirectory(testcpp)
endif()
if(CERF_C)
    add_subdirectory(libc)
    add_subdirectory(testc)
endif()
include("cmake/installCerf.cmake")

if(LIB_RUN)
    add_subdirectory(run)
endif()

if(LIB_MAN)
    add_subdirectory(man)
endif()
