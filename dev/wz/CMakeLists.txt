##  ************************************************************************************************
##
##  Project:   Code generating code for Taylor expansion on square tiles
##
##  Reference: Joachim Wuttke,
##             Code generation for computing an analytical function with near machine precision
##             on square tiles, with application to the Faddeeva function
##             (in preparation, available upon request)
##
##  File:      CMakeLists.txt
##
##  License:   MIT License (file LICENSE in toplevel directory)
##  Copyright: Forschungszentrum Jülich GmbH 2025
##  Author:    Joachim Wuttke <j.wuttke@fz-juelich.de>
##
##  ************************************************************************************************

cmake_minimum_required(VERSION 3.20)

project(cerf_flint CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)

find_package(cerf CONFIG REQUIRED COMPONENTS shared CXX)
get_property(Cerf_LIBRARIES TARGET cerf::cerfcpp PROPERTY LOCATION)
message(STATUS "cerf: version = ${cerf_VERSION}, lib=${Cerf_LIBRARIES}")

add_executable(ref_val ref_val.cpp wn.cpp hp_fref.cpp)
add_executable(d2grid d2grid.cpp terms.cpp wn.cpp hp_fref.cpp)
add_executable(runat runat.cpp terms.cpp wn.cpp hp_fref.cpp)
add_executable(cover cover.cpp terms.cpp wn.cpp hp_fref.cpp)
add_executable(worst_per_tile worst_per_tile.cpp wn.cpp hp_fref.cpp)
add_executable(statistics statistics.cpp relerr.cpp wn.cpp hp_fref.cpp)
#add_executable(wminat wminat.cpp terms.cpp wn.cpp hp_fref.cpp)

target_link_libraries(ref_val PRIVATE flint)
target_link_libraries(d2grid PRIVATE flint OpenMP::OpenMP_CXX)
target_link_libraries(cover PRIVATE flint)
target_link_libraries(runat PRIVATE flint)
target_link_libraries(worst_per_tile PRIVATE flint ${Cerf_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(statistics PRIVATE flint ${Cerf_LIBRARIES} OpenMP::OpenMP_CXX)
#target_link_libraries(wminat PRIVATE flint)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)
