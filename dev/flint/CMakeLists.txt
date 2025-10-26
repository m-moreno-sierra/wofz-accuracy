cmake_minimum_required(VERSION 3.20)

project(cerf_flint CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)

find_package(cerf CONFIG REQUIRED COMPONENTS shared C)
get_property(Cerf_LIBRARIES TARGET cerf::cerf PROPERTY LOCATION)
message(STATUS "cerf: version = ${cerf_VERSION}, lib=${Cerf_LIBRARIES}")

add_executable(d2grid d2grid.cpp terms.cpp wn.cpp)
add_executable(runat runat.cpp terms.cpp wn.cpp)
add_executable(cover cover.cpp terms.cpp wn.cpp)
add_executable(const_arg const_arg.cpp terms.cpp wn.cpp)
add_executable(worst_per_tile worst_per_tile.cpp wn.cpp)
add_executable(statistics statistics.cpp relerr.cpp wn.cpp old_w_of_z.cpp)
#add_executable(wminat wminat.cpp terms.cpp wn.cpp)

target_link_libraries(d2grid PRIVATE flint OpenMP::OpenMP_CXX)
target_link_libraries(cover PRIVATE flint)
target_link_libraries(runat PRIVATE flint)
target_link_libraries(const_arg PRIVATE flint)
target_link_libraries(worst_per_tile PRIVATE flint ${Cerf_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(statistics PRIVATE flint ${Cerf_LIBRARIES} OpenMP::OpenMP_CXX)
#target_link_libraries(wminat PRIVATE flint)
