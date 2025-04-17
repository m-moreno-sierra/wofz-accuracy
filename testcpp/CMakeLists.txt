file(GLOB tests ../test/*test.c)

add_compile_definitions("CERF_AS_CPP=ON")

foreach(src ${tests})
    set_property(SOURCE ${src} PROPERTY LANGUAGE CXX)
    get_filename_component(test ${src} NAME_WE)
    add_executable(${test}_cpp ${src})
    target_include_directories(${test}_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
    set_property(TARGET ${test}_cpp PROPERTY LINKER_LANGUAGE CXX)
    target_link_libraries(${test}_cpp ${cerfcpp_LIBRARY})
    add_test(NAME ${test}_cpp COMMAND ${test}_cpp)
endforeach()
