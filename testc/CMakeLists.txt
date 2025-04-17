file(GLOB tests ../test/*test.c)

foreach(src ${tests})
    get_filename_component(test ${src} NAME_WE)
    add_executable(${test}_c ${src})
    target_include_directories(${test}_c PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
    target_link_libraries(${test}_c ${cerf_LIBRARY})
    add_test(NAME ${test}_c COMMAND ${test}_c)
endforeach()
