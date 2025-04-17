set(apps
    run_erfcx
    run_fwhm
    run_imwofx
    run_voigt
    run_wofz
    tabulate
    bigloop
    )

if (CERF_C)
    foreach(app ${apps})
        add_executable(${app} ${app}.c)
        target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
        target_link_libraries(${app} ${cerf_LIBRARY})
    endforeach()

else()
    add_compile_definitions("CERF_AS_CPP=ON")
    foreach(app ${apps})
        set_property(SOURCE ${app}.c PROPERTY LANGUAGE CXX)
        add_executable(${app} ${app}.c)
        target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
        set_property(TARGET ${app} PROPERTY LINKER_LANGUAGE CXX)
        target_link_libraries(${app} ${cerfcpp_LIBRARY})
    endforeach()

endif()
