add_input_file_links()

add_dumux_test(test_multidim test_multidim test_multidim.cc)

install(FILES
        test_multidim.cc
        test_multidimproblem.hh
        test_multidimspatialparams.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test/modelcoupling/2pve2pfullmono)
