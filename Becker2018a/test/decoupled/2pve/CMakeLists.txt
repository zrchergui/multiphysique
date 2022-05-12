add_input_file_links()

add_dumux_test(test_2pve test_2pve test_2pve.cc)

install(FILES
        mystart.hh
        test_2pve.cc
        test_2pveproblem.hh
        test_2pvespatialparams.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/test/decoupled/2pve)
