add_input_file_links()

add_dumux_test(test_dec2p test_dec2p test_dec2p.cc)

install(FILES
        test_dec2p.cc
        test_dec2pproblem.hh
        test_dec2pspatialparams.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/appl/energyStorage)
