include_directories(
  SYSTEM
  )
  

# ~~ Programs

if (hdf5 AND mpi)

  # ~~ hdf5 programs

  set( program_sources
    #lattice_rpa_1d_semiconductor.cpp
    #lattice_rpa_1d_semiconductor_optical.cpp
    lattice_rpa_1d_semiconductor_optical_abinitio.cpp
    #lattice_rpa_1d_semiconductorDEB.cpp
    #test_hdf5.cpp
    #lattice_rpa_1d_semimetal.cpp
    #lattice_hubbard.cpp
    #fkm_bethe_quench.cpp
    #fkm_bethe_floquet.cpp
    #test_adam.cpp
    #rotate.cpp
    #lattice_rpa_1d_holstein.cpp
  )

  foreach(source_file ${program_sources})

    STRING(REGEX REPLACE ".cpp" ".ex" executable_file "${source_file}")  
    add_executable(${executable_file} ${source_file})
    target_link_libraries(${executable_file} cntr cntr_hdf5 ${COMMON_LIBS})
    install(TARGETS ${executable_file} DESTINATION bin)

  endforeach(source_file)

endif (hdf5 AND mpi)
