include_directories(
  SYSTEM
  )
  

# ~~ Programs

if (hdf5 AND mpi)

  # ~~ hdf5 programs

  set( program_sources
    exciton_photo.cpp
  )

  foreach(source_file ${program_sources})

    STRING(REGEX REPLACE ".cpp" ".ex" executable_file "${source_file}")  
    add_executable(${executable_file} ${source_file})
    target_link_libraries(${executable_file} cntr cntr_hdf5 ${COMMON_LIBS})
    install(TARGETS ${executable_file} DESTINATION bin)

  endforeach(source_file)

endif (hdf5 AND mpi)