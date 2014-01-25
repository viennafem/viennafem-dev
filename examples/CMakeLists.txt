set(EXAMPLES_TESTDATA
   data/cube1_hex.mesh
   data/cube8_hex.mesh
   data/cube343_hex.mesh
   data/cube3072.mesh
   data/line4.mesh
   data/line8.mesh
   data/square8.mesh
   data/square16_rect.mesh
   data/square128.mesh
   data/square224.mesh
   data/square512.mesh
   data/sshape2d.mesh
   data/sshape3d.mesh
   data/sshape3d-pimped.mesh
   )
foreach(f IN LISTS EXAMPLES_TESTDATA)
   configure_file(${f} "${CMAKE_CURRENT_BINARY_DIR}/${f}" COPYONLY)
endforeach()

add_subdirectory(tutorials)

