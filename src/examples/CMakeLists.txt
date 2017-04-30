add_executable(srrg_sdf_map_example srrg_sdf_map_example.cpp)

target_link_libraries(srrg_sdf_map_example
  srrg_system_utils_library
  srrg_txt_io_library
  srrg_sdf_map_library
  ${catkin_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARY}
)
