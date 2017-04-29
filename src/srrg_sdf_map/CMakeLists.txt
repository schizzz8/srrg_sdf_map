add_library(srrg_sdf_map_library SHARED
	sdf_map.h sdf_map.cpp
)

target_link_libraries(srrg_sdf_map_library
	${catkin_LIBRARIES}
)
