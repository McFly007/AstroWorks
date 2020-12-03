cmake_minimum_required (VERSION 3.11)

project (fastbvh)

add_library (fastbvh)
target_sources (fastbvh PRIVATE src/BBox.cpp src/BVH.cpp)
target_include_directories (fastbvh PUBLIC include)
target_include_directories (fastbvh PRIVATE src)

find_package (pybind11 QUIET)
if (${pybind11_FOUND})
  pybind11_add_module (pybvh MODULE src/pybvh.cpp)
  target_link_libraries (pybvh PRIVATE fastbvh)
  target_include_directories (pybvh PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endif ()
