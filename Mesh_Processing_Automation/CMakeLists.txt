cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(SIMPLIFY_MESH)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF) # Disable GNU extensions (e.g., std=gnu++14).
set(CMAKE_BUILD_TYPE Release) #bb added

# Eigen.
find_package(PkgConfig REQUIRED)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories( ${EIGEN3_INCLUDE_DIRS} )

# PCL.
find_package(PCL 1.8 REQUIRED COMPONENTS common io visualization)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# # LibIGL.
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# option(LIBIGL_WITH_OPENGL      "Use OpenGL" ON)
# option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"   ON)
# option(LIBIGL_WITH_CGAL        "Use CGAL"   ON) #bbb added
# find_package(LIBIGL REQUIRED QUIET)

# CGAL
find_package(CGAL)

# Build the executable.
file(GLOB SRCFILES src/*.cpp)
# This should set everything in src/ to the SRCFILES
add_executable(simplify_mesh ${SRCFILES})
target_link_libraries(simplify_mesh CGAL::CGAL)
# target_link_libraries(simplify_mesh igl::core igl::opengl igl::opengl_glfw igl::cgal)
target_link_libraries(simplify_mesh ${PCL_LIBRARIES})
# igl_copy_cgal_dll(${PROJECT_NAME})
