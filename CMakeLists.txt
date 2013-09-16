cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(CoGL)

set (USE_CUDA true CACHE BOOL "Use CUDA")
set (ENABLE_INTEROP ON CACHE BOOL "Use interop")
set (PISTON_DIR "" CACHE PATH "PISTON directory")

if(APPLE)
  include_directories(/System/Library/Frameworks)
  find_library(OPENGL_LIBRARIES OpenGL)
  find_library(GLUT_LIBRARIES GLUT)
  mark_as_advanced(OPENGL_LIBRARIES GLUT_LIBRARIES)
else ()
  find_package(OpenGL)
  find_package(GLUT)
endif ()

if (USE_CUDA)
  find_package(CUDA)
  include_directories(${CMAKE_SOURCE_DIR})
  include_directories(${CUDA_TOOLKIT_ROOT_DIR}/include)
else ()
  set (THRUST_DIR "" CACHE PATH "Thrust directory")
  include_directories(${CMAKE_SOURCE_DIR})
  include_directories(${THRUST_DIR})
endif ()

include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PISTON_DIR})

set (DATA_DIR "" CACHE PATH "Data directory")

if (ENABLE_INTEROP)
  set(GPU_FLAGS "-DUSE_INTEROP")
endif ()

if (USE_CUDA)
  cuda_add_executable(coglGPU main.cu CoGLRender.cu ${PISTON_DIR}/piston/util/quaternion.cpp OPTIONS "-DDATA_DIRECTORY=${DATA_DIR}" ${GPU_FLAGS})
  target_link_libraries(coglGPU pthread GLEW ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
endif ()

add_executable(coglOMP main.cpp CoGLRender.cpp ${PISTON_DIR}/piston/util/quaternion.cpp)
set_target_properties(coglOMP PROPERTIES COMPILE_FLAGS "-fopenmp -DDATA_DIRECTORY=${DATA_DIR} -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP")
target_link_libraries(coglOMP pthread gomp GLEW ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
