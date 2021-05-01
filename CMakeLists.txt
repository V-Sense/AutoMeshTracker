CMAKE_MINIMUM_REQUIRED(VERSION 3.4.1)


SET(PROJECTNAME "MeshTracker")
PROJECT(${PROJECTNAME} C CXX)

set (CMAKE_BUILD_TYPE "Release")

set (CMAKE_CXX_STANDARD 17)

include_directories( ${CMAKE_CURRENT_BINARY_DIR} )

# Add CGAL
FIND_PACKAGE(CGAL REQUIRED)
IF (CGAL_FOUND)
  message(CGAL_FOUND)
  include_directories(${CGAL_INCLUDE_DIRS})
  link_directories(${CGAL_LIBRARY_DIRS})
  include(${CGAL_USE_FILE})
ENDIF (CGAL_FOUND)

# Add Eigen3
FIND_PACKAGE(Eigen3 REQUIRED)
IF (EIGEN3_FOUND)
  message(EIGEN3_FOUND)
	INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
ENDIF (EIGEN3_FOUND)

# Add SuiteSparse
find_package(SuiteSparse REQUIRED COMPONENTS CHOLMOD)
if (SUITESPARSE_FOUND)
  message(SUITESPARSE_FOUND)
  add_library (SuiteSparse::Cholmod INTERFACE IMPORTED)
  target_link_libraries (SuiteSparse::Cholmod INTERFACE ${SuiteSparse_CHOLMOD_LIBRARY_RELEASE})
  target_include_directories (SuiteSparse::Cholmod INTERFACE ${SuiteSparse_INCLUDE_DIRS})
ENDIF(SUITESPARSE_FOUND)

# Add OpenCV
find_package(OpenCV REQUIRED)
if (OPENCV_FOUND)
  message(OPENCV_FOUND)
  INCLUDE_DIRECTORIES(${OPENCV_INCLUDE_DIRS})
ENDIF (OPENCV_FOUND)

# Add Point Cloud Library
find_package(PCL 1.8 REQUIRED)
if (PCL_FOUND)
  message(PCL_FOUND)
  list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")
  add_definitions(-DPCL_NEW=1)
  include_directories(${PCL_INCLUDE_DIRS})
  link_directories(${PCL_LIBRARY_DIRS})
  add_definitions(${PCL_DEFINITIONS})
endif (PCL_FOUND)

# Add Boost
FIND_PACKAGE(Boost COMPONENTS program_options log REQUIRED)
IF (Boost_FOUND)
  message(BOOST_FOUND)
  ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF (Boost_FOUND)

#Add VTK
FIND_PACKAGE(VTK REQUIRED)
IF (VTK_FOUND)
  message(VTK_FOUND)
  INCLUDE (${VTK_USE_FILE})
ENDIF (VTK_FOUND)

find_package(OpenMP)
if (OPENMP_FOUND)
  message(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(meshtracker)

install(TARGETS ${PROJECTNAME} DESTINATION bin)