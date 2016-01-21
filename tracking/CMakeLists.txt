cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project (BODIES_RAS_2015)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_definitions (-std=c++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

find_package (Boost REQUIRED filesystem system)
link_directories (${Boost_LIBRARY_DIRS})

find_package (Eigen REQUIRED)
include_directories (${EIGEN_INCLUDE_DIR})

find_package (PCL 1.7 REQUIRED common io search registration surface visualization)
include_directories (${PCL_INCLUDE_DIRS})
link_directories (${PCL_LIBRARY_DIRS})
add_definitions (${PCL_DEFINITIONS})

find_package (OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

if ( "${OpenCV_VERSION}" VERSION_GREATER 2.10 )
	add_definitions( -DOPENCV_3_0 )
endif  ( "${OpenCV_VERSION}" VERSION_GREATER 2.10 )

find_package (FBX REQUIRED)
include_directories (${FBX_INCLUDE_DIR})

find_package (GLUT REQUIRED)
find_package (GLEW REQUIRED)
find_package (QGLViewer REQUIRED)

FIND_PACKAGE (Qt4 COMPONENTS QtCore QtGui QtOpenGL QtXml QtNetwork REQUIRED)
INCLUDE (${QT_USE_FILE})
ADD_DEFINITIONS( ${QT_DEFINITIONS})
SET (QT_USE_QTOPENGL TRUE)



include_directories (include)


set (INCS include/common.h
		      include/mesh.h
		      include/blendshape_model.h
		      include/bodies_utils.h
		      include/skeleton.h
          include/skeleton_mesh.h
          include/SkeletonIOFbx.h
          include/contours.h
          include/visualization_utils.h)
set (IMPLS include/impl/mesh.hpp
		   include/impl/blendshape_model.hpp)
set (SRCS src/common.cpp
		      src/mesh.cpp
		      src/blendshape_model.cpp
		      src/skeleton.cpp
          src/skeleton_mesh.cpp
          src/bodies_utils.cpp
          src/contours.cpp
          src/visualization_utils.cpp)

SET (libs_to_link ${QT_LIBRARIES} ${QGLVIEWER_LIBRARY} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${FBX_LIBRARY})


add_library (af_bodies_lib ${INCS} ${IMPLS} ${SRCS})
target_link_libraries (af_bodies_lib ${libs_to_link})

add_executable (build_body_pose_pca_model_eigen src/main_build_body_pose_pca_model_eigen.cpp)
target_link_libraries (build_body_pose_pca_model_eigen af_bodies_lib ${libs_to_link})

add_executable (tracking_modeling_offline src/main_tracking_modeling_offline.cpp)
target_link_libraries (tracking_modeling_offline af_bodies_lib ${libs_to_link})

add_executable (tracking_modeling_online src/main_tracking_modeling_online.cpp)
target_link_libraries (tracking_modeling_online af_bodies_lib ${libs_to_link})

add_executable (process_results src/main_process_results.cpp)
target_link_libraries (process_results af_bodies_lib ${libs_to_link})
