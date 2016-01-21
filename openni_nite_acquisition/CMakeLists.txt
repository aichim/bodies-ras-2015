cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project (nite2_body_tracker)


add_definitions (-std=c++11)


set(OPENNI2_DEFINITIONS ${PC_OPENNI_CFLAGS_OTHER})
FIND_LIBRARY( OPENNI2_LIBRARY
             NAMES OpenNI2
             HINTS ${PC_OPENNI2_LIBDIR} ${PC_OPENNI2_LIBRARY_DIRS} /usr/lib /home/alex/Work/Libs/OpenNI2/Bin/x64-Release/
             PATHS "$ENV{PROGRAMFILES}/OpenNI2/Lib${OPENNI2_SUFFIX}" "$ENV{PROGRAMW6432}/OpenNI2/ Lib${OPENNI2_SUFFIX}" "$ENV{PROGRAMW6432}/OpenNI2" "/home/alex/Work/Libs/OpenNI2/Lib${OPENNI2_SUFFIX}"
             PATH_SUFFIXES lib lib64
)
FIND_PATH( OPENNI2_INCLUDE_DIR OpenNI.h
          HINTS ${PC_OPENNI2_INCLUDEDIR} ${PC_OPENNI2_INCLUDE_DIRS} /usr/include/openni2 /usr/include/ni2 /home/alex/Work/Libs/OpenNI2/Include
          PATHS "$ENV{PROGRAMFILES}/OpenNI2/include" "$ENV{PROGRAMW6432}/OpenNI2/include" "/home/alex/Work/Libs/OpenNI2/Include"
          PATH_SUFFIXES openni2 ni2)

FIND_LIBRARY( NITE2_LIBRARY
             NAMES NiTE2
             HINTS ${PC_OPENNI2_LIBDIR} ${PC_OPENNI2_LIBRARY_DIRS} /usr/lib /home/alex/Work/Libs/NiTE-Linux-x64-2.2/Redist
             PATHS "$ENV{PROGRAMFILES}/PrimeSense/NiTE2/lib${OPENNI2_SUFFIX}" "$ENV{PROGRAMW6432}/PrimeSense/NiTE2/lib${OPENNI2_SUFFIX}" "/home/alex/Work/Libs/NiTE-Linux-x64-2.2/lib${OPENNI2_SUFFIX}"
             PATH_SUFFIXES lib
)
FIND_PATH( NITE2_INCLUDE_DIR NiTE.h
          HINTS ${PC_OPENNI2_INCLUDEDIR} ${PC_OPENNI2_INCLUDE_DIRS} /usr/include/openni2 /usr/include/nite2 /home/alex/Work/Libs/NiTE-Linux-x64-2.2/Include
          PATHS "$ENV{PROGRAMFILES}/PrimeSense/NiTE2/include" "$ENV{PROGRAMW6432}/PrimeSense/NiTE2/include" "/home/alex/Work/Libs/NiTE-Linux-x64-2.2/Include"
          PATH_SUFFIXES openni2         
)

message ("openni library: " ${OPENNI2_LIBRARY})
message ("nite 2 library: " ${NITE2_LIBRARY})
message ("nite2 include dir: " ${NITE2_INCLUDE_DIR})
include_directories(${NITE2_INCLUDE_DIR})
include_directories(${OPENNI2_INCLUDE_DIR})


FIND_PACKAGE (OpenCV REQUIRED)


ADD_EXECUTABLE (nite2_body_tracker main.cpp)
TARGET_LINK_LIBRARIES (nite2_body_tracker ${OPENNI2_LIBRARY} ${NITE2_LIBRARY} ${OpenCV_LIBS})