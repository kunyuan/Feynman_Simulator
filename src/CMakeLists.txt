cmake_minimum_required (VERSION 2.6)
project (FeynmanSimulator)
set(CMAKE_CXX_COMPILER icc)
#set(CMAKE_CXX_FLAGS "-lz")
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

file(GLOB_RECURSE SRCS *.cpp)
file(GLOB_RECURSE HDRS *.h)
ADD_EXECUTABLE(FeynmanSimulator  ${SRCS} ${HDRS})

# $ZLIB_DIR is an environment variable that would
# correspond to the ./configure --prefix=$ZLIB_DIR
# used in building zlib.
FIND_PATH(ZLIB_INCLUDE_DIR zlib.h
    $ENV{ZLIB_DIR}/include
    $ENV{ZLIB_DIR}
    /usr/local/include
    /usr/include
    /opt/local/include # DarwinPorts
    /opt/include
)

FIND_LIBRARY(ZLIB_LIBRARY
    NAMES z libz zlib
    PATHS
    $ENV{ZLIB_DIR}/lib
    $ENV{ZLIB_DIR}
    /usr/local/lib
    /usr/lib
    /opt/local/lib
    /opt/lib
)

message(${ZLIB_INCLUDE_DIR})
message(${ZLIB_LIBRARY})
IF(ZLIB_LIBRARY)
    message("zlib is found.")
    target_link_libraries(FeynmanSimulator ${ZLIB_LIBRARY})
ELSE()
    message("zlib is not found!")
ENDIF()

