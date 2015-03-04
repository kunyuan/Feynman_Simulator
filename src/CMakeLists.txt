cmake_minimum_required (VERSION 2.6)
project (FeynmanSimulator)
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_CXX_COMPILER clang++)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(${FeynmanSimulator_SOURCE_DIR})

file(GLOB_RECURSE SRCS *.cpp)
file(GLOB_RECURSE HDRS *.h)
ADD_EXECUTABLE(simulator.exe  ${SRCS} ${HDRS})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(NumPy REQUIRED)
find_package(PythonLibsNew REQUIRED)
message("python include dir:" ${PYTHON_INCLUDE_DIR})
message("python libraries:" ${PYTHON_LIBRARIES})
message("numpy version:" ${NUMPY_VERSION})
message("numpy include dir:" ${NUMPY_INCLUDE_DIRS})

# Require version >= 1.6
#if(NUMPY_VERSION_DECIMAL LESS 10600)
    #message(FATAL_ERROR,
        #"requires NumPy >= 1.6")
#endif()

include_directories(
    ${PYTHON_INCLUDE_DIR}
    ${NUMPY_INCLUDE_DIRS}
    )
target_link_libraries(simulator.exe ${PYTHON_LIBRARIES})
install (TARGETS simulator.exe DESTINATION ${PROJECT_SOURCE_DIR}/..)
