cmake_minimum_required(VERSION 3.10)

project(d3aist CXX)

SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

# This is to prevent a warning due to boost libraries in CMake
if(POLICY CMP0167) # Use upstream's FindBoost.cmake to find Boost libraries.
  cmake_policy(SET CMP0167 NEW) # Use upstream's FindBoost.cmake to find Boost libraries.
endif()

# This should not be needed, but CMake could not find my ITKConfig.cmake file unless the path is explictly stated
# If FIND_PACKAGE(ITK) fails in your system, set the correct directory for it and uncomment next line
# SET(CMAKE_PREFIX_PATH "/usr/lib64/cmake/InsightToolkit")
FIND_PACKAGE(ITK)

# This is for my particular installation directories. Change it to fit your system.
# If not set, it defaults to /usr/local
# SET(CMAKE_INSTALL_PREFIX "/home/RAID5" CACHE STRING "Directory for local instalation" FORCE)

ADD_DEFINITIONS(-Wall -Wno-dev -O2)

include(${ITK_USE_FILE})
add_executable(csv2bin csv2bin.cpp auxbin.cpp ProcessCSV.cpp Gene.cpp Genepanel.cpp Cell.cpp Cellset.cpp Transcript.cpp Transcriptset.cpp)
add_executable(buildim buildim.cpp auxbin.cpp ProcessCSV.cpp Cell.cpp Cellset.cpp Gene.cpp Genepanel.cpp Transcript.cpp Transcriptset.cpp SpPseudoImage.cpp SpRData.cpp)
add_executable(imfuse imfuse.cpp)

TARGET_LINK_LIBRARIES(buildim ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(imfuse ${ITK_LIBRARIES})

INSTALL(TARGETS csv2bin buildim imfuse RUNTIME DESTINATION bin)

