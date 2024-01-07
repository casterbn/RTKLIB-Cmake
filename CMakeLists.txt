cmake_minimum_required(VERSION 3.0)
project(RTKLIB-CMake)

# set build flags. 
set(CMAKE_CXX_FLAGS "-std=c++20" )
set(CMAKE_CXX_FLAGS "-fpermissive")
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -fsee -fomit-frame-pointer -fno-signed-zeros -fno-math-errno -funroll-loops")
endif()

# set output dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/app/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/app/lib)

# RTKLIB
add_definitions(-DENAGLO -DENACMP -DENAGAL -DNFREQ=3 -DNEXOBS=3 -DDLL)
add_subdirectory(rtklib)

# VT
add_library(vt vt/vt.c)
target_link_libraries(vt rtklib)
target_include_directories(vt PUBLIC vt)

# executable
add_executable(RNX2RTKP app/rnx2rtkp.c)
target_link_libraries(RNX2RTKP rtklib)

add_executable(RTKRCV app/rtkrcv.c)
target_link_libraries(RTKRCV rtklib vt)

add_executable(STR2STR app/str2str.c)
target_link_libraries(STR2STR rtklib)

add_executable(CONVBIN app/convbin.c)
target_link_libraries(CONVBIN rtklib)

add_executable(POS2KML app/pos2kml.c)
target_link_libraries(POS2KML rtklib)