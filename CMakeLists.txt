cmake_minimum_required(VERSION 3.10)
project(IVPs)

set(CMAKE_BUILD_TYPE "DEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall ")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")  
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")	
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")	

aux_source_directory(./src DIR_SRCS)

add_executable(IVPs  main.cpp  ${DIR_SRCS})

target_link_libraries(IVPs lapacke blas)

target_include_directories(IVPs PUBLIC ${PROJECT_SOURCE_DIR}/include)

project(Test1)

add_executable(Test1  Test1.cpp  ${DIR_SRCS})

target_link_libraries(Test1 lapacke blas)

target_include_directories(Test1 PUBLIC ${PROJECT_SOURCE_DIR}/include)
