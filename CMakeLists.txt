cmake_minimum_required(VERSION 3.10)
project(main_1)

set(CMAKE_BUILD_TYPE "RELEASE")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall ")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")  
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")	
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")	

aux_source_directory(./src DIR_SRCS)

add_executable(main_1  ./main_test/main_1.cpp  ${DIR_SRCS})

target_link_libraries(main_1 lapacke blas)

target_include_directories(main_1 PUBLIC ${PROJECT_SOURCE_DIR}/include)

project(main_2)

add_executable(main_2  ./main_test/main_2.cpp  ${DIR_SRCS})

target_link_libraries(main_2 lapacke blas)

target_include_directories(main_2 PUBLIC ${PROJECT_SOURCE_DIR}/include)

project(Test_1)

add_executable(Test_1  ./main_test/Test_1.cpp  ${DIR_SRCS})

target_link_libraries(Test_1 lapacke blas)

target_include_directories(Test_1 PUBLIC ${PROJECT_SOURCE_DIR}/include)

project(Test_2)

add_executable(Test_2  ./main_test/Test_2.cpp  ${DIR_SRCS})

target_link_libraries(Test_2 lapacke blas)

target_include_directories(Test_2 PUBLIC ${PROJECT_SOURCE_DIR}/include)

