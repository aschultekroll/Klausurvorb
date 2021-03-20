
cmake_minimum_required(VERSION 3.1.3 FATAL_ERROR)
project(Klausur VERSION 1.0 LANGUAGES CXX)
#set(DCMAKE_CXX_FLAGS="-O2") #sets Compiler flags

##Find Packages####


find_package(OpenMP)
find_package(Boost 1.56.0 REQUIRED)
find_package(MPI REQUIRED)

find_package(Boost COMPONENTS mpi)
##Add Executables###

add_executable(Main main.cpp)
add_executable(Main2 main2.cpp)
##Add Librarys###
 #Bsp: add_library(merge mergesort2.cpp)

##Link Librarys###



target_link_libraries(Main2 MPI::MPI_CXX)
#target_link_libraries(Main Boost::mpi)

target_link_libraries(Main2 Boost::boost OpenMP::OpenMP_CXX)
target_link_libraries(Main MPI::MPI_CXX)
#target_link_libraries(Main Boost::mpi)

target_link_libraries(Main Boost::boost OpenMP::OpenMP_CXX) #links other librarys to library
#target_link_libraries(Main library) #link librarys to executable
#target_link_libraries(Test library)

##Add Tests###
add_test(NAME DER_TEST COMMAND Test)
