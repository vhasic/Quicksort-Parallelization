# CMakeList.txt : CMake project for PRS projekat, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
# Kompajliranje sa OpenMP
set(GCC_OPENMP_COMPILE_FLAGS "-fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_OPENMP_COMPILE_FLAGS}" )

#veličina steka; negdje je overridano ovo pa nema nikakvu ulogu
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--stack,500000000")
#set(CMAKE_CXX_STACK_SIZE "50000000")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -foo" )
#add_compile_options(/F 50000000)

# Add source to this project's executable.
add_executable (CMakeTarget "PRS projekat.cpp" )

# TODO: Add tests and install targets if needed.
