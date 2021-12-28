CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(Gsolver)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" -std=c++11")
include_directories(cmake-build-debug/armadillo-10.7.5/include)
link_directories(/cmake-build-debug/armadillo-10.7.5/examples/lib_win64)
add_executable(Cyliner2D main.cpp)
#target_link_libraries(Cylinder2D /cmake-build-debug/armadillo-10.7.5/examples/lib_win64/libopenblas.lib)