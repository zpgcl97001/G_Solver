CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
project(Gsolver)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" -std=c++11")
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(cmake-build-debug/armadillo-10.7.5/include)
link_directories(/cmake-build-debug/armadillo-10.7.5/examples/lib_win64)
add_executable(Cylinder2D main.cpp)
target_link_libraries(Cylinder2D ${LAPACK_LIBRARIES} ${BLAS_LIBARIES} ${ARMADILLO_LIBRARIES})
