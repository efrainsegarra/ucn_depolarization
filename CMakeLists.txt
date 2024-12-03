# Set up
cmake_minimum_required(VERSION 3.11)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
set(CMAKE_VERBOSE_MAKEFILE ON)

#find_package(ROOT REQUIRED)
find_package(CMinpack REQUIRED)
project (gradient_fit)


set(INCLUDES ${CMAKE_SOURCE_DIR}/include /xdata/CMINPACK/include/cminpack-1 )
include_directories(${INCLUDES})


add_executable(gsl_trial gsl_trial.cpp)
target_link_libraries(gsl_trial cminpack::cminpack )
