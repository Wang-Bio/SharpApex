cmake_minimum_required(VERSION 3.10) # Ensure compatibility

# Define the library
add_library(wLib src/wSystem.cpp)

# Specify where the library's headers are
target_include_directories(wLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)