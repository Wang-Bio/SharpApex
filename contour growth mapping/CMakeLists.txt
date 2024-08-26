cmake_minimum_required(VERSION 3.10)

project(2dv)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenCV REQUIRED)
include_directories(
    ${OpenCV_INCLUDE_DIRS}
)


add_executable(
    ${PROJECT_NAME}
    src/main.cpp
    src/plot.cpp
    src/initialization.cpp
    src/growth.cpp
    src/geo.cpp
)

target_link_libraries(2dv ${OpenCV_LIBS})