cmake_minimum_required(VERSION 3.10)

project(2dv)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(
    ${PROJECT_NAME}
    src/main.cpp
    src/wFile.cpp
    src/wContour.cpp
    src/wMode.cpp
    src/wBMP.cpp
    src/wAnalyzer.cpp
)