cmake_minimum_required(VERSION 3.10)

project(2dv)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED)
message(STATUS "VTK_DIR: ${VTK_DIR}")

find_package(SDL2 REQUIRED)

include_directories(
    ${OpenCV_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/wlib/include
    ${PROJECT_SOURCE_DIR}/include
)

add_subdirectory(${PROJECT_SOURCE_DIR}/wlib)

include(${VTK_USE_FILE})

add_executable(
    ${PROJECT_NAME}
    src/main.cpp
    src/parameter2dv.cpp
    src/force2dv.cpp
    src/division.cpp
    
    src/geo2dv.cpp
    src/wangMath.cpp
    src/IO2dv.cpp
    src/wangSystem.cpp
    src/Initial2dv.cpp
    src/Analysis2dv.cpp

    src/Organ2dv.cpp
    src/Cell2dv.cpp
    src/Line2dv.cpp
    src/Vertex2dv.cpp

    src/gnuplot.cpp
    src/vtk2dv.cpp
    src/image2dv.cpp
)

target_link_libraries(2dv wLib ${OpenCV_LIBS}   ${VTK_LIBRARIES} ${SDL2_LIBRARIES})

# Include the SDL2 include directories
target_include_directories(2dv PRIVATE ${SDL2_INCLUDE_DIRS})