cmake_minimum_required(VERSION 2.8.4)
project(flash_math)

set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
        CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

set(INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS})

include_directories(${INCLUDE_DIRS})

set(SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/EulerAngles.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Mat4.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Orientation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Quaternion.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Vec4.cpp"
    "${INCLUDE_DIRS}/EulerAngles.h"
    "${INCLUDE_DIRS}/Mat4.h"
    "${INCLUDE_DIRS}/Orientation.h"
    "${INCLUDE_DIRS}/Quaternion.h"
    "${INCLUDE_DIRS}/Vec4.h")

add_library(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
