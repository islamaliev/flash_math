cmake_minimum_required(VERSION 2.8.4)
project(flash_math)

set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
        CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

include_directories(${${PROJECT_NAME}_INCLUDE_DIRS})

set(TARGET_NAME fmath)

set(${PROJECT_NAME}_TARGET_NAME ${TARGET_NAME}
        CACHE INTERNAL "${PROJECT_NAME}: Lib Binary Name" FORCE)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(${TARGET_NAME} ${SOURCES})
set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE CXX)
