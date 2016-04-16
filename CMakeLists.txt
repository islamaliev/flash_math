cmake_minimum_required(VERSION 2.8.4)
project(flash_math)

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(${PROJECT_NAME}_INCLUDE_DIRS ${SOURCE_DIR}/core ${SOURCE_DIR}/orientation
        CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

set(TARGET_NAME fmath)

set(${PROJECT_NAME}_TARGET_NAME ${TARGET_NAME}
        CACHE INTERNAL "${PROJECT_NAME}: Lib Binary Name" FORCE)

file(GLOB SOURCES ${SOURCE_DIR}/**/*.cpp)
add_library(${TARGET_NAME} ${SOURCES})
set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE CXX)