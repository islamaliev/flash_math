cmake_minimum_required(VERSION 2.8.4)
project(flash_math)

set(SOURCE_DIR src)

set(SOURCES
    "${SOURCE_DIR}/EulerAngles.cpp"
    "${SOURCE_DIR}/EulerAngles.h"
    "${SOURCE_DIR}/Mat4.cpp"
    "${SOURCE_DIR}/Mat4.h"
    "${SOURCE_DIR}/Orientation.cpp"
    "${SOURCE_DIR}/Orientation.h"
    "${SOURCE_DIR}/Quaternion.cpp"
    "${SOURCE_DIR}/Quaternion.h"
    "${SOURCE_DIR}/Vec4.cpp"
    "${SOURCE_DIR}/Vec4.h")

set(CMAKE_DEBUG_POSTFIX d)

add_library(${PROJECT_NAME} ${SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(${PROJECT_NAME} PUBLIC ${SOURCE_DIR})

function(install_package)
	set(PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/package)
	install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${PACKAGE_DIR}/lib ARCHIVE DESTINATION ${PACKAGE_DIR}/lib CONFIGURATIONS)
	install(DIRECTORY ${SOURCE_DIR}/ DESTINATION ${PACKAGE_DIR}/include FILES_MATCHING PATTERN "*.h")
	add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -D BUILD_TYPE=$<CONFIG> -P cmake_install.cmake)
endfunction()

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    install_package()
endif()
