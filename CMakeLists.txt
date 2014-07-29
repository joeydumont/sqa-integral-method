# Name of project
project(2d-sword)
set (2d-sword_VERSION_MAJOR 0)
set (2d-sword_VERSION_MINOR 1)
set (2d-sword_VERSION_RELEASE 1)

# CMake config
cmake_minimum_required(VERSION 2.8)
set (CMAKE_INSTALL_PREFIX /usr/)

# Compiler config
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -march=native")

# Included files
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Source files
aux_source_directory(./src SRC_LIST )

# Build a shared library
add_library(${PROJECT_NAME} SHARED ${SRC_LIST}) 
TARGET_LINK_LIBRARIES(${PROJECT_NAME} armadillo)

SET_TARGET_PROPERTIES(${PROJECT_NAME}
    PROPERTIES
        VERSION ${2d-sword_VERSION_MAJOR}.${2d-sword_VERSION_MINOR}.${2d-sword_VERSION_RELEASE}
        SOVERSION ${2d-sword_VERSION_MAJOR}.${2d-sword_VERSION_MINOR}.${2d-sword_VERSION_RELEASE})

# Install directories
install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
