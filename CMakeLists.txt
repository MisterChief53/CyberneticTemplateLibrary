cmake_minimum_required(VERSION 3.10)

# Set project name
project(CyberneticTemplateLibrary)

# Set the source directories for different parts of the project
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add an executable target
add_executable(${PROJECT_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

# Specify the include directories for the target
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIR})