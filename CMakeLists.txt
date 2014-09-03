cmake_minimum_required(VERSION 2.8)

project(DislScatANAThreadingMisfit)

#Check if important libraries are available.

#Gnu Scientific library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules") # findGSL module is not shipped with cmake
find_package(GSL REQUIRED)

#LAPACK 
find_package(LAPACK REQUIRED)

#BLAS
find_package(BLAS REQUIRED)

#some libraries are build as external projects
include(ExternalProject)

#config library
set(config_project_prefix "${CMAKE_BINARY_DIR}")
ExternalProject_Add(
    config_project
    PREFIX ${config_project_prefix}
    URL http://www.hyperrealm.com/libconfig/libconfig-1.4.9.tar.gz
    SOURCE_DIR ${config_project_prefix}/src/config_project
    CONFIGURE_COMMAND ${config_project_prefix}/src/config_project/configure --prefix=${config_project_prefix}/src/config_project-build
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
)

#misfit interfaces library
ExternalProject_Add(
    minterface_project
    PREFIX ${CMAKE_BINARY_DIR}
    GIT_REPOSITORY https://github.com/vifactor/AnalyticalMisfitInterface.git
    GIT_TAG origin
)

#compile project itself
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src/minterface_project-build)
include_directories(${CMAKE_BINARY_DIR}/src/config_project-build/include)

link_directories(${CMAKE_BINARY_DIR}/src/minterface_project-build)
link_directories(${CMAKE_BINARY_DIR}/src/config_project-build/lib)

file(GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
set(EXEC_NAME DislScatCubFITMisfitThreading)
add_executable(${EXEC_NAME} ${CPP_SOURCES})

add_dependencies(${EXEC_NAME} config_project)
add_dependencies(${EXEC_NAME} minterface_project)

target_link_libraries(${EXEC_NAME} interface config++  gsl blas lapack)
