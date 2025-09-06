cmake_minimum_required(VERSION 3.20)
project(hpsx VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific options
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/W4)
endif()

# Create header from source
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/hpsx/src.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/hpsx.h
    COPYONLY
)

# Header-only library target
add_library(hpsx INTERFACE)
add_library(hpsx::hpsx ALIAS hpsx)

target_include_directories(hpsx
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(hpsx INTERFACE cxx_std_20)

# Example executable
add_executable(hpsx_example ${CMAKE_CURRENT_SOURCE_DIR}/hpsx/src.cpp)
target_link_libraries(hpsx_example PRIVATE hpsx::hpsx)

# Install rules
include(GNUInstallDirs)

install(TARGETS hpsx
        EXPORT hpsxTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hpsx.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT hpsxTargets
        FILE hpsxTargets.cmake
        NAMESPACE hpsx::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hpsx)

# Create config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    hpsxConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/hpsxConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/hpsxConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hpsx
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hpsxConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hpsxConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/hpsx
)
