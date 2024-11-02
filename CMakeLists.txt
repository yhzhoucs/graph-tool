cmake_minimum_required(VERSION 3.25.0)
project(gtool
    VERSION 0.1.2
    DESCRIPTION "A simple graph processing tool"
    LANGUAGES CXX)

option(ENABLE_BUILD_DEMO "enable build demos" OFF)
option(ENABLE_BUILD_TEST "enable test" OFF)
option(BUILD_WITH_RESTORE "build gtool with restore" OFF)

add_library(gtool INTERFACE)
add_library(gtool::gtool ALIAS gtool)
target_compile_features(gtool INTERFACE cxx_std_11)
target_include_directories(gtool INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
if (BUILD_WITH_RESTORE)
    target_compile_definitions(gtool INTERFACE BUILD_WITH_RESTORE)
endif ()

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP Found")
    target_link_libraries(gtool INTERFACE OpenMP::OpenMP_CXX)
else ()
    message(STATUS "OpenMP Not Found")
endif ()

if (ENABLE_BUILD_DEMO)
    add_subdirectory(src)
endif ()

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND ENABLE_BUILD_TEST)
    include(CTest)
    add_subdirectory(tests)
endif ()

# install and pack
install(TARGETS gtool
        EXPORT gtool
        DESTINATION lib/cmake/gtool)

install(EXPORT gtool
        FILE gtool.cmake
        NAMESPACE gtool::
        DESTINATION lib/cmake/gtool)

install(DIRECTORY include/ DESTINATION include
        FILES_MATCHING PATTERN "*.h")

include(CMakePackageConfigHelpers)
configure_package_config_file(
        cmake/gtoolConfig.cmake.in ${PROJECT_BINARY_DIR}/gtoolConfig.cmake
        INSTALL_DESTINATION lib/cmake/gtool)

write_basic_package_version_file(
        gtoolConfigVersion.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY SameMajorVersion)

install(FILES "${PROJECT_BINARY_DIR}/gtoolConfig.cmake"
        "${PROJECT_BINARY_DIR}/gtoolConfigVersion.cmake"
        DESTINATION lib/cmake/gtool)

set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
include(CPack)