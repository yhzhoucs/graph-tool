cmake_minimum_required(VERSION 3.12)
project(gtool
    VERSION 0.0.1
    DESCRIPTION "A simple graph processing tool"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(gtool INTERFACE)
target_include_directories(gtool INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

add_executable(demo src/main.cpp)
target_link_libraries(demo PRIVATE gtool)
target_compile_definitions(demo PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP Found")
    target_link_libraries(demo PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(gtool INTERFACE OpenMP::OpenMP_CXX)
endif()

install(TARGETS gtool
        EXPORT gtool
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(DIRECTORY include/ DESTINATION include/gtool
        FILES_MATCHING PATTERN "*.h")

install(
        EXPORT gtool
        FILE gtool.cmake
        NAMESPACE gtool::
        DESTINATION lib/cmake/gtool
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
        gtoolConfig.cmake.in ${PROJECT_BINARY_DIR}/gtoolConfig.cmake
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