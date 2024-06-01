cmake_minimum_required(VERSION 3.12)
project(Experiment)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(gtool INTERFACE)
target_include_directories(gtool INTERFACE include)

add_executable(demo src/main.cpp)
target_link_libraries(demo PRIVATE gtool)
target_compile_definitions(demo PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP Found")
    target_link_libraries(demo PUBLIC OpenMP::OpenMP_CXX)
endif()
