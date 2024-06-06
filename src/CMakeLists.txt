add_executable(demo1 demo1.cpp)
add_executable(demo2 demo2.cpp)

target_compile_definitions(demo1 PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")
target_compile_definitions(demo2 PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")

target_link_libraries(demo1 PRIVATE gtool)
target_link_libraries(demo2 PRIVATE gtool)