Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.6.0
)
FetchContent_MakeAvailable(Catch2)

add_executable(test_gtool test.cpp)
target_compile_definitions(test_gtool PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")
target_link_libraries(test_gtool PRIVATE gtool Catch2::Catch2WithMain)
add_test(NAME testatest COMMAND test_gtool)