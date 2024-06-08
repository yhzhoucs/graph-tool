Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.6.0
)
FetchContent_MakeAvailable(Catch2)

set(TESTS test_static test_stream)
foreach (t IN LISTS TESTS)
    add_executable(${t} ${t}.cpp)
    set_property(TARGET ${t} PROPERTY CXX_STANDARD 20)
    target_compile_definitions(${t} PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")
    target_link_libraries(${t} PRIVATE gtool Catch2::Catch2WithMain)
    add_test(NAME "test${t}" COMMAND ${t})
endforeach ()