set(DEMOS demo1 demo2)

foreach (t IN LISTS DEMOS)
    add_executable(${t} ${t}.cpp)
    set_property(TARGET ${t} PROPERTY CXX_STANDARD 20)
    target_compile_definitions(${t} PRIVATE DATASET_PATH="${PROJECT_SOURCE_DIR}/dataset")
    target_link_libraries(${t} PRIVATE gtool)
endforeach ()
