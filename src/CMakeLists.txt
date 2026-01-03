# ライブラリの定義（テストでリンクするため）
add_library(mphhc_lib core.cpp)
target_include_directories(mphhc_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

# 実行ファイルの定義
add_executable(mphhc main.cpp)
target_link_libraries(mphhc PRIVATE mphhc_lib)

