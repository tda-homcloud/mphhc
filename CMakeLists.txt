cmake_minimum_required(VERSION 3.14)
project(mphhc VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- GoogleTest の自動ダウンロード設定 ---
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/heads/main.zip
)
# For Windows: GTest が MSVC のランタイムライブラリ設定と矛盾しないようにする
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# --- 各ディレクトリの追加 ---
add_subdirectory(src)

# テストの有効化
enable_testing()
add_subdirectory(tests)
