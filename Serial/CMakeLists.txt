cmake_minimum_required(VERSION 3.10.1)

project(project1)

# Include headers
include_directories(include)

add_library(burger src/burger.c include/burger.h)

# Adding all sources
file(GLOB SOURCES "src/*.c")
#set(SOURCES src/test.c src/burger.c)

add_executable(main-c src/main.c)
add_executable(test-c src/test.c)

target_link_libraries(main-c PUBLIC burger)
target_link_libraries(test-c PUBLIC burger)
