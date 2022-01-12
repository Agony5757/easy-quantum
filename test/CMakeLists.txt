cmake_minimum_required(VERSION 3.1)

file(GLOB Test_SRC_CPP ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
message("Find Test : ${Test_SRC_CPP}")
add_executable(testmain testmain.cpp ${Test_SRC_CPP})
target_link_libraries(testmain EasyQuantumStatic)