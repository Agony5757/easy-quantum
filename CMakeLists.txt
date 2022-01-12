cmake_minimum_required(VERSION 3.1)

project(EasyQuantum)

include_directories(${PROJECT_SOURCE_DIR}/easy_quantum/include
                    ${PROJECT_SOURCE_DIR}/ThirdParty/fmt-8.0.1/include
                    )

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_CXX_STANDARD 17)

if(UNIX)
    add_compile_options(-fPIC -fpermissive)
endif()

# find_package(PythonInterp 3 REQUIRED)
# find_package(PythonLibs 3 REQUIRED)
# include_directories(${PYTHON_INCLUDE_DIRS})
# include_directories(${PROJECT_SOURCE_DIR}/ThirdParty/pybind11)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

add_subdirectory(easy_quantum)
add_subdirectory(test)