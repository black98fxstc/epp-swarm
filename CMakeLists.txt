cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_PREFIX_PATH "/usr/local")
set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
set(BUILD_SHARED_LIBS ON)

project(EPPSwarm LANGUAGES CXX VERSION 0.1)
include_directories(include)
find_package(AWSSDK REQUIRED COMPONENTS s3)
find_package(nlohmann_json 3.2.0 REQUIRED)

# add the executable
add_executable(EPPWorker worker.cxx include/client.h client.cxx sample.cxx)
# add_executable(Sample sample.cxx include/client.h)
target_link_libraries(EPPWorker curl ${AWSSDK_LINK_LIBRARIES})
target_include_directories(EPPWorker PUBLIC "${PROJECT_BINARY_DIR}")