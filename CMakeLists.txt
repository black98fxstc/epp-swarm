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
include_directories(/usr/local/include)


# add the executable
add_library(epp STATIC include/client.h client.cxx sample.cxx subset.cxx include/work.h include/modal.h modal.cxx include/boundary.h)
# add_executable(EPPWorker worker.cxx include/client.h include/work.h include/modal.h include/boundary.h)
# target_link_libraries(EPPWorker epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPDummy dummy.cxx include/client.h include/work.h)
# target_link_libraries(EPPDummy epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
add_executable(EPPtest testing.cxx include/boundary.h include/modal.h)
target_link_libraries(EPPtest epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )