cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE PATH "toolchain file")
endif()
set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")

set(BUILD_SHARED_LIBS ON)

project(epp-swarm LANGUAGES CXX VERSION 0.3.0)

find_package(Threads REQUIRED)
find_package(FFTW3 CONFIG REQUIRED)
find_package(FFTW3f CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

set(EPP_H constants.h client.h worker.h boundary.h modal.h pursuit.h)

add_library(epp STATIC client.cxx polygon.cxx json.cxx stream.cxx "${EPP_H}")

add_executable(EPPtest testing.cxx constants.h client.h MATLAB.h)
target_link_libraries(EPPtest PRIVATE epp "${CMAKE_THREAD_LIBS_INIT}")
target_link_libraries(EPPtest PRIVATE FFTW3::fftw3 FFTW3::fftw3f)
target_link_libraries(EPPtest PRIVATE nlohmann_json::nlohmann_json)

add_executable(EPPcli cli.cxx constants.h client.h MATLAB.h)
target_link_libraries(EPPcli PRIVATE epp "${CMAKE_THREAD_LIBS_INIT}")
target_link_libraries(EPPcli PRIVATE FFTW3::fftw3 FFTW3::fftw3f)
target_link_libraries(EPPcli PRIVATE nlohmann_json::nlohmann_json)

# add_executable(EPPworker worker.cxx "${EPP_H}")
# target_link_libraries(EPPworker PRIVATE epp "${FFTW_LIBRARIES}" "${CMAKE_THREAD_LIBS_INIT}")
# target_link_libraries(EPPworker PRIVATE nlohmann_json::nlohmann_json)

# find_package(AWSSDK REQUIRED COMPONENTS s3)
# add_library(epp STATIC include/client.h client.cxx sample.cxx subset.cxx)
# target_link_libraries(EPPWorker epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPDummy dummy.cxx include/client.h include/work.h)
# target_link_libraries(EPPDummy epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPtest testing.cxx pursuit.cxx constants.h pursuit.h pursuer.h boundary.h modal.h)
# target_link_libraries(EPPtest fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )