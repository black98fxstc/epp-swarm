cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE PATH "toolchain file")
endif()
set(VCPKG "vcpkg_installed/x86-linux")
set(CMAKE_PREFIX_PATH "/usr/local;${VCPKG}")

project(EPPSwarm LANGUAGES CXX VERSION 0.1)

set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
set(BUILD_SHARED_LIBS ON)

find_package(Threads REQUIRED)
find_package(FFTW3 CONFIG REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_library(MY_LIBRARIES NAMES fftw3 fftw3f PATHS "${VCPKG}/lib")

if(NOT MSVC)
    include_directories(/usr/local/include)
    link_directories(/usr/local/lib)
endif(NOT MSVC)

find_library(FFTW3_LIBRARY fftw3)
find_library(FFTW3F_LIBRARY fftw3f)

set(EPP_H constants.h client.h worker.h boundary.h modal.h pursuit.h)

add_library(epp STATIC client.cxx polygon.cxx json.cxx stream.cxx "${EPP_H}")

add_executable(EPPtest testing.cxx constants.h client.h MATLAB.h)
target_link_libraries(EPPtest PRIVATE epp "${FFTW3_LIBRARY}" "${FFTW3F_LIBRARY}" "${CMAKE_THREAD_LIBS_INIT}")
target_link_libraries(EPPtest PRIVATE nlohmann_json::nlohmann_json)

add_executable(EPPcli cli.cxx constants.h client.h MATLAB.h)
target_link_libraries(EPPcli PRIVATE epp "${FFTW3_LIBRARY}" "${FFTW3F_LIBRARY}" "${CMAKE_THREAD_LIBS_INIT}")
target_link_libraries(EPPcli PRIVATE nlohmann_json::nlohmann_json)

if (WAM)
    find_package(nlohmann_json 3.2.0 REQUIRED)
    # add_executable(EPPworker worker.cxx "${EPP_H}")
    # target_link_libraries(EPPworker PRIVATE epp "${FFTW_LIBRARIES}" "${CMAKE_THREAD_LIBS_INIT}")
    # target_link_libraries(EPPworker PRIVATE nlohmann_json::nlohmann_json)
endif(WAM)

# find_package(AWSSDK REQUIRED COMPONENTS s3)
# add_library(epp STATIC include/client.h client.cxx sample.cxx subset.cxx)
# target_link_libraries(EPPWorker epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPDummy dummy.cxx include/client.h include/work.h)
# target_link_libraries(EPPDummy epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPtest testing.cxx pursuit.cxx constants.h pursuit.h pursuer.h boundary.h modal.h)
# target_link_libraries(EPPtest fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )