cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(WAM ON)
if (WAM)
    SET(VCPKG "C:\\Users\\wmoore\\source\\repos\\vcpkg\\packages")
endif(WAM)

set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
set(BUILD_SHARED_LIBS ON)

project(EPPSwarm LANGUAGES CXX VERSION 0.1)

# include_directories(include)
include_directories(/usr/local/include)

find_package(Threads REQUIRED)
if (MSVC)
    set(CMAKE_PREFIX_PATH "/usr/local;${VCPKG}")
    find_package(FFTW3 CONFIG REQUIRED)
else(MSVC)
    SET(FFTW_LIBRARIES fftw3 fftw3f)
endif(MSVC)

if (APPLE)
    link_directories(/usr/local/lib)
endif(APPLE)
set(EPP_H constants.h client.h boundary.h modal.h pursuit.h)
add_library(epp STATIC MATLAB.cxx polygon.cxx json.cxx stream.cxx "${EPP_H}")
add_executable(EPPtest testing.cxx constants.h client.h)
target_link_libraries(EPPtest PRIVATE epp "${FFTW_LIBRARIES}" "${CMAKE_THREAD_LIBS_INIT}")

if (WAM)
    find_package(nlohmann_json 3.2.0 REQUIRED)
    add_executable(EPPworker worker.cxx "${EPP_H}")
    target_link_libraries(EPPworker PRIVATE epp "${FFTW_LIBRARIES}" "${CMAKE_THREAD_LIBS_INIT}")
    target_link_libraries(EPPworker PRIVATE nlohmann_json::nlohmann_json)
ENDIF(WAM)

# find_package(AWSSDK REQUIRED COMPONENTS s3)
# add_library(epp STATIC include/client.h client.cxx sample.cxx subset.cxx)
# target_link_libraries(EPPWorker epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPDummy dummy.cxx include/client.h include/work.h)
# target_link_libraries(EPPDummy epp curl fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )
# add_executable(EPPtest testing.cxx pursuit.cxx constants.h pursuit.h pursuer.h boundary.h modal.h)
# target_link_libraries(EPPtest fftw3f fftw3 ${AWSSDK_LINK_LIBRARIES} )