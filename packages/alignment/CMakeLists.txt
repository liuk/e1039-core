# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(alignment CXX C Fortran)

# source code
file(GLOB          sources ${PROJECT_SOURCE_DIR}/*.cxx)
file(GLOB     dist_headers ${PROJECT_SOURCE_DIR}/*.h)
file(GLOB non_dist_headers ${PROJECT_SOURCE_DIR}/AlignmentLinkDef.h)
list(APPEND sources ${PROJECT_SOURCE_DIR}/millepede.f)
list(REMOVE_ITEM dist_headers ${non_dist_headers})

# ROOT dict generation - do I need to remove SQMillepede.* from this list?
add_custom_command (
  OUTPUT AlignmentDict.cxx
  COMMAND rootcint
  ARGS -f AlignmentDict.cxx -noIncludePaths -inlineInputHeader -c -p -I${PROJECT_SOURCE_DIR}/ -I${ROOT_PREFIX}/include/ ${dist_headers} ${non_dist_headers}
) 

# ROOT
find_program(ROOTCONF "root-config")
if(ROOTCONF)
  message("-- Detecting ROOT:    found at ${ROOTCONF}")
else()
  message(FATAL_ERROR "-- Detecting ROOT:    not found")
endif()
execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOT_PREFIX  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs   OUTPUT_VARIABLE ROOT_LINK    OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS}")

include_directories("${PROJECT_SOURCE_DIR}/")

add_library(alignment SHARED ${sources} AlignmentDict.cxx)
target_link_libraries(alignment -lfun4all -lphool ${ROOT_LINK} -lgeom_svc -linterface_main -lktracker_interface -lgfortran)

install(TARGETS alignment     DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${dist_headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${CMAKE_PROJECT_NAME}/)

add_custom_target(install_pcm ALL COMMAND mkdir -p ${CMAKE_INSTALL_PREFIX}/lib COMMAND cp -up *_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)
add_dependencies(install_pcm alignment)

# enable/disable test code build
option(BUILDTEST "Build the test code for the package" OFF)
if(BUILDTEST)
  file(GLOB tests ${PROJECT_SOURCE_DIR}/test/test*.cxx)
  foreach(test ${tests})
    string(REPLACE "${PROJECT_SOURCE_DIR}/test/" "" test ${test})
    string(REPLACE ".cxx" "" test ${test})
    message(STATUS "Adding tests to the compile list: ${test}")

    add_executable(${test} ${PROJECT_SOURCE_DIR}/test/${test}.cxx)
    target_link_libraries(${test} ${PROJECT_BINARY_DIR}/libalignment.so -lktracker -lktracker_interface -lgeom_svc -ldb_svc -lphool "${ROOT_LINK}")
    add_dependencies(${test} alignment)
    install(TARGETS ${test} DESTINATION ${PROJECT_SOURCE_DIR}/test)
  endforeach(test)
endif()
