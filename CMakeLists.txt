cmake_minimum_required (VERSION 3.24)
# cmake_policy(VERSION 3.0)
project (1d_nodal_dg_project)
enable_language(CXX)

# Version
set(VERSION_MAJOR 2)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(PROJECT_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
# Undefine
set(VERSION_MAJOR)
set(VERSION_MINOR)
set(VERSION_PATCH)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# Select the build type:
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES ON)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
# reuse already downloaded content without updating:
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
# # for offline builds if dependencies are present
# set(FETCHCONTENT_FULLY_DISCONNECTED ON)

# -------------------------

if(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
  MESSAGE(STATUS "Using custom installation prefix: ${CMAKE_BINARY_DIR}/install")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# -------------------------
# ----- find packages -----
# -------------------------

set(DEPENDENCIES)
set(TEST_DEPENDENCIES)

# -------------------------

set(BOOST_INCLUDE_LIBRARIES program_options thread system filesystem)
set(BOOST_ENABLE_CMAKE ON)

FetchContent_Declare(
  Boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-1.84.0
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  FIND_PACKAGE_ARGS COMPONENTS program_options thread system filesystem
)
FetchContent_MakeAvailable(Boost)
LIST(APPEND DEPENDENCIES Boost::program_options Boost::thread Boost::system Boost::filesystem)

# -------------------------

FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(fmt)
find_package(fmt)
if(NOT ${fmt_FOUND})
 message(ERROR "fmt is required but was not found!")
endif()
LIST(APPEND DEPENDENCIES fmt::fmt)

# -------------------------

set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_BUILD_DOC OFF)

# target version is 3.5.0, but it was not available at the time of writing
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG "master"
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  FIND_PACKAGE_ARGS NO_MODULE
)
FetchContent_MakeAvailable(eigen)
if(NOT ${Eigen3_FOUND})
 message(ERROR "eigen3 is required but was not found!")
endif()
LIST(APPEND DEPENDENCIES Eigen3::Eigen)

# -------------------------

FetchContent_Declare(Catch2
  GIT_REPOSITORY "https://github.com/catchorg/Catch2.git"
  GIT_TAG "v3.5.3"
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  FIND_PACKAGE_ARGS 3
)
FetchContent_MakeAvailable(Catch2)
if(NOT ${catch2_FOUND})
 message(ERROR "catch2 is required but was not found!")
endif()
LIST(APPEND TEST_DEPENDENCIES Catch2::Catch2WithMain)

# -------------------------

find_package(OpenMP)
if(${OPENMP_FOUND})
 LIST(APPEND DEPENDENCIES OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP was not found!")
endif()

# -------------------------

# Compiler Flags:
IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # clang
  set(CMAKE_CXX_COMPILER_IS_CLANG TRUE)

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -DDEBUG -O0 -g -fdiagnostics-color=auto -fsanitize=undefined -fsanitize=address")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -DNDEBUG -O3 -flto -march=native -fdiagnostics-color=auto")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O2 -g -DNDEBUG")

ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # gcc
  set(CMAKE_CXX_COMPILER_IS_GCC TRUE)

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -DDEBUG -O0 -g -fdiagnostics-color=auto -fsanitize=undefined -fsanitize=address")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wall -Wfatal-errors -DNDEBUG -O3 -flto -march=native -fdiagnostics-color=auto")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O2 -g -DNDEBUG")

endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  MESSAGE(STATUS "Compiler flags = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif()
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  MESSAGE(STATUS "Compiler flags = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()
if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  MESSAGE(STATUS "Compiler flags = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

# -------------------------

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in ${CMAKE_BINARY_DIR}/include/${PROJECT_NAME}/config.hpp)

# -------------------------

# get global source files:
file(GLOB_RECURSE Program_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/*.cpp")
#file( GLOB APP_SOURCES progs/*.cpp )
foreach(srcfile ${Program_SOURCES})
 get_filename_component(program_name ${srcfile} NAME_WE)
 #string(REPLACE ".cpp" "" program_name ${srcfile} )
 set(program_name "${program_name}.exe")

 MESSAGE(STATUS "Add executable '${program_name}' from file '${srcfile}' ")
 add_executable(${program_name} ${srcfile})
 target_link_libraries(${program_name} PUBLIC ${DEPENDENCIES})

 target_include_directories(${program_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

 # for config.hpp.in
 target_include_directories(${program_name} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include)

 INSTALL(PROGRAMS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${program_name}" DESTINATION bin)
endforeach()

# -------------------------

message(STATUS "INCLUDE_DIRECTORIES:")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "    include_dir='${dir}'")
endforeach()

message(STATUS "SOURCES:")
foreach(src ${Program_SOURCES})
 message(STATUS "    src='${src}'")
endforeach()

MESSAGE(STATUS "LIBRARIES:")
foreach(lib ${DEPENDENCIES})
  message(STATUS "    lib='${lib}'")
endforeach()

# -------------------------

# add a target to generate documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND)

# -------------------------

include(CMakeDependentOption)
cmake_dependent_option(BUILD_TESTING "Enable creation of tests." ON "PROJECT_IS_TOP_LEVEL" OFF)
option(NDG_BUILD_TESTING "Enable creation of NDG tests." ${BUILD_TESTING})
if(NDG_BUILD_TESTING)
  add_subdirectory(tests)
  enable_testing()
endif()
