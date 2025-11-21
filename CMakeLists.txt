cmake_minimum_required(VERSION 3.11)

project(raven VERSION 1.5.0
              LANGUAGES CXX
              DESCRIPTION "Raven is a de novo genome assembler for long uncorrected reads.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)
include(ExternalProject)

# ---------------------------------------------------------
# Options
# ---------------------------------------------------------
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(raven_main_project ON)
endif ()
option(raven_build_tests "Build unit tests" ${raven_main_project})

# ---------------------------------------------------------
# Dependencies
# ---------------------------------------------------------
find_package(biosoup 0.10.0 QUIET)
if (NOT biosoup_FOUND)
  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)
  FetchContent_GetProperties(biosoup)
  if (NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(${biosoup_SOURCE_DIR} ${biosoup_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(bioparser 3.0.13 QUIET)
if (NOT bioparser_FOUND)
  FetchContent_Declare(
    bioparser
    GIT_REPOSITORY https://github.com/rvaser/bioparser
    GIT_TAG 3.0.13)
  FetchContent_GetProperties(bioparser)
  if (NOT bioparser_POPULATED)
    FetchContent_Populate(bioparser)
    add_subdirectory(${bioparser_SOURCE_DIR} ${bioparser_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(cereal 1.3.0 QUIET)
if (NOT cereal_FOUND)
  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal
    GIT_TAG v1.3.0)
  FetchContent_GetProperties(cereal)
  if (NOT cereal_POPULATED)
    FetchContent_Populate(cereal)
    add_subdirectory(${cereal_SOURCE_DIR} ${cereal_BINARY_DIR} EXCLUDE_FROM_ALL)
    add_library(cereal::cereal ALIAS cereal)
  endif ()
endif ()

find_package(ram 1.0.0 QUIET)
if (NOT ram_FOUND)
  FetchContent_Declare(
    ram
    GIT_REPOSITORY https://github.com/FilipTomas/ram_count
    GIT_TAG master
    GIT_SUBMODULES "")
  FetchContent_GetProperties(ram)
  if (NOT ram_POPULATED)
    FetchContent_Populate(ram)
    add_subdirectory(${ram_SOURCE_DIR} ${ram_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif ()
endif ()

find_package(edlib 1.2.7 QUIET)
if (NOT edlib_FOUND)
  FetchContent_Declare(
    edlib
    GIT_REPOSITORY https://github.com/martinsos/edlib
    GIT_TAG v1.2.7)
  FetchContent_GetProperties(edlib)
  if (NOT edlib_POPULATED)
    FetchContent_Populate(edlib)
    add_subdirectory(${edlib_SOURCE_DIR} ${edlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif ()
endif ()

# ---------------------------------------------------------
# CatBoost (header-only, ALWAYS ON)
# ---------------------------------------------------------
set(RAVEN_CATBOOST_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/third_party/catboost_model/catboost_model.hpp")

if (NOT EXISTS "${RAVEN_CATBOOST_HEADER}")
  message(FATAL_ERROR "CatBoost header '${RAVEN_CATBOOST_HEADER}' not found.\n"
    "Export from Python first:\n"
    "  model.save_model('third_party/catboost_model/catboost_model.hpp', format='cpp')")
endif()

add_library(catboost_header_only INTERFACE)
target_include_directories(catboost_header_only INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/catboost_model)
target_compile_definitions(catboost_header_only INTERFACE RAVEN_CATBOOST_HEADER=1)
add_library(catboost::model ALIAS catboost_header_only)

# # ---------------------------------------------------------
# # KMC integration (fast k-mer counter, k > 31) — ALWAYS ON
# # ---------------------------------------------------------
# FetchContent_Declare(
#   KMC
#   GIT_REPOSITORY https://github.com/refresh-bio/KMC.git
#   GIT_TAG v3.2.2
#   GIT_SHALLOW TRUE
# )
# FetchContent_GetProperties(KMC)
# if (NOT KMC_POPULATED)
#   FetchContent_Populate(KMC)
# endif()

# # Note: older CMake exposes lowercase variables for FetchContent:
# #  kmc_SOURCE_DIR / kmc_BINARY_DIR
# if (NOT EXISTS "${kmc_SOURCE_DIR}/kmc_api")
#   message(FATAL_ERROR
#     "KMC API sources not found in '${kmc_SOURCE_DIR}/kmc_api'. "
#     "Resolved kmc_SOURCE_DIR='${kmc_SOURCE_DIR}'. "
#     "Check that the KMC repository was fetched correctly.")
# endif()

# # Build KMC (kmc, kmc_tools) via its Makefile
# find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)
# add_custom_target(kmc_build
#   COMMAND ${MAKE_EXECUTABLE} -C ${kmc_SOURCE_DIR}
#   WORKING_DIRECTORY ${kmc_SOURCE_DIR}
#   COMMENT "Building KMC (kmc, kmc_tools) via Makefile"
#   VERBATIM)

# # Build KMC C++ API as a static library
# file(GLOB KMC_API_SRC "${kmc_SOURCE_DIR}/kmc_api/*.cpp")
# add_library(kmc_api STATIC ${KMC_API_SRC})
# add_dependencies(kmc_api kmc_build)
# target_include_directories(kmc_api PUBLIC "${kmc_SOURCE_DIR}/kmc_api")

# # Link zlib + pthread
# find_package(ZLIB REQUIRED)
# if (UNIX)
#   target_link_libraries(kmc_api PUBLIC ZLIB::ZLIB pthread)
# else()
#   target_link_libraries(kmc_api PUBLIC ZLIB::ZLIB)
# endif()

# # Expose path to KMC binaries for runtime helpers
# target_compile_definitions(kmc_api PUBLIC KMC_BIN_DIR="${kmc_SOURCE_DIR}/bin")

# ---------------------------------------------------------
# Your targets
# ---------------------------------------------------------
add_library(raven
  src/graph.cpp
  src/pile.cpp
  src/node.cpp
  src/edge.cpp
  src/parser.cpp
  src/graph_constructor/graph_constructor.cpp
  src/graph_assembler/graph_assembler.cpp
  src/option_manager.cpp
  src/graph_constructor/overlap_parser.cpp
  src/graph_constructor/ram_overlap.cpp
  #src/extended_overlap.cpp
  #src/read_preprocessing.cpp
)

target_link_libraries(raven PUBLIC
  cereal::cereal
  biosoup::biosoup
  bioparser::bioparser
  ram::ram
  edlib::edlib
  #kmc_api
  catboost::model
)

target_include_directories(raven PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

add_executable(raven_exe
  src/main.cpp
)

target_link_libraries(raven_exe PUBLIC raven)

target_compile_definitions(raven PUBLIC RAVEN_VERSION="${PROJECT_VERSION}")
set_property(TARGET raven_exe PROPERTY OUTPUT_NAME raven)

install(TARGETS raven_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

# ---------------------------------------------------------
# Tests
# ---------------------------------------------------------
if (raven_build_tests)
  find_package(GTest 1.14.0 QUIET)
  if (NOT GTest_FOUND)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG v1.14.0)
    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()

  enable_testing()
  add_executable(raven_test test/raven_test.cpp)
  target_link_libraries(raven_test PUBLIC raven bioparser::bioparser GTest::Main)
  target_compile_definitions(raven_test PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/")
  gtest_add_tests(TARGET raven_test SOURCES test/raven_test.cpp)
endif ()
