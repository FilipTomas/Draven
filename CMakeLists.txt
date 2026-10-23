cmake_minimum_required(VERSION 3.11)

project(raven VERSION 1.5.0
              LANGUAGES CXX C
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
    GIT_REPOSITORY https://github.com/FilipTomas/biosoup
    GIT_TAG master)
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
    GIT_TAG internal-count
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

FetchContent_Declare(
  wfa2
  GIT_REPOSITORY https://github.com/smarco/WFA2-lib
  GIT_TAG v2.3.5)
FetchContent_GetProperties(wfa2)
if (NOT wfa2_POPULATED)
  FetchContent_Populate(wfa2)
  file(GLOB WFA2_WAVEFRONT_SOURCES  ${wfa2_SOURCE_DIR}/wavefront/*.c)
  file(GLOB WFA2_UTILS_SOURCES      ${wfa2_SOURCE_DIR}/utils/*.c)
  file(GLOB WFA2_SYSTEM_SOURCES     ${wfa2_SOURCE_DIR}/system/*.c)
  file(GLOB WFA2_ALIGNMENT_SOURCES  ${wfa2_SOURCE_DIR}/alignment/*.c)
  file(GLOB WFA2_BINDINGS_SOURCES   ${wfa2_SOURCE_DIR}/bindings/cpp/*.cpp)
  add_library(wfa2 STATIC
    ${WFA2_WAVEFRONT_SOURCES}
    ${WFA2_UTILS_SOURCES}
    ${WFA2_SYSTEM_SOURCES}
    ${WFA2_ALIGNMENT_SOURCES}
    ${WFA2_BINDINGS_SOURCES}
  )
  target_include_directories(wfa2 PUBLIC ${wfa2_SOURCE_DIR})
  set_target_properties(wfa2 PROPERTIES LINKER_LANGUAGE CXX)
endif ()

# ---------------------------------------------------------
# LibTorch (TorchScript C++ inference) — fetched via FetchContent
# ---------------------------------------------------------
option(RAVEN_FETCH_LIBTORCH "Download LibTorch automatically (CPU, shared, C++11 ABI)" ON)

if (RAVEN_FETCH_LIBTORCH)
  FetchContent_Declare(
    libtorch
    URL https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.5.1%2Bcpu.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )

  FetchContent_GetProperties(libtorch)
  if (NOT libtorch_POPULATED)
    message(STATUS "Fetching LibTorch (CPU, shared-with-deps, cxx11-abi, latest)...")
    FetchContent_Populate(libtorch)
    message(STATUS "LibTorch extracted root: ${libtorch_SOURCE_DIR}")
  endif ()

  file(GLOB_RECURSE _TORCH_CONFIG_CANDIDATES
    RELATIVE "${libtorch_SOURCE_DIR}"
    "${libtorch_SOURCE_DIR}/**/TorchConfig.cmake"
  )

  list(LENGTH _TORCH_CONFIG_CANDIDATES _TORCH_CONFIG_COUNT)
  if (_TORCH_CONFIG_COUNT EQUAL 0)
    message(FATAL_ERROR
      "LibTorch fetched, but TorchConfig.cmake was not found under:\n"
      "  ${libtorch_SOURCE_DIR}\n"
      "Sanity check: the zip should contain share/cmake/Torch/TorchConfig.cmake somewhere.")
  endif ()

  list(GET _TORCH_CONFIG_CANDIDATES 0 _TORCH_CONFIG_REL)
  set(_TORCH_CONFIG_ABS "${libtorch_SOURCE_DIR}/${_TORCH_CONFIG_REL}")

  get_filename_component(_TORCH_DIR "${_TORCH_CONFIG_ABS}" DIRECTORY)

  set(Torch_DIR "${_TORCH_DIR}" CACHE PATH "Torch CMake package directory" FORCE)
  message(STATUS "TorchConfig.cmake found at: ${_TORCH_CONFIG_ABS}")
  message(STATUS "Torch_DIR set to: ${Torch_DIR}")
endif ()

find_package(Torch REQUIRED CONFIG)
message(STATUS "Found Torch version: ${Torch_VERSION}")

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
  src/graph_constructor/candidate_spill.cpp
  src/graph_constructor/cnn_preprocessing.cpp
  src/graph_constructor/overlap_selection.cpp
)

target_link_libraries(raven PUBLIC
  cereal::cereal
  biosoup::biosoup
  bioparser::bioparser
  ram::ram
  edlib::edlib
  wfa2
  catboost::model
  ${TORCH_LIBRARIES}
)

target_compile_options(raven PUBLIC ${TORCH_CXX_FLAGS})

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