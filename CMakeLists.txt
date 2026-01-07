cmake_minimum_required(VERSION 3.20)
project(montecarlo_1 CXX)

# --- Set C++ standard ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- MuParserX build options ---
# Disable building examples and samples to avoid compilation errors
set(MUPARSERX_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(MUPARSERX_BUILD_SAMPLES OFF CACHE BOOL "" FORCE)

# --- Add MuParserX library ---
# It is located in src/external/muparserx
add_subdirectory(src/external/muparserx EXCLUDE_FROM_ALL)
include_directories(src/external/muparserx/parser)

add_executable(montecarlo_1
        src/apps/main.cpp
        src/apps/benchmarks.cpp
)

add_library(montecarlo_optimizers
        src/montecarlo/optimizers/PSO.cpp
        src/montecarlo/optimizers/GA.cpp
)

# Ora si può fare include <montecarlo/...>
target_include_directories(montecarlo_1 PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# --- Find Boost libraries ---
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system)

# --- Link libraries ---
target_link_libraries(montecarlo_1 PRIVATE
        montecarlo_optimizers
        Boost::iostreams
        Boost::filesystem
        Boost::system
)