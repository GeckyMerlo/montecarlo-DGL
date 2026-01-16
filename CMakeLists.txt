cmake_minimum_required(VERSION 3.20)
project(montecarlo_1 CXX)

# --- Suppress all compiler warnings ---
add_compile_options(-w)

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
        src/apps/benchmarks/integration_benchmarks.cpp
        src/apps/benchmarks/pso_benchmarks.cpp
        src/apps/benchmarks/ga_benchmarks.cpp
)

add_executable(drone_optimization 
        src/apps/drone_optimization.cpp
)

add_executable(drone_optimization_ellipsoid_strength 
        src/apps/drone_optimization_ellipsoid_strength.cpp
)

add_library(montecarlo_optimizers
        src/montecarlo/optimizers/PSO.cpp
        src/montecarlo/optimizers/GA.cpp
)

# Ora si può fare include <montecarlo/...>
target_include_directories(montecarlo_1 PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(drone_optimization PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(drone_optimization_ellipsoid_strength PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# --- Find Boost libraries ---
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system)

find_package(OpenMP REQUIRED)

# --- Link libraries ---
target_link_libraries(montecarlo_1 PRIVATE
        OpenMP::OpenMP_CXX
        montecarlo_optimizers
        Boost::iostreams
        Boost::filesystem
        Boost::system
        muparserx
)

target_link_libraries(drone_optimization PRIVATE
        OpenMP::OpenMP_CXX
        montecarlo_optimizers
        Boost::iostreams
        Boost::filesystem
        Boost::system
        muparserx
)

target_link_libraries(drone_optimization_ellipsoid_strength PRIVATE
        OpenMP::OpenMP_CXX
        montecarlo_optimizers
        Boost::iostreams
        Boost::filesystem
        Boost::system
        muparserx
)