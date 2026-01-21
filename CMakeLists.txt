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

# Test for RNG reproducibility
add_executable(test_rng_reproducibility 
        src/apps/test_rng_reproducibility.cpp
)

add_library(montecarlo_optimizers
        src/montecarlo/rng/rng_global.cpp
        src/montecarlo/rng/rng_factory.cpp
        src/montecarlo/optimizers/PSO.cpp
        src/montecarlo/optimizers/GA.cpp
)

add_executable(wind_farm_simulator 
        src/apps/wind_farm_simulator.cpp
)

# montecarlo_optimizers needs OpenMP for omp_get_thread_num() in rng_factory.cpp
target_link_libraries(montecarlo_optimizers PRIVATE OpenMP::OpenMP_CXX)

# Include path for montecarlo_optimizers
target_include_directories(montecarlo_optimizers PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Ora si può fare include <montecarlo/...>
target_include_directories(montecarlo_1 PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(drone_optimization PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(wind_farm_simulator PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(test_rng_reproducibility PRIVATE
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

# Test for RNG reproducibility - only needs montecarlo_optimizers
target_link_libraries(test_rng_reproducibility PRIVATE
        OpenMP::OpenMP_CXX
        montecarlo_optimizers
)

target_link_libraries(wind_farm_simulator PRIVATE
        OpenMP::OpenMP_CXX
        montecarlo_optimizers
        Boost::iostreams    
        Boost::filesystem
        Boost::system
        muparserx           
)