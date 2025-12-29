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

# --- Project sources ---
# Collect all .cpp files in src/ and subfolders
file(GLOB_RECURSE SOURCES src/*.cpp)

# === FIX START ===
# Remove anything inside src/external from the list.
# The library is already built by add_subdirectory; we don't want to compile
# its sources (or its documentation files) into our executable manually.
list(FILTER SOURCES EXCLUDE REGEX ".*/src/external/.*")
# === FIX END ===

# Optional: print sources for debugging
message(STATUS "Project sources: ${SOURCES}")

# --- Main executable ---
add_executable(montecarlo_1 ${SOURCES})

# --- Find Boost libraries ---
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system)

# --- Link libraries ---
target_link_libraries(montecarlo_1
        PRIVATE
        Boost::iostreams
        Boost::filesystem
        Boost::system
        muparserx
)
