cmake_minimum_required(VERSION 3.12.0)
project(lab05 LANGUAGES CXX C)

include(../common/cmake-common.cmake)


file(GLOB SOURCES src/*.cpp)
add_executable(progetto ${SOURCES})

# Trova Boost (moduli necessari a gnuplot-iostream)
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system)

# Linka le librerie Boost al tuo eseguibile
target_link_libraries(progetto
    PRIVATE
        Boost::iostreams
        Boost::filesystem
        Boost::system
)
