# Project <proj> main cmake file

cmake_minimum_required(VERSION 3.1)

project(ncov2019 VERSION 0.0.1)

# Option to install targets to env location
if(DEFINED ENV{NCOV2019_HOME})
     message("-- Will install TEMPLATE binaries and libraries to location: " $ENV{NCOV2019_HOME})
	 set(CMAKE_INSTALL_PREFIX $ENV{NCOV2019_HOME})
else(DEFINED ENV{NCOV2019_HOME})
     message("-- Will install TEMPLATE binaries and libraries to location: " ${CMAKE_CURRENT_BINARY_DIR})
	 set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif(DEFINED ENV{NCOV2019_HOME})


# Disable some windows specific warnings for using more modern c functions
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

# Options
option(WITH_testing "Allows tests to be run" ON)
option(WITH_examples "Builds the examples" ON)

if(WITH_testing MATCHES ON)
   enable_testing()
endif(WITH_testing MATCHES ON)

if(WITH_examples MATCHES ON)
   add_subdirectory(examples)
endif(WITH_examples MATCHES ON)

add_executable(covid19 covid19.cpp)

install(TARGETS covid19 DESTINATION bin)

# Install codiv19 gencome
install(FILES data/COVID-19-genome.txt DESTINATION dat)