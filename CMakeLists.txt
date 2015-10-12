cmake_minimum_required (VERSION 2.6)

project (ompl_test)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS "-g -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-g -std=c++11")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ompl_test_SOURCE_DIR})

include_directories(/usr/include/eigen3 
                    /usr/inclue/boost
                    /usr/include/python2.7)

# set the path to the library folder
link_directories(/usr/lib /usr/local/lib)

# link the libraries to the executable


add_library(libintegrate SHARED integrate.cpp)
set_property(TARGET libintegrate PROPERTY CXX_STANDARD 11)
target_link_libraries (libintegrate
                       boost_system
                       boost_python)
                       
#file(COPY "ompl_test" DESTINATION "..")
