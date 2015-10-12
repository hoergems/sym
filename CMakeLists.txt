cmake_minimum_required (VERSION 2.6)

project (integrator)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS "-g -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-g -std=c++11")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${integrator_SOURCE_DIR})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${integrator_SOURCE_DIR})

include_directories(/usr/include/eigen3 
                    /usr/inclue/boost
                    /usr/include/python2.7)

# set the path to the library folder
link_directories(/usr/lib /usr/local/lib)

# link the libraries to the executable

add_library(integrate SHARED integrate.cpp)

target_link_libraries (integrate
                       boost_system
                       boost_python
                       python2.7)
set_property(TARGET integrate PROPERTY CXX_STANDARD 11)

#add_executable(inte main.cpp)
#set_property(TARGET inte PROPERTY CXX_STANDARD 11)
#target_link_libraries (inte
#                       integrate)



                       
#file(COPY "ompl_test" DESTINATION "..")
