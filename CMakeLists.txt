#Made in VS code with ctrl+shift+p and cmake:quick
cmake_minimum_required(VERSION 3.0.0)

#> Variables 
# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -O0 -g3 -pthread \ #TODO: check the flags.. do we need pthread ?
-lsdsl -ldivsufsort -ldivsufsort64")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra -O3 -g0 -pthread \
#-lsdsl -ldivsufsort -ldivsufsort64")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O9 -std=c++11 -DNODEBUG \
-lsdsl -ldivsufsort -ldivsufsort64")
#< Variables 

project(Ring VERSION 0.1.0)

include(CTest)
enable_testing()

find_package(Boost 1.71.0 REQUIRED)
#If the above library is found, then several result variables are set, such as Boost_INCLUDE_DIRS. 
#I.e. :https://cmake.org/cmake/help/latest/module/FindBoost.html
include_directories(${Boost_INCLUDE_DIRS})

#> SDSL Specifics
include_directories(/home/fabrizio/include) #compiler
link_directories(/home/fabrizio/lib) #linker
#< SDSL Specific

#Retrieve recursively all cpp files within 'src' https://cmake.org/cmake/help/latest/command/file.html#glob-recurse
file(GLOB_RECURSE SRCS src/*.cpp)
#And then set them as 'sources'
set(sources ${SRCS})
# Include ring's 'sources' as a static library (used later by the linker)
add_library(ringlib STATIC ${sources})

#Setting up the source code
set(BUILD_TARGETS
    build-index
    query-index
    query-index-orden-promisorio
)

foreach(target ${BUILD_TARGETS})
    add_executable(${target} src/${target}.cpp)
    target_link_libraries(${target} ringlib boost_program_options boost_system sdsl divsufsort divsufsort64)
endforeach(target)

#Setting up tests
set(TEST_TARGETS
    test_colored_range_counting_wm
    test_colored_range_counting_wm2
    test_initializable_bitmap
)

add_test(NAME test_vm COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/build/test_wm)
foreach(target ${TEST_TARGETS})
    add_executable(${target} src/tests/${target}.cpp)
    target_link_libraries(${target} ringlib boost_program_options boost_system sdsl divsufsort divsufsort64)
endforeach(target)

#Others: todo investigate what cpack does.
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
