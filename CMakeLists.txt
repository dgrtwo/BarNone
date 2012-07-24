cmake_minimum_required(VERSION 2.6)

# files to compile
set(APPLICATION_EXEC_SRC 
        src/BarNone.cpp
        src/SeqParser.cpp
        src/BarcodeCache.cpp
        src/BarcodeCounter.cpp
        src/BkTree.cpp
)

# where to look for header files
include_directories (
    .
    ./src/flamingo-4.1/src/
    ./src/flamingo-4.1/src/filtertree/
    ./src/flamingo-4.1/src/filtertree/src/
    ./src/flamingo-4.1/src/stringmap/src/
    ${CMAKE_SOURCE_DIR}/../
    include
    lib    
)

# where to look for dependent libraries
link_directories(
    ${CMAKE_SOURCE_DIR}/src/flamingo-4.1/src/common/build/
    ${CMAKE_SOURCE_DIR}/src/flamingo-4.1/src/util/build/
    ${CMAKE_SOURCE_DIR}/src/flamingo-4.1/src/listmerger/build/
    ${CMAKE_SOURCE_DIR}/src/flamingo-4.1/src/filtertree/build/
    ${CMAKE_SOURCE_DIR}/src/flamingo-4.1/src/stringmap/build/
) 

# have cmake also build the filtertree module (if not built already)
# the filtertree module will build util, common and listmerger
add_subdirectory(./src/flamingo-4.1/src/filtertree/ ./flamingo-4.1/src/filtertree/)
#add_subdirectory(./src/flamingo-4.1/src/stringmap/  ./flamingo-4.1/src/stringmap/)

# GCC command line args
add_definitions(-Wall -O3 -DDEBUG_TIMER_FANCY -DDEBUG_STAT -DED_MATRIX_DIM=2000)

# create executable
add_executable(BarNone ${APPLICATION_EXEC_SRC})
add_dependencies(BarNone wrappers-lib filtertree-lib common-lib util-lib listmerger-lib)
target_link_libraries(BarNone wrappers-lib filtertree-lib common-lib util-lib listmerger-lib)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/")
set(CMAKE_EXE_LINKER_FLAGS -ggdb -g)
