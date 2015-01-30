### Main CMakeLists.txt for project link_chameleon

cmake_minimum_required(VERSION 2.8)
project(LINK_CHAMELEON Fortran C CXX)

# Add extra cmake module path and initialize morse cmake modules
# --------------------------------------------------------------
set( MORSE_DISTRIB_DIR "" CACHE PATH "Directory of MORSE distribution")

if (MORSE_DISTRIB_DIR)
    set( MORSE_CMAKE_MODULE_DIR "${MORSE_DISTRIB_DIR}/cmake_modules/morse" CACHE PATH
         "Directory where to find MORSE CMake modules (cmake_modules/morse)")
    list(APPEND CMAKE_MODULE_PATH "${MORSE_CMAKE_MODULE_DIR}")
    list(APPEND CMAKE_MODULE_PATH "${MORSE_CMAKE_MODULE_DIR}/find")
    include(MorseInit)

    # Detect CHAMELEON
    #find_package(CHAMELEON REQUIRED COMPONENTS QUARK)
    find_package(CHAMELEON REQUIRED)
    link_directories(${CHAMELEON_LIBRARY_DIRS})
    include_directories(${CHAMELEON_INCLUDE_DIRS})

    # link_chameleon exe
    add_executable(link_chameleon link_chameleon.c)
    target_link_libraries(link_chameleon ${CHAMELEON_LIBRARIES})
else()
    message(STATUS "MORSE_DISTRIB_DIR is not set")
    message(STATUS "Please indicate where is located your MORSE distribution directory."
    " This is necessary to find cmake_modules.")
endif()

###
### END CMakeLists.txt
###
