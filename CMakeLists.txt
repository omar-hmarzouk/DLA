###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file CMakeLists.txt
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###
cmake_minimum_required(VERSION 2.8)
project(CHAMELEON Fortran C CXX)

# directly make an error if in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.\n"
   "Please create a build directory first and execute cmake configuration from "
   "this directory. Example: mkdir build && cd build && cmake ..")
endif()

# set project version number
set(CHAMELEON_VERSION_MAJOR 0)
set(CHAMELEON_VERSION_MINOR 9)
set(CHAMELEON_VERSION_PATCH 0)

# Add extra cmake module path and initialize morse cmake modules
# --------------------------------------------------------------
set( MORSE_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse )
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
list(APPEND CMAKE_MODULE_PATH ${MORSE_CMAKE_MODULE_PATH})
include(MorseInit)
include(GenPkgConfig)

#############################################
#                                           #
#        Compilation of CHAMELEON           #
#                                           #
#############################################


###############################################################################
# Parameters/Options #
######################

# Add define for Fortran Mangling (should be defined somewhere else)
# ------------------------------------------------------------------
add_definitions(-DADD_)


# Misc options
# ------------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define precision supported by CHAMELEON
# -----------------------------------------
set( RP_CHAMELEON_DICTIONNARY ${MORSE_CMAKE_MODULE_PATH}/precision_generator/subs.py )
set( RP_CHAMELEON_PRECISIONS  "s;d;c;z" )
include(RulesPrecisions)

# Options to select the runtime
# -----------------------------

# Create a list of possible runtime
set(CHAMELEON_SCHED_list "STARPU" "QUARK")

# Initially StarPU runtime is enabled
option(CHAMELEON_SCHED_STARPU
    "Enable StarPU scheduler as the default runtime
    (Conflict with other CHAMELEON_SCHED_* options)" ON)
option(CHAMELEON_SCHED_QUARK
    "Enable Quark scheduler as the default runtime
    (Conflict with other CHAMELEON_SCHED_* options)" OFF)

if (CHAMELEON_SCHED_QUARK)
    set(CHAMELEON_SCHED_STARPU OFF)
elseif (CHAMELEON_SCHED_STARPU)
    set(CHAMELEON_SCHED_QUARK OFF)
endif()

# Check that one, and only one, SCHED option is set to ON
# count number of runtime sets to ON
math(EXPR number_of_active_runtime 0)
foreach (runtime ${CHAMELEON_SCHED_list})
    if (CHAMELEON_SCHED_${runtime})
        math(EXPR number_of_active_runtime "${number_of_active_runtime}+1")
    endif()
endforeach()
if (NOT number_of_active_runtime STREQUAL 1)
    message(FATAL_ERROR
            "Number of active runtime is ${number_of_active_runtime}, "
            "the user should activate one (and only one) runtime. ")
endif()

# Warning if not coherent options
if (CHAMELEON_USE_MPI AND NOT CHAMELEON_SCHED_STARPU)
    message(WARNING "You have activated MPI but CHAMELEON_SCHED_STARPU is OFF.\n"
    "You should use StarPU Runtime system if you want to benefit from MPI.\n"
    "Use -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON at cmake \n"
    "configure to do so.")
endif()

if(CHAMELEON_USE_CUDA AND NOT CHAMELEON_SCHED_STARPU)
    message(WARNING "You have activated CUDA but CHAMELEON_SCHED_STARPU is OFF.\n"
    "You should use StarPU Runtime system if you want to benefit from GPUs.\n"
    "Use -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON at cmake \n"
    "configure to do so.")
endif()

if(CHAMELEON_USE_MAGMA AND NOT CHAMELEON_SCHED_STARPU)
    message(WARNING "You have activated MAGMA but CHAMELEON_SCHED_STARPU is OFF.\n"
    "You should use StarPU Runtime system if you want to benefit from GPUs.\n"
    "Use -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON at cmake \n"
    "configure to do so.")
endif()

if(CHAMELEON_USE_FXT AND NOT CHAMELEON_SCHED_STARPU)
    message(WARNING "You have activated FXT but CHAMELEON_SCHED_STARPU is OFF.\n"
    "You should use StarPU Runtime system if you want to benefit from FxT.\n"
    "Use -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON at cmake \n"
    "configure to do so.")
endif()
        
# Additional options
# ------------------

# Enable the distributed interface (allowed only when StarPU is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
cmake_dependent_option(CHAMELEON_USE_MPI
                    "Enable distributed memory through MPI" OFF
                    "CHAMELEON_SCHED_STARPU" OFF)

# Enable Cuda kernels if StarPU (only if StarPU runtime is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
cmake_dependent_option(CHAMELEON_USE_CUDA
                    "Enable CUDA kernels" OFF
                    "CHAMELEON_SCHED_STARPU" OFF)

# Enable MAGMA advanced kernels if CUDA is enabled
cmake_dependent_option(CHAMELEON_USE_MAGMA
                    "Enable MAGMA Cuda kernels" OFF
                    "CHAMELEON_USE_CUDA" OFF)

# Enable FXT if StarPU (only if StarPU runtime is enabled)
cmake_dependent_option(CHAMELEON_USE_FXT
                    "Enable FXT with STARPU" OFF
                    "CHAMELEON_SCHED_STARPU" OFF)

# Options to enable/disable testings and timings
# ----------------------------------------------
option(CHAMELEON_ENABLE_DOCS     "Enable documentation build"  ON)
option(CHAMELEON_ENABLE_EXAMPLE  "Enable examples build"       ON)
option(CHAMELEON_ENABLE_TESTING  "Enable testings build"       ON)
option(CHAMELEON_ENABLE_TIMING   "Enable timings build"        ON)

# Option to activate or not simulation mode (use Simgrid through StarPU)
# ----------------------------------------------------------------------
if(CHAMELEON_SCHED_STARPU)
    option(CHAMELEON_SIMULATION "Enable simulation mode using Simgrid through StarPU" OFF)
endif()

# Initially we need to generate files for different precisions
# TODO: use this option to avoid generate multiple precisions each time we launch cmake
#option(CHAMELEON_GEN_PREC "Generate source files precisions" ON)
#------------------------------------------------------------------------------


###############################################################################
# Look for dependencies #
#########################

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if( THREADS_FOUND )
    list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
set(M_LIBRARIES "")
if(UNIX OR WIN32)
    find_library(
        M_m_LIBRARY
        NAMES m
        )
    mark_as_advanced(M_m_LIBRARY)
    if (M_m_LIBRARY)
        list(APPEND M_LIBRARIES "${M_m_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${M_m_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find libm on your system."
            "Are you sure to a have a C compiler installed?")        
    endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library) 
# on Unix systems except Apple ones because it does not exist on it
set(RT_LIBRARIES "")
if(UNIX AND NOT APPLE)
    find_library(
        RT_rt_LIBRARY
        NAMES rt
        )
    mark_as_advanced(RT_rt_LIBRARY)
    if (RT_rt_LIBRARY)
        list(APPEND RT_LIBRARIES "${RT_rt_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${RT_rt_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find librt on your system")        
    endif()
endif()

# If simulation we don't enter in kernel functions so that we don't need to 
# link with concerned libraries
if(NOT CHAMELEON_SIMULATION)

    # CHAMELEON depends on CBLAS
    #----------------------------
    set(BLAS_DIR "" CACHE PATH "Root directory of BLAS library")
    find_package(BLASEXT REQUIRED)
    set(HAVE_BLAS ${BLAS_FOUND})

    set(CBLAS_DIR "" CACHE PATH "Root directory of CBLAS library")
    find_package(CBLAS REQUIRED)
    set(HAVE_CBLAS ${CBLAS_FOUND})
    if(CBLAS_FOUND)
        include_directories(${CBLAS_INCLUDE_DIRS})
    endif()

    #---------------------------------------------
    # CHAMELEON depends on LAPACKE (via coreblas)
    set(LAPACK_DIR "" CACHE PATH "Root directory of LAPACK library")
    find_package(LAPACKEXT REQUIRED)
    set(HAVE_LAPACK ${LAPACK_FOUND})

    # standalone version of lapacke seems useless for now
    # let the comment in case we meet some problems of non existing lapacke 
    # functions in lapack library such as mkl, acml, ...
    #set(LAPACKE_STANDALONE TRUE)
    set(LAPACKE_DIR "" CACHE PATH "Root directory of LAPACKE library")
    find_package(LAPACKE REQUIRED)
    set(HAVE_LAPACKE ${LAPACKE_FOUND})
    if(LAPACKE_FOUND)
        include_directories(${LAPACKE_INCLUDE_DIRS})
    endif()

    # CHAMELEON depends on TMG
    #--------------------------
    if (CHAMELEON_ENABLE_TESTING OR CHAMELEON_ENABLE_TIMING)
        set(TMG_DIR "" CACHE PATH "Root directory of TMG library")
        find_package(TMG REQUIRED)
        set(HAVE_TMG ${TMG_FOUND})
    endif()

    # CHAMELEON depends on CUDA/CUBLAS
    #----------------------------------
    if (CHAMELEON_USE_CUDA)

        find_package(CUDA REQUIRED)

        set(HAVE_CUDA ${CUDA_FOUND})
        if (CUDA_FOUND)
            add_definitions(-DCHAMELEON_USE_CUDA)
            if(CUDA_VERSION VERSION_LESS "3.0")
                set(CUDA_HOST_COMPILATION_CPP OFF)
            endif(CUDA_VERSION VERSION_LESS "3.0")
            set(CUDA_BUILD_EMULATION OFF)
            include_directories(${CUDA_INCLUDE_DIRS})
            list(APPEND EXTRA_LIBRARIES ${CUDA_LIBRARIES}) # need CUDA libs to link.
            list(APPEND EXTRA_SOURCES src/gpu_data.c)
            set(saved_include "${CMAKE_REQUIRED_INCLUDES}")
            set(saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
            set(CMAKE_REQUIRED_INCLUDES  "${CMAKE_REQUIRED_INCLUDES};${CUDA_INCLUDE_DIRS}")
            set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CUDA_LIBRARIES}")
            if(CUDA_VERSION VERSION_LESS "4.0")
                set(CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS 0)
            else()
                check_function_exists(cuDeviceCanAccessPeer CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS)
            endif()
            set(CMAKE_REQUIRED_INCLUDES  "${saved_include}")
            set(CMAKE_REQUIRED_LIBRARIES "${saved_libs}")
        endif (CUDA_FOUND)

        # CHAMELEON depends on MAGMA gpu kernels
        # call our cmake module to test (in cmake_modules)
        # change this call position if not appropriated
        #----------------------------------------------
        if ( CUDA_FOUND )
            if ( CHAMELEON_USE_MAGMA )
                set(MAGMA_DIR "" CACHE PATH "Root directory of MAGMA library")
                find_package(MAGMA REQUIRED COMPONENTS CBLAS LAPACK CUDA)

                set(HAVE_MAGMA ${MAGMA_FOUND})
                if ( MAGMA_FOUND )
                    add_definitions(-DHAVE_CUBLAS)
                    add_definitions(-DCHAMELEON_USE_MAGMA)
                    list(APPEND EXTRA_LIBRARIES ${CUDA_CUBLAS_LIBRARIES})
                    include_directories(${MAGMA_INCLUDE_DIRS})
                endif( MAGMA_FOUND )
            endif()
        endif()

    endif(CHAMELEON_USE_CUDA)
    
    # CHAMELEON depends on MPI
    #--------------------------
    if (CHAMELEON_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure 
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        find_package(MPI REQUIRED)

        set(HAVE_MPI ${MPI_C_FOUND})
        if (MPI_C_FOUND)
            add_definitions(-DCHAMELEON_USE_MPI)
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            # tests for intel mpi
            #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
            #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
        endif (MPI_C_FOUND)

    endif (CHAMELEON_USE_MPI)
    
else (NOT CHAMELEON_SIMULATION)

    message("${BoldBlue}Simulation mode of CHAMELEON is activated (CHAMELEON_SIMULATION=ON)."
        "\nThis mode allows you to simulate execution of algorithms with StarPU compiled with SimGrid."
        "\nTo do so, we provide some perfmodels in the simucore/perfmodels/ directory of CHAMELEON sources."
        "\nThe resulting directory where to find perfmodels is: ${CMAKE_CURRENT_SOURCE_DIR}/simucore/perfmodels."
        "\nTo use these perfmodels, please set your STARPU_HOME environment variable to ${CMAKE_CURRENT_SOURCE_DIR}/simucore/perfmodels."
        "\nFinally, you need to set your STARPU_HOSTNAME environment variable to the name of the machine to simulate."
        "For example: STARPU_HOSTNAME=mirage."
        "\nNote that only POTRF kernels with block sizes of 320 or 960 (simple and double precision) on mirage machine are available for now."
        "\nDatabase of models is subject to change, it should be enrich in a near future.${ColourReset}")

    if (CHAMELEON_USE_CUDA)
        set(CHAMELEON_USE_CUDA OFF)
        message("${BoldBlue}CHAMELEON_USE_CUDA is set to ON but we turn it OFF."
        "Because we are compiling the simulation mode (CHAMELEON_SIMULATION=ON),"
        "we do not need to link with CUDA.${ColourReset}")
    endif ()
    if (CHAMELEON_USE_MAGMA)
        set(CHAMELEON_USE_MAGMA OFF)
        message("${BoldBlue}CHAMELEON_USE_MAGMA is set to ON but we turn it OFF."
        "Because we are compiling the simulation mode (CHAMELEON_SIMULATION=ON),"
        "we do not need to link with MAGMA.${ColourReset}")
    endif ()
    if (CHAMELEON_USE_MPI)
        set(CHAMELEON_USE_MPI OFF)
        message("${BoldBlue}CHAMELEON_USE_MPI is set to ON but we turn it OFF."
        "Because we are compiling the simulation mode (CHAMELEON_SIMULATION=ON),"
        "we do not need to link with MPI.${ColourReset}")
    endif ()
    if (CHAMELEON_ENABLE_TESTING)
        set(CHAMELEON_ENABLE_TESTING OFF)
        message("${BoldBlue}CHAMELEON_ENABLE_TESTING is set to ON but we turn it OFF."
        "Because we are compiling the simulation mode (CHAMELEON_SIMULATION=ON),"
        "there is no sense in compiling testing drivers that are used to check"
        "numerical correctness of algorithms and kernels.${ColourReset}")
    endif ()

endif(NOT CHAMELEON_SIMULATION)


# CHAMELEON depends on a runtime
# -------------------------------

# TODO: discuss on this with Manu
# try to find starpu runtime
if( CHAMELEON_SCHED_STARPU )

    set(CHAMELEON_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")

    # create list of components in order to make a single call to find_package(starpu...)
    if(NOT CHAMELEON_SIMULATION)
        set(STARPU_COMPONENT_LIST "HWLOC")
        if(CHAMELEON_USE_MPI)
            list(APPEND STARPU_COMPONENT_LIST "MPI")
        endif()
        if(CHAMELEON_USE_CUDA)
            list(APPEND STARPU_COMPONENT_LIST "CUDA")
        endif()
    endif()
    if(CHAMELEON_USE_FXT)
        list(APPEND STARPU_COMPONENT_LIST "FXT")
    endif()

    set(HWLOC_DIR "" CACHE PATH "Root directory of HWLOC library")
    set(FXT_DIR "" CACHE PATH "Root directory of FXT library")
    set(STARPU_DIR "" CACHE PATH "Root directory of STARPU library")
    find_package(STARPU ${CHAMELEON_STARPU_VERSION} REQUIRED
                 COMPONENTS ${STARPU_COMPONENT_LIST})

    # Is this useful??
    set(HAVE_STARPU ${STARPU_FOUND})

    # Add definition and include_dir if found
    if(STARPU_FOUND)
        add_definitions(-DCHAMELEON_SCHED_STARPU) # TODO: remove that
        include_directories(${STARPU_INCLUDE_DIRS})
        set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS}")
        foreach(libdir ${STARPU_LIBRARY_DIRS})
            list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
        endforeach()
        string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
        if(CHAMELEON_USE_FXT)
            # check if fxt profiling is accessible in starpu and activate it in chameleon
            set(CMAKE_REQUIRED_LIBRARIES "${STARPU_SHM_LIBRARIES}")
            unset(STARPU_FXT_START_PROFILING_FOUND CACHE)
            check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING_FOUND)
            if ( STARPU_FXT_START_PROFILING_FOUND )
                message(STATUS "Set HAVE_STARPU_FXT_PROFILING - Activate FxT profiling through StarPU")
                add_definitions(-DHAVE_STARPU_FXT_PROFILING)
            endif()
        endif()
        if (CHAMELEON_USE_MPI AND STARPU_MPI_LIBRARIES)
            # Check if a specific function exist
            set(CMAKE_REQUIRED_LIBRARIES "${STARPU_MPI_LIBRARIES}")
            unset(STARPU_MPI_DATA_REGISTER_FOUND CACHE)
            check_function_exists(starpu_mpi_data_register STARPU_MPI_DATA_REGISTER_FOUND)
            if ( STARPU_MPI_DATA_REGISTER_FOUND )
                add_definitions(-DHAVE_STARPU_MPI_DATA_REGISTER)
            endif()
        endif()
    endif()

endif()

# try to find quark runtime
if( CHAMELEON_SCHED_QUARK )

    set(HWLOC_DIR "" CACHE PATH "Root directory of HWLOC library")
    set(QUARK_DIR "" CACHE PATH "Root directory of QUARK library")
    find_package(QUARK REQUIRED COMPONENTS HWLOC)
    set(HAVE_QUARK ${QUARK_FOUND})

    # Add definition and include_dir if found
    if(QUARK_FOUND)
        add_definitions(-DCHAMELEON_SCHED_QUARK) # TODO: remove that
        include_directories(${QUARK_INCLUDE_DIRS})
    endif()

endif()

#------------------------------------------------------------------------------
# Save all dependencies (all required links) 
set(CHAMELEON_DEP -lchameleon)
if(CHAMELEON_SCHED_STARPU)
    if ( CHAMELEON_USE_MPI )
        list(APPEND CHAMELEON_DEP
          -lchameleon_starpu
          ${STARPU_MPI_LIBRARIES}
        )
    else()
        list(APPEND CHAMELEON_DEP
          -lchameleon_starpu
          ${STARPU_SHM_LIBRARIES}
        )
    endif()
elseif(CHAMELEON_SCHED_QUARK)
    list(APPEND CHAMELEON_DEP
      -lchameleon_quark
      ${QUARK_LIBRARIES}
      )
    link_directories(${QUARK_LIBRARY_DIRS})
endif()
if(NOT CHAMELEON_SIMULATION)
    if(CHAMELEON_USE_CUDA)
        list(APPEND CHAMELEON_DEP
        ${CUDA_LIBRARIES}
        )
        link_directories(${CUDA_LIBRARY_DIRS})
    endif()
    if(CHAMELEON_USE_MAGMA)
        list(APPEND CHAMELEON_DEP
        ${MAGMA_LIBRARIES}
        )
        link_directories(${MAGMA_LIBRARY_DIRS})
    endif()

    list(APPEND CHAMELEON_DEP
    -lcoreblas
    ${LAPACKE_LIBRARIES}
    ${CBLAS_LIBRARIES}    
    ${LAPACK_SEQ_LIBRARIES}
    ${BLAS_SEQ_LIBRARIES}
    ${HWLOC_LIBRARIES}
    ${EXTRA_LIBRARIES}
    )
else()
    list(APPEND CHAMELEON_DEP
    -lcoreblas
    -lsimulapacke    
    -lsimucblas
    ${HWLOC_LIBRARIES}
    ${EXTRA_LIBRARIES}
    )
endif()
list(REMOVE_DUPLICATES CHAMELEON_DEP)

# Check for the subdirectories
# ----------------------------

# include headers
# ---------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/coreblas/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/coreblas/include)

###############################################################################
# simucore libraries #
######################

if ( CHAMELEON_SIMULATION )
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/simucore/simucblas)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/simucore/simucblas)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/simucore/simulapacke)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/simucore/simulapacke)
    add_subdirectory(simucore)
    add_definitions(-DCHAMELEON_SIMULATION)
endif()
#------------------------------------------------------------------------------


###############################################################################
# Coreblas library (kernels for CPUs) #
#######################################

add_subdirectory(coreblas)
#------------------------------------------------------------------------------


###############################################################################
# Main library #
################

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/control)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/control)
add_subdirectory(include)
add_subdirectory(control)
add_subdirectory(compute)
#------------------------------------------------------------------------------


###############################################################################
# Runtime library #
###################

add_subdirectory(runtime)
#------------------------------------------------------------------------------


###############################################################################
# Executables and Tests #
#########################

enable_testing()
include(CTest)
#enable_testing()

# Necessary to compile executables (temporary)
# should use SystemDetection.cmake module?
string(REGEX MATCH ".*ifort$" _match_ifort ${CMAKE_Fortran_COMPILER})
if(_match_ifort)
    #message(STATUS "Add -nofor_main to the Fortran linker")
    set(CMAKE_Fortran_FLAGS "-nofor_main")
endif(_match_ifort)

# Testing executables
if(CHAMELEON_ENABLE_EXAMPLE)
  add_subdirectory(example)
endif(CHAMELEON_ENABLE_EXAMPLE)

# Testing executables
if(CHAMELEON_ENABLE_TESTING)
  add_subdirectory(testing)
endif(CHAMELEON_ENABLE_TESTING)

# Timing executables
if(CHAMELEON_ENABLE_TIMING)
  add_subdirectory(timing)
endif(CHAMELEON_ENABLE_TIMING)
#------------------------------------------------------------------------------


###############################################################################
# Documentation #
#################

if(CHAMELEON_ENABLE_DOCS)
    add_subdirectory(docs)
endif()
#------------------------------------------------------------------------------


###############################################################################
# Config files (pkg_config) #
#############################

# Create file morse_starpu.pc
# ---------------------------
if(CHAMELEON_SCHED_STARPU)
    generate_pkgconfig_file("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/chameleon-starpu.pc.in")
endif()

# # Create file morse_quark.pc
# # --------------------------
if(CHAMELEON_SCHED_QUARK)
    generate_pkgconfig_file("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/chameleon-quark.pc.in")
endif()
#------------------------------------------------------------------------------


###############################################################################
# Print Options #
#################
include(PrintOpts)

###
### END CMakeLists.txt
###
