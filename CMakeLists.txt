###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
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
set(MORSE_VERSION_MAJOR 1)
set(MORSE_VERSION_MINOR 0)
set(MORSE_VERSION_MICRO 0)
# set subproject version number
set(CHAMELEON_VERSION_MAJOR 0)
set(CHAMELEON_VERSION_MINOR 9)
set(CHAMELEON_VERSION_MICRO 1)

set( MORSE_DISTRIB_DIR "" CACHE PATH "Directory of MORSE distribution")

if (MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/")

    # Add extra cmake module path and initialize morse cmake modules
    # --------------------------------------------------------------
    if(MORSE_DISTRIB_DIR)
        list(APPEND CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/cmake_modules)
        list(APPEND CMAKE_MODULE_PATH "${MORSE_DISTRIB_DIR}/cmake_modules/morse")
        set(MORSE_CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/cmake_modules/morse )
    elseif(EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/")
        list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
        list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/morse")
        set(MORSE_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse )
    endif()
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

    set(CHAMELEON_DEFINITIONS_LIST "")

    # Set the RPATH config
    # --------------------

    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    # the RPATH to be used when installing
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

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
    if (CHAMELEON_SCHED_STARPU)
        message("-- ${BoldGreen}CHAMELEON_SCHED_STARPU is set to ON: CHAMELEON uses StarPU runtime${ColourReset}")
        message("-- ${BoldGreen}To use CHAMELEON with QUARK runtime: set CHAMELEON_SCHED_QUARK"
        "to ON (CHAMELEON_SCHED_STARPU will be disabled)${ColourReset}")
    elseif(CHAMELEON_SCHED_QUARK)
        message("-- ${BoldGreen}CHAMELEON_SCHED_QUARK is set to ON: CHAMELEON uses QUARK runtime${ColourReset}")
        message("-- ${BoldGreen}To use CHAMELEON with StarPU runtime: set CHAMELEON_SCHED_STARPU"
        "to ON (CHAMELEON_SCHED_QUARK will be disabled)${ColourReset}")
    endif()

    if (CHAMELEON_USE_MAGMA)
        set(CHAMELEON_USE_CUDA ON)
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
    if (NOT CHAMELEON_USE_MPI)
        message("-- ${BoldGreen}CHAMELEON_USE_MPI is set to OFF, turn it ON to use MPI (only with StarPU)${ColourReset}")
    endif()

    # Enable Cuda kernels if StarPU (only if StarPU runtime is enabled)
    # TODO: Default should be changed to ON/OFF when it will be ok
    cmake_dependent_option(CHAMELEON_USE_CUDA
                        "Enable CUDA kernels" OFF
                        "CHAMELEON_SCHED_STARPU" OFF)
    if (NOT CHAMELEON_USE_CUDA)
        message("-- ${BoldGreen}CHAMELEON_USE_CUDA is set to OFF, turn it ON to use CUDA (only with StarPU)${ColourReset}")
    endif()

    # Enable MAGMA advanced kernels if CUDA is enabled
    cmake_dependent_option(CHAMELEON_USE_MAGMA
                        "Enable MAGMA Cuda kernels" OFF
                        "CHAMELEON_USE_CUDA" OFF)
    if (NOT CHAMELEON_USE_MAGMA)
        message("-- ${BoldGreen}CHAMELEON_USE_MAGMA is set to OFF, turn it ON to use MAGMA (only with StarPU)${ColourReset}")
    endif()

    # Enable FXT if StarPU (only if StarPU runtime is enabled)
    cmake_dependent_option(CHAMELEON_USE_FXT
                        "Enable FXT with STARPU" OFF
                        "CHAMELEON_SCHED_STARPU" OFF)
    if (NOT CHAMELEON_USE_FXT)
        message("-- ${BoldGreen}CHAMELEON_USE_FXT is set to OFF, turn it ON to use FxT (only with StarPU)${ColourReset}")
    endif()

    # Options to enable/disable testings and timings
    # ----------------------------------------------
    option(CHAMELEON_ENABLE_DOCS     "Enable documentation build"  OFF)
    if (CHAMELEON_ENABLE_DOCS)
        message("-- ${BoldGreen}CHAMELEON_ENABLE_DOCS is set to ON, turn it OFF to avoid building docs${ColourReset}")
    endif()
    option(CHAMELEON_ENABLE_EXAMPLE  "Enable examples build"       ON)
    if (CHAMELEON_ENABLE_EXAMPLE)
        message("-- ${BoldGreen}CHAMELEON_ENABLE_EXAMPLE is set to ON, turn it OFF to avoid building examples${ColourReset}")
    endif()
    option(CHAMELEON_ENABLE_TESTING  "Enable testings build"       ON)
    if (CHAMELEON_ENABLE_TESTING)
        message("-- ${BoldGreen}CHAMELEON_ENABLE_TESTING is set to ON, turn it OFF to avoid building testing${ColourReset}")
    endif()
    option(CHAMELEON_ENABLE_TIMING   "Enable timings build"        ON)
    if (CHAMELEON_ENABLE_TIMING)
        message("-- ${BoldGreen}CHAMELEON_ENABLE_TIMING is set to ON, turn it OFF to avoid building timing${ColourReset}")
    endif()

    # Option to activate or not simulation mode (use Simgrid through StarPU)
    # ----------------------------------------------------------------------
    if(CHAMELEON_SCHED_STARPU)
        option(CHAMELEON_SIMULATION "Enable simulation mode using Simgrid through StarPU" OFF)
        if (NOT CHAMELEON_SIMULATION)
            message("-- ${BoldGreen}CHAMELEON_SIMULATION is set to OFF, turn it ON to use"
                "SIMULATION mode (only with StarPU compiled with SimGrid)${ColourReset}")
        endif()
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
        #---------------------------
        find_package(CBLAS COMPONENTS BLASEXT)
        if(BLAS_FOUND)
            if (BLAS_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
            endif()
            if(BLAS_LINKER_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
            endif()
        else()
            message(FATAL_ERROR "BLAS library has not been found")
        endif()
        if(CBLAS_FOUND)
            include_directories(${CBLAS_INCLUDE_DIRS})
            if(CBLAS_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS}")
            endif()
        else()
            if(MORSE_VERBOSE_FIND_PACKAGE)
                if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
                    if (NOT CBLAS_cblas.h_DIRS)
                        Print_Find_Header_Status(cblas cblas.h)
                    endif ()
                    if (NOT CBLAS_cblas_LIBRARY)
                        Print_Find_Library_Status(cblas libcblas)
                    endif ()
                endif()
            else()
                message(WARNING "CBLAS library has not been found and MORSE_VERBOSE_FIND_PACKAGE is set to OFF"
                "Try to activate MORSE_VERBOSE_FIND_PACKAGE option (-DMORSE_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
            endif()
            message(FATAL_ERROR "A CBLAS library is required but has not been found")
        endif()

        # CHAMELEON depends on LAPACKE
        #-----------------------------
        # standalone version of lapacke seems useless for now
        # let the comment in case we meet some problems of non existing lapacke
        # functions in lapack library such as mkl, acml, ...
        #set(LAPACKE_STANDALONE TRUE)
        find_package(LAPACKE COMPONENTS LAPACKEXT)
        if(LAPACK_FOUND AND LAPACK_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${LAPACK_LIBRARY_DIRS}")
        else()
            message(FATAL_ERROR "A LAPACK library is required but has not been found")
        endif()
        if(LAPACKE_FOUND)
            include_directories(${LAPACKE_INCLUDE_DIRS})
            if(LAPACKE_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS}")
            endif()
            if(LAPACKE_LINKER_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${LAPACKE_LINKER_FLAGS}")
            endif()
        else()
            if(MORSE_VERBOSE_FIND_PACKAGE)
                if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
                    if (NOT LAPACKE_lapacke.h_DIRS)
                        Print_Find_Header_Status(lapacke lapacke.h)
                    endif ()
                    if (NOT LAPACKE_lapacke_LIBRARY)
                        Print_Find_Library_Status(lapacke liblapacke)
                    endif ()
                endif()
            else()
                message(WARNING "LAPACKE library has not been found and MORSE_VERBOSE_FIND_PACKAGE is set to OFF"
                "Try to activate MORSE_VERBOSE_FIND_PACKAGE option (-DMORSE_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
            endif()
            message(FATAL_ERROR "A LAPACKE library is required but has not been found")
        endif()

        # CHAMELEON depends on TMG
        #-------------------------
        if (CHAMELEON_ENABLE_TESTING OR CHAMELEON_ENABLE_TIMING)
            find_package(TMG)
            if(TMG_FOUND AND TMG_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${TMG_LIBRARY_DIRS}")
            else()
                message(FATAL_ERROR "A TMG library is required but has not been found")
            endif()
            if(TMG_LINKER_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${TMG_LINKER_FLAGS}")
            endif()
        endif()

        # CHAMELEON depends on CUDA/CUBLAS
        #----------------------------------
        if (CHAMELEON_USE_CUDA)

            find_package(CUDA REQUIRED)

            if (CUDA_FOUND)
                message("-- ${Blue}Add definition CHAMELEON_USE_CUDA"
                " - Activate CUDA in Chameleon${ColourReset}")
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
                # Add cublas if found
                if (CUDA_CUBLAS_LIBRARIES)
                    message("-- ${Blue}Add definition CHAMELEON_USE_CUBLAS"
                    " - Use GPU kernels from cuBLAS${ColourReset}")
                    add_definitions(-DCHAMELEON_USE_CUBLAS)
                    list(APPEND EXTRA_LIBRARIES ${CUDA_CUBLAS_LIBRARIES})
                endif()
                #FindHeader(CUDA cublas_v2.h)
                #if(CUDA_cublas_v2.h_DIRS)
                    #add_definitions(-DCHAMELEON_USE_CUBLAS_V2)
                #endif()
            endif (CUDA_FOUND)

            # CHAMELEON depends on MAGMA gpu kernels
            # call our cmake module to test (in cmake_modules)
            # change this call position if not appropriated
            #----------------------------------------------
            if ( CUDA_FOUND AND CHAMELEON_USE_MAGMA )
                set(CHAMELEON_MAGMA_VERSION "1.4" CACHE STRING "oldest MAGMA version desired")
                find_package(MAGMA ${CHAMELEON_MAGMA_VERSION})
                if ( MAGMA_FOUND )
                    message("-- ${Blue}Add definition CHAMELEON_USE_MAGMA"
                    " - Use GPU kernels from MAGMA${ColourReset}")
                    add_definitions(-DCHAMELEON_USE_MAGMA)
                    include_directories(${MAGMA_INCLUDE_DIRS})
                    if(MAGMA_LIBRARY_DIRS)
                        # the RPATH to be used when installing
                        list(APPEND CMAKE_INSTALL_RPATH "${MAGMA_LIBRARY_DIRS}")
                    endif()
                else( MAGMA_FOUND )
                    if(MORSE_VERBOSE_FIND_PACKAGE)
                        if (NOT MAGMA_magma.h_DIRS)
                            Print_Find_Header_Status(magma magma.h)
                        endif ()
                        if (NOT MAGMA_magma_LIBRARY)
                            Print_Find_Library_Status(magma libmagma)
                        endif ()
                    else()
                        message(WARNING "MAGMA library has not been found and MORSE_VERBOSE_FIND_PACKAGE is set to OFF"
                        "Try to activate MORSE_VERBOSE_FIND_PACKAGE option (-DMORSE_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
                    endif()
                    message(FATAL_ERROR "MAGMA library is required but has not been found")
                endif( MAGMA_FOUND )
            endif()

        endif(CHAMELEON_USE_CUDA)

        # CHAMELEON depends on MPI
        #-------------------------
        if (CHAMELEON_USE_MPI)

            # allows to use an external mpi compilation by setting compilers with
            # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
            # at cmake configure
            if(NOT MPI_C_COMPILER)
                set(MPI_C_COMPILER mpicc)
            endif()
            find_package(MPI REQUIRED)

            if (MPI_C_FOUND)
                message("-- ${Blue}Add definition CHAMELEON_USE_MPI"
                " - Activate MPI in Chameleon${ColourReset}")
                add_definitions(-DCHAMELEON_USE_MPI)
                list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
                include_directories( ${MPI_C_INCLUDE_PATH} )
                # tests for intel mpi
                #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
                #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
                if(MPI_C_LINK_FLAGS)
                    list(APPEND CMAKE_EXE_LINKER_FLAGS "${MPI_C_LINK_FLAGS}")
                endif()
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

        find_package(STARPU ${CHAMELEON_STARPU_VERSION} REQUIRED
                     COMPONENTS ${STARPU_COMPONENT_LIST})

        # Add definition and include_dir if found
        if ( STARPU_FOUND )
            message("-- ${Blue}Add definition CHAMELEON_SCHED_STARPU"
            " - Activate StarPU in Chameleon${ColourReset}")
            add_definitions(-DCHAMELEON_SCHED_STARPU) # TODO: remove that
            include_directories(${STARPU_INCLUDE_DIRS})
            if(STARPU_LINKER_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARPU_LINKER_FLAGS}")
            endif()
            set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS}")
            foreach(libdir ${STARPU_LIBRARY_DIRS})
                list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
            endforeach()
            string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
            set(CMAKE_REQUIRED_LIBRARIES "${STARPU_SHM_LIBRARIES}")
            check_function_exists(starpu_data_idle_prefetch_on_node STARPU_IDLE_PREFETCH_FOUND)
            if ( STARPU_IDLE_PREFETCH_FOUND )
                add_definitions(-DCHAMELEON_USE_STARPU_IDLE_PREFETCH)
            endif()
            check_function_exists(starpu_data_wont_use STARPU_DATA_WONT_USE_FOUND)
            if ( STARPU_DATA_WONT_USE_FOUND )
                add_definitions(-DCHAMELEON_USE_STARPU_DATA_WONT_USE)
            endif()
            if(CHAMELEON_USE_FXT)
                # check if fxt profiling is accessible in starpu and activate it in chameleon
                set(CMAKE_REQUIRED_LIBRARIES "${STARPU_SHM_LIBRARIES}")
                unset(STARPU_FXT_START_PROFILING_FOUND CACHE)
                check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING_FOUND)
                if ( STARPU_FXT_START_PROFILING_FOUND )
                    message("-- ${Blue}Add definition HAVE_STARPU_FXT_PROFILING"
                    " - Activate FxT profiling through StarPU${ColourReset}")
                    add_definitions(-DHAVE_STARPU_FXT_PROFILING)
                else()
                    message("-- ${Red}Looking for starpu with fxt"
                    " - starpu_fxt_start_profiling() test fails in StarPU${ColourReset}")
                    message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
                endif()
            endif()
            if (CHAMELEON_USE_MPI AND STARPU_MPI_LIBRARIES)
                # Check if a specific function exist
                set(CMAKE_REQUIRED_LIBRARIES "${STARPU_MPI_LIBRARIES}")
                unset(STARPU_MPI_DATA_REGISTER_FOUND CACHE)
                check_function_exists(starpu_mpi_data_register_comm STARPU_MPI_DATA_REGISTER_FOUND)
                if ( STARPU_MPI_DATA_REGISTER_FOUND )
                    message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_REGISTER - Activate"
                    "use of starpu_mpi_data_register() in Chameleon with StarPU${ColourReset}")
                    add_definitions(-DHAVE_STARPU_MPI_DATA_REGISTER)
                else()
                    message("-- ${Red}Looking for starpu with starpu_mpi_data_register"
                    " - starpu_mpi_data_register() test fails in StarPU${ColourReset}")
                    message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
                endif()
            endif()
            if(HWLOC_FOUND AND HWLOC_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${HWLOC_LIBRARY_DIRS}")
            endif()
            if(FXT_FOUND AND FXT_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${FXT_LIBRARY_DIRS}")
            endif()
            if(STARPU_FOUND AND STARPU_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${STARPU_LIBRARY_DIRS}")
            endif()
        else ( STARPU_FOUND )
            if(MORSE_VERBOSE_FIND_PACKAGE)
                if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                    if (NOT HWLOC_hwloc.h_DIRS)
                        Print_Find_Header_Status(hwloc hwloc.h)
                    endif ()
                    if (NOT HWLOC_hwloc_LIBRARY)
                        Print_Find_Library_Status(hwloc libhwloc)
                    endif ()
                endif()
                if(CHAMELEON_USE_FXT AND (NOT FXT_FOUND OR NOT FXT_LIBRARIES))
                    if (NOT FXT_fxt.h_DIRS)
                        Print_Find_Header_Status(fxt fxt.h)
                    endif ()
                    if (NOT FXT_fxt_LIBRARY)
                        Print_Find_Library_Status(fxt libfxt)
                    endif ()
                endif()
                if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
                    ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
                  )
                    foreach(starpu_hdr ${STARPU_hdrs_to_find})
                        if (NOT STARPU_${starpu_hdr}_INCLUDE_DIRS)
                            Print_Find_Header_Status(starpu ${starpu_hdr})
                        endif ()
                    endforeach()
                    if(STARPU_VERSION_STRING)
                        foreach(starpu_lib ${STARPU_libs_to_find})
                            if (NOT STARPU_${starpu_lib}_LIBRARY)
                                Print_Find_Library_Status(starpu ${starpu_lib})
                            endif ()
                        endforeach()
                    endif()
                endif()
            else(MORSE_VERBOSE_FIND_PACKAGE)
                message(WARNING "StarPU library has not been found and MORSE_VERBOSE_FIND_PACKAGE is set to OFF"
                "Try to activate MORSE_VERBOSE_FIND_PACKAGE option (-DMORSE_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
            endif(MORSE_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                message(FATAL_ERROR "hwloc library is required but has not been found")
            endif()
            if(CHAMELEON_USE_FXT AND (NOT FXT_FOUND OR NOT FXT_LIBRARIES))
                message(FATAL_ERROR "FxT library is required but has not been found")
            endif()
            if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
                ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
              )
                message(FATAL_ERROR "StarPU library is required but has not been found")
            endif()
        endif ( STARPU_FOUND )

    endif( CHAMELEON_SCHED_STARPU )

    # try to find quark runtime
    if( CHAMELEON_SCHED_QUARK )

        find_package(QUARK COMPONENTS HWLOC)

        # Add definition and include_dir if found
        if(QUARK_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_SCHED_QUARK"
            " - Activate QUARK in Chameleon${ColourReset}")
            add_definitions(-DCHAMELEON_SCHED_QUARK) # TODO: remove that
            include_directories(${QUARK_INCLUDE_DIRS})
            if(QUARK_LIBRARY_DIRS)
                # the RPATH to be used when installing
                list(APPEND CMAKE_INSTALL_RPATH "${QUARK_LIBRARY_DIRS}")
            endif()
        else(QUARK_FOUND)
            if(MORSE_VERBOSE_FIND_PACKAGE)
                if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                    if (NOT HWLOC_hwloc.h_DIRS)
                        Print_Find_Header_Status(hwloc hwloc.h)
                    endif ()
                    if (NOT HWLOC_hwloc_LIBRARY)
                        Print_Find_Library_Status(hwloc libhwloc)
                    endif ()
                endif()
                if (NOT QUARK_quark.h_DIRS)
                    Print_Find_Header_Status(quark quark.h)
                endif ()
                if (NOT QUARK_quark_LIBRARY)
                    Print_Find_Library_Status(quark libquark)
                endif ()
            else()
                message(WARNING "QUARK library has not been found and MORSE_VERBOSE_FIND_PACKAGE is set to OFF"
                "Try to activate MORSE_VERBOSE_FIND_PACKAGE option (-DMORSE_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
            endif()
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                message(FATAL_ERROR "hwloc library is required but has not been found")
            endif()
            message(FATAL_ERROR "QUARK library is required but has not been found")
        endif(QUARK_FOUND)

    endif()

    list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
    string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

    # Add define for Fortran Mangling (should be defined somewhere else)
    # ------------------------------------------------------------------
    message("-- ${Blue}Add definition ADD_"
            " - For Fortran mangling${ColourReset}")
    add_definitions(-DADD_)

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
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

    ###############################################################################
    # simucore libraries #
    ######################

    if ( CHAMELEON_SIMULATION )
        add_subdirectory(simucore)
        message("-- ${Blue}Add definition CHAMELEON_SIMULATION"
        " - Activate simulation mode in Chameleon (to use with StarPU+SimGrid)${ColourReset}")
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
        list(APPEND CMAKE_Fortran_FLAGS "-nofor_main")
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
    get_directory_property( CHAMELEON_DEFINITIONS_LIST DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
    include(PrintOpts)

    if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files")
        message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} which is a default system path."
        "You may want to change it: set the CMAKE_INSTALL_PREFIX variable to do so${ColourReset}")
    else()
        message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}${ColourReset}")
    endif()

else(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/")

    message(STATUS "MORSE_DISTRIB_DIR is not set")
    message(STATUS "Please indicate where is located your MORSE distribution directory."
    " This is necessary to find cmake_modules.")
    message(STATUS "If you are a MORSE developer, you normally have access to the svn repository of the MORSE_DISTRIB_DIR: "
    "svn checkout svn+ssh://yourlogin@scm.gforge.inria.fr/svnroot/morse/trunk/morse_distrib")

endif(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/")

###
### END CMakeLists.txt
###
