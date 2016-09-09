###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# @file CMakeLists.c
#
#  MORSE example routines
#  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
#  University of Bordeaux, Bordeaux INP
#
#  @version 1.0.0
#  @author Florent Pruvost
#  @date 2016-08-23
#
###


include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# list of sources
set(OOC_SOURCES
    out_of_core.c
   )

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_ooc)
list(APPEND libs_for_ooc
     chameleon
     chameleon_starpu
     ${STARPU_LIBRARIES_DEP}
)
link_directories(${STARPU_LIBRARY_DIRS})


if(NOT CHAMELEON_SIMULATION)

    if(CHAMELEON_USE_CUDA OR CHAMELEON_USE_MAGMA)
        list(APPEND libs_for_ooc
        cudablas)
    endif()
    if(CHAMELEON_USE_CUDA)
        list(APPEND libs_for_ooc
             ${CUDA_LIBRARIES}
        )
        link_directories(${CUDA_LIBRARY_DIRS})
    endif()
    if(CHAMELEON_USE_MAGMA)
        list(APPEND libs_for_ooc
             ${MAGMA_LIBRARIES}
        )
        link_directories(${MAGMA_LIBRARY_DIRS})
    endif()

    list(APPEND libs_for_ooc
         coreblas
         ${LAPACKE_LIBRARIES}
         ${CBLAS_LIBRARIES}
         ${LAPACK_SEQ_LIBRARIES}
         ${BLAS_SEQ_LIBRARIES}
         ${HWLOC_LIBRARIES}
         ${EXTRA_LIBRARIES}
    )

    link_directories(${LAPACKE_LIBRARY_DIRS})
    link_directories(${LAPACK_LIBRARY_DIRS})
    link_directories(${CBLAS_LIBRARY_DIRS})
    link_directories(${BLAS_LIBRARY_DIRS})

else()

    list(APPEND libs_for_ooc
         coreblas
         simulapacke
         simucblas
         ${HWLOC_LIBRARIES}
         ${EXTRA_LIBRARIES}
    )

endif()

link_directories(${HWLOC_LIBRARY_DIRS})


# message(STATUS "libs examples: ${libs_for_ooc}")
foreach(_ooc ${OOC_SOURCES})
    get_filename_component(_name_exe ${_ooc} NAME_WE)
    add_executable(${_name_exe} ${_ooc})
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -undefined dynamic_lookup" )
    endif()
    target_link_libraries(${_name_exe} ${libs_for_ooc})
    install(TARGETS ${_name_exe}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/example/out_of_core)
endforeach()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
