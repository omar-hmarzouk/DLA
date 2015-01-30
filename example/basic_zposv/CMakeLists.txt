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
# @file CMakeLists.c
#
#  MORSE example routines
#  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
#  University of Bordeaux, Bordeaux INP
#
#  @version 1.0.0
#  @author Florent Pruvost
#  @date 2014-10-13
#
###

set(EXAMPLES_HDR
    basic_posv.h
    posv_morse_functions.h
    posv_users_functions.h
   )
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Generate the morse testing sources for all possible precisions
# --------------------------------------------------------------
set(EXAMPLES "")
set(ZSRC
    zposv_morse_functions.c
    zposv_users_functions.c
   )

precisions_rules_py(EXAMPLES "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_examples)
list(APPEND libs_for_examples
     chameleon
)

if(CHAMELEON_SCHED_STARPU)
    if ( CHAMELEON_USE_MPI )
        list(APPEND libs_for_examples
             chameleon_starpu
             ${STARPU_MPI_LIBRARIES}
        )
    else()
        list(APPEND libs_for_examples
             chameleon_starpu
             ${STARPU_SHM_LIBRARIES}
        )
    endif()
    link_directories(${STARPU_LIBRARY_DIRS})
elseif(CHAMELEON_SCHED_QUARK)
    list(APPEND libs_for_examples
         chameleon_quark
         ${QUARK_LIBRARIES}
      )
    link_directories(${QUARK_LIBRARY_DIRS})
endif()

if(NOT CHAMELEON_SIMULATION)

    if(CHAMELEON_USE_CUDA)
        list(APPEND libs_for_examples
             ${CUDA_LIBRARIES}
        )
        link_directories(${CUDA_LIBRARY_DIRS})
    endif()
    if(CHAMELEON_USE_MAGMA)
        list(APPEND libs_for_examples
             ${MAGMA_LIBRARIES}
        )
        link_directories(${MAGMA_LIBRARY_DIRS})
    endif()

    list(APPEND libs_for_examples
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

    list(APPEND libs_for_examples
         coreblas
         simulapacke
         simucblas
         ${HWLOC_LIBRARIES}
         ${EXTRA_LIBRARIES}
    )

endif()

link_directories(${HWLOC_LIBRARY_DIRS})

# message(STATUS "libs timings: ${libs_for_examples}")
foreach(_example ${EXAMPLES})
    get_filename_component(_name_exe ${_example} NAME_WE)
    add_executable(${_name_exe} ${_example})
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_examples})
    install(TARGETS ${_name_exe}
            DESTINATION lib/chameleon/example/basic_zposv)
endforeach()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
