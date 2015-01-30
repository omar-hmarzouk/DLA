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


include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# compilation for other sources step1 and >
set(LTM_SOURCES 
    step1.c
    step2.c
    step3.c
    step4.c
    step5.c
    step6.c
   )

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_ltm)
list(APPEND libs_for_ltm
     chameleon
)

if(CHAMELEON_SCHED_STARPU)
    if ( CHAMELEON_USE_MPI )
        list(APPEND libs_for_ltm
             chameleon_starpu
             ${STARPU_MPI_LIBRARIES}
        )
    else()
        list(APPEND libs_for_ltm
             chameleon_starpu
             ${STARPU_SHM_LIBRARIES}
        )
    endif()
    link_directories(${STARPU_LIBRARY_DIRS})
elseif(CHAMELEON_SCHED_QUARK)
    list(APPEND libs_for_ltm
         chameleon_quark
         ${QUARK_LIBRARIES}
      )
    link_directories(${QUARK_LIBRARY_DIRS})
endif()


# specific compilation for step0 because we potentially want to use 
# multithreaded BLAS and LAPACK libraries for this step
unset(libs_for_step0)

if(NOT CHAMELEON_SIMULATION)

    if(CHAMELEON_USE_CUDA)
        list(APPEND libs_for_ltm
             ${CUDA_LIBRARIES}
        )
        link_directories(${CUDA_LIBRARY_DIRS})
    endif()
    if(CHAMELEON_USE_MAGMA)
        list(APPEND libs_for_ltm
             ${MAGMA_LIBRARIES}
        )
        link_directories(${MAGMA_LIBRARY_DIRS})
    endif()

    list(APPEND libs_for_step0 ${libs_for_ltm})

    list(APPEND libs_for_ltm
         coreblas
         ${LAPACKE_LIBRARIES}
         ${CBLAS_LIBRARIES}
         ${LAPACK_SEQ_LIBRARIES}
         ${BLAS_SEQ_LIBRARIES}
         ${HWLOC_LIBRARIES}
         ${EXTRA_LIBRARIES}
    )

    if( BLA_VENDOR MATCHES "Intel10_64lp*" OR BLA_VENDOR MATCHES "ACML*")
        if(BLAS_PAR_LIBRARIES)
            set(CBLAS_LIBRARIES "${BLAS_PAR_LIBRARIES}")
        endif()
        if(LAPACK_PAR_LIBRARIES)
            set(LAPACKE_LIBRARIES "${LAPACK_PAR_LIBRARIES}")
        endif()
    endif()
    if (BLAS_PAR_LIBRARIES AND LAPACK_PAR_LIBRARIES)
        list(APPEND libs_for_step0
             coreblas
             ${LAPACKE_LIBRARIES}
             ${CBLAS_LIBRARIES} 
             ${LAPACK_PAR_LIBRARIES}
             ${BLAS_PAR_LIBRARIES}
             ${HWLOC_LIBRARIES}
             ${EXTRA_LIBRARIES}
        )

    else()
        list(APPEND libs_for_step0
             coreblas
             ${LAPACKE_LIBRARIES}
             ${CBLAS_LIBRARIES}
             ${LAPACK_SEQ_LIBRARIES}
             ${BLAS_SEQ_LIBRARIES}
             ${HWLOC_LIBRARIES}
             ${EXTRA_LIBRARIES}
        )
    endif ()

    link_directories(${LAPACKE_LIBRARY_DIRS})
    link_directories(${LAPACK_LIBRARY_DIRS})
    link_directories(${CBLAS_LIBRARY_DIRS})
    link_directories(${BLAS_LIBRARY_DIRS})

else()

    list(APPEND libs_for_ltm
         coreblas
         simulapacke
         simucblas
         ${HWLOC_LIBRARIES}
         ${EXTRA_LIBRARIES}
    )
    list(APPEND libs_for_step0
         coreblas
         simulapacke
         simucblas
         ${HWLOC_LIBRARIES}
         ${EXTRA_LIBRARIES}
    )

endif()

link_directories(${HWLOC_LIBRARY_DIRS})


# message(STATUS "libs examples: ${libs_for_ltm}")
foreach(_ltm ${LTM_SOURCES})
    get_filename_component(_name_exe ${_ltm} NAME_WE)
    add_executable(${_name_exe} ${_ltm})
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_ltm})
    install(TARGETS ${_name_exe}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/example/lapack_to_morse)
endforeach()

add_executable(step0 step0.c)
set_property(TARGET step0 PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(step0 ${libs_for_step0})
install(TARGETS step0
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/example/lapack_to_morse)

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
