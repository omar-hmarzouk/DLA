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
#  @file CMakeLists.txt
#
#  @project MORSE
#  MORSE is a software package provided by:
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
#  @date 13-07-2012
#
###

# Generate the morse auxiliary headers for all possible precisions
# ----------------------------------------------------------------
set(TIMING_AUX_HDRS_GENERATED "")
set(ZHDR
    timing_zauxiliary.h
   )

precisions_rules_py(TIMING_AUX_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "s;d;c;z;ds;zc" )

set(TIMING_AUX_HDRS
    timing.h
    timing.c
    ${TIMING_AUX_HDRS_GENERATED}
   )

# Force generation of headers
# ---------------------------
add_custom_target(timing_include ALL SOURCES ${TIMING_AUX_HDRS})

if (NOT CHAMELEON_SIMULATION)
    # Generate the morse auxiliary sources for all possible precisions
    # ----------------------------------------------------------------
    set(TIMING_AUX_SRCS_GENERATED "")
    set(ZSRC
        timing_zauxiliary.c
       )

    precisions_rules_py(TIMING_AUX_SRCS_GENERATED "${ZSRC}"
                        PRECISIONS "${CHAMELEON_PRECISION}")

    set(TIMING_SRCS
        ${TIMING_AUX_SRCS_GENERATED}
       )

    # Create libchameleon_timing.a
    # -----------------------------
    add_library(chameleon_timing STATIC ${TIMING_SRCS})
    set_property(TARGET chameleon_timing PROPERTY LINKER_LANGUAGE Fortran)
    add_dependencies(chameleon_timing timing_include)
endif()

# Generate the morse testing sources for all possible precisions
# --------------------------------------------------------------
set(TIMINGS "")
if (NOT CHAMELEON_SIMULATION)
    set(ZSRC
        ##################
        # BLAS 3
        ##################
        time_zgemm.c
        time_zgemm_tile.c
        time_ztrsm.c
        ##################
        # LAPACK
        ##################
        time_zgels.c
        time_zgels_tile.c
        time_zgeqrf.c
        time_zgeqrf_hqr.c
        time_zgeqrf_hqr_tile.c
        time_zgeqrf_tile.c
        time_zgelqf.c
        time_zgelqf_tile.c
        time_zgeqrs_tile.c
        time_zgetrf_incpiv.c
        time_zgetrf_incpiv_tile.c
        time_zgetrf_nopiv.c
        time_zgetrf_nopiv_tile.c
        time_zgetrs_incpiv_tile.c
        time_zgetrs_nopiv_tile.c
        #time_zgetrf.c
        #time_zgetrf_tile.c
        time_zposv.c
        time_zposv_tile.c
        time_zpotrf.c
        time_zpotrf_tile.c
        time_zpotrs_tile.c
        time_zsytrf_tile.c
        time_zpotri_tile.c
        ##################
        # MIXED PRECISION
        ##################
        #time_zcgesv.c
        #time_zcgesv_tile.c
        #time_zcposv.c
        #time_zcposv_tile.c
        ##################
        # OTHERS
        ##################
        time_zlange.c
        time_zlange_tile.c
        #time_zgebrd_tile.c
        #time_zgecfi.c
        time_zgesvd_tile.c
        #time_zgetrf_reclap.c
        #time_zgetrf_rectil.c
        time_zheevd_tile.c
        #time_zheev_tile.c
        #time_zhegv_tile.c
        #time_zlapack2tile.c
        #time_zgetri_tile.c
        #time_zgesv.c
        #time_zgesv_tile.c
        time_zgesv_incpiv.c
        time_zgesv_incpiv_tile.c
        time_zgesv_nopiv.c
        time_zgesv_nopiv_tile.c
       )
else(NOT CHAMELEON_SIMULATION)
    set(ZSRC
        ##################
        # BLAS 3
        ##################
        time_zgemm_tile.c
        ##################
        # LAPACK
        ##################
        time_zgels_tile.c
        time_zgeqrf_tile.c
        time_zgelqf_tile.c
        time_zgeqrs_tile.c
        time_zgetrf_incpiv_tile.c
        time_zgetrf_nopiv_tile.c
        time_zgetrs_incpiv_tile.c
        time_zgetrs_nopiv_tile.c
        #time_zgetrf_tile.c
        time_zposv_tile.c
        time_zpotrf_tile.c
        time_zpotrs_tile.c
        time_zsytrf_tile.c
        time_zpotri_tile.c
        ##################
        # MIXED PRECISION
        ##################
        #time_zcgesv_tile.c
        #time_zcposv_tile.c
        ##################
        # OTHERS
        ##################
        time_zlange_tile.c
        #time_zgebrd_tile.c
        #time_zgecfi.c
        #time_zgesvd_tile.c
        #time_zgetrf_reclap.c
        #time_zgetrf_rectil.c
        #time_zheevd_tile.c
        #time_zheev_tile.c
        #time_zhegv_tile.c
        #time_zlapack2tile.c
        #time_zgetri_tile.c
        #time_zgesv.c
        #time_zgesv_tile.c
        time_zgesv_incpiv_tile.c
        time_zgesv_nopiv_tile.c
       )
endif(NOT CHAMELEON_SIMULATION)

precisions_rules_py(TIMINGS "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

# Add include and link directories
# --------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Define what libraries we have to link with
# ------------------------------------------
set(libs_for_timings)
if(NOT CHAMELEON_SIMULATION)
    list(APPEND libs_for_timings chameleon_timing)
endif()

list(APPEND libs_for_timings chameleon)

if(CHAMELEON_SCHED_STARPU)
    list(APPEND libs_for_timings
      chameleon_starpu
      ${STARPU_LIBRARIES_DEP}
    )
    link_directories(${STARPU_LIBRARY_DIRS})
elseif(CHAMELEON_SCHED_PARSEC)
    list(APPEND libs_for_timings
      chameleon_parsec
      ${PARSEC_LIBRARIES_DEP}
      )
    link_directories(${PARSEC_LIBRARY_DIRS_DEP})
elseif(CHAMELEON_SCHED_QUARK)
    list(APPEND libs_for_timings
      chameleon_quark
      ${QUARK_LIBRARIES}
      )
    link_directories(${QUARK_LIBRARY_DIRS})
endif()

if(NOT CHAMELEON_SIMULATION)

    if(CHAMELEON_USE_CUDA)
        list(APPEND libs_for_timings
        cudablas)
    endif()
    if(CHAMELEON_USE_CUDA)
        link_directories(${CUDA_LIBRARY_DIRS})
    endif()

    list(APPEND libs_for_timings
    coreblas
    ${LAPACKE_LIBRARIES}
    ${TMG_LIBRARIES}
    ${CBLAS_LIBRARIES}
    ${LAPACK_SEQ_LIBRARIES}
    ${BLAS_SEQ_LIBRARIES}
    ${LIBHQR_LIBRARIES}
    ${HWLOC_LIBRARIES}
    ${EXTRA_LIBRARIES}
    )

    link_directories(${LAPACKE_LIBRARY_DIRS})
    link_directories(${TMG_LIBRARY_DIRS})
    link_directories(${LAPACK_LIBRARY_DIRS})
    link_directories(${CBLAS_LIBRARY_DIRS})
    link_directories(${BLAS_LIBRARY_DIRS})
    link_directories(${LIBHQR_LIBRARY_DIRS})

else()

    list(APPEND libs_for_timings${EXTRA_LIBRARIES})

endif()

# message("BLAS_SEQ_LIBRARIES: ${BLAS_SEQ_LIBRARIES}")
# message("CBLAS_LIBRARIES: ${CBLAS_LIBRARIES}")
# message("LAPACK_SEQ_LIBRARIES: ${LAPACK_SEQ_LIBRARIES}")
# message("LAPACKE_LIBRARIES: ${LAPACKE_LIBRARIES}")
# message("HWLOC_LIBRARIES: ${HWLOC_LIBRARIES}")

# message("LAPACKE_LIBRARY_DIRS: ${LAPACKE_LIBRARY_DIRS}")
# message("LAPACK_LIBRARY_DIRS: ${LAPACK_LIBRARY_DIRS}")
# message("CBLAS_LIBRARY_DIRS: ${CBLAS_LIBRARY_DIRS}")
# message("BLAS_LIBRARY_DIRS: ${BLAS_LIBRARY_DIRS}")
# message("HWLOC_LIBRARY_DIRS: ${HWLOC_LIBRARY_DIRS}")

# message(STATUS "libs timings: ${libs_for_timings}")
foreach(_timing ${TIMINGS})
    get_filename_component(_name_exe ${_timing} NAME_WE)
    add_executable(${_name_exe} ${_timing})
    add_dependencies(${_name_exe} timing_include)
    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_timings})
    install(TARGETS ${_name_exe}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/timing)
endforeach()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
