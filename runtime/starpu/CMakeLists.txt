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


# check if magma_dgetrf_incpiv_gpu is accessible in libmagma and activate it in chameleon
if ( CBLAS_FOUND AND LAPACKE_FOUND AND LAPACK_FOUND AND CUDA_FOUND AND CUDA_CUBLAS_LIBRARIES AND MAGMA_FOUND )
    set(CMAKE_REQUIRED_LIBRARIES "${CBLAS_LIBRARIES};${LAPACKE_LIBRARIES};${LAPACK_SEQ_LIBRARIES};${CUDA_LIBRARIES};${CUDA_CUBLAS_LIBRARIES};${MAGMA_LIBRARIES};${COREBLAS_LIBRARIES}")
    unset(MAGMA_DGETRF_INCPIV_GPU_FOUND CACHE)
    check_function_exists(magma_dgetrf_incpiv_gpu MAGMA_DGETRF_INCPIV_GPU_FOUND)
    if ( MAGMA_DGETRF_INCPIV_GPU_FOUND )
        message(STATUS "Set HAVE_MAGMA_GETRF_INCPIV_GPU")
        add_definitions(-DHAVE_MAGMA_GETRF_INCPIV_GPU)
    endif()
endif()

# Generate headers for all possible precisions
# --------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
set(ZHDR
  include/runtime_codelet_z.h
)

precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "s;d;c;z;ds;zc"
                    TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
    include/morse_starpu.h
    include/runtime_codelet_profile.h
    include/runtime_codelets.h
    include/runtime_profiling.h
    include/runtime_workspace.h
    ${RUNTIME_HDRS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(runtime_starpu_include ALL SOURCES ${RUNTIME_HDRS})

# installation
# ------------
#install(FILES ${RUNTIME_HDRS}
#        DESTINATION include)

# Generate the morse common for all possible precisions
# -----------------------------------------------------
set(RUNTIME_COMMON_GENERATED "")
set(ZSRC
    control/runtime_zprofiling.c
    control/runtime_zlocality.c
    )

precisions_rules_py(RUNTIME_COMMON_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}"
                    TARGETDIR "control")

set(RUNTIME_COMMON
    control/runtime_control.c
    control/runtime_context.c
    control/runtime_descriptor.c
    control/runtime_async.c
    control/runtime_options.c
    control/runtime_profiling.c
    control/runtime_workspace.c
    codelets/codelet_dataflush.c
    ${RUNTIME_COMMON_GENERATED}
    )

set(flags_to_add "")
foreach(_prec ${CHAMELEON_PRECISION})
    set(flags_to_add "${flags_to_add} -DPRECISION_${_prec}")
endforeach()
set_source_files_properties(control/runtime_profiling.c PROPERTIES COMPILE_FLAGS "${flags_to_add}")

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
    codelets/codelet_zcallback.c
    codelets/codelet_ztile_zero.c
    codelets/codelet_zasum.c
    ##################
    # BLAS 1
    ##################
    codelets/codelet_zaxpy.c
    ##################
    # BLAS 3
    ##################
    codelets/codelet_zgemm.c
    codelets/codelet_zhemm.c
    codelets/codelet_zher2k.c
    codelets/codelet_zherk.c
    codelets/codelet_zsymm.c
    codelets/codelet_zsyr2k.c
    codelets/codelet_zsyrk.c
    codelets/codelet_ztrmm.c
    codelets/codelet_ztrsm.c
    ##################
    # LAPACK
    ##################
    codelets/codelet_zgeadd.c
    codelets/codelet_zgelqt.c
    codelets/codelet_zgeqrt.c
    codelets/codelet_zgessm.c
    codelets/codelet_zgessq.c
    codelets/codelet_zgetrf.c
    codelets/codelet_zgetrf_incpiv.c
    codelets/codelet_zgetrf_nopiv.c
    codelets/codelet_zhessq.c
    codelets/codelet_zlacpy.c
    codelets/codelet_zlange.c
    codelets/codelet_zlanhe.c
    codelets/codelet_zlansy.c
    codelets/codelet_zlantr.c
    codelets/codelet_zlaset2.c
    codelets/codelet_zlaset.c
    codelets/codelet_zlauum.c
    codelets/codelet_zplghe.c
    codelets/codelet_zplgsy.c
    codelets/codelet_zplrnt.c
    codelets/codelet_zplssq.c
    codelets/codelet_zpotrf.c
    codelets/codelet_zssssm.c
    codelets/codelet_zsyssq.c
    codelets/codelet_zsytrf_nopiv.c
    codelets/codelet_ztrasm.c
    codelets/codelet_ztrssq.c
    codelets/codelet_ztrtri.c
    codelets/codelet_ztslqt.c
    codelets/codelet_ztsmlq.c
    codelets/codelet_ztsmqr.c
    codelets/codelet_ztsqrt.c
    codelets/codelet_ztstrf.c
    codelets/codelet_zttlqt.c
    codelets/codelet_zttmlq.c
    codelets/codelet_zttmqr.c
    codelets/codelet_zttqrt.c
    codelets/codelet_zunmlq.c
    codelets/codelet_zunmqr.c
    )
list(REMOVE_DUPLICATES ZSRC)

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}"
                    TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  )

# Add library
# -----------
add_library(chameleon_starpu ${RUNTIME_SRCS})
set_property(TARGET chameleon_starpu PROPERTY LINKER_LANGUAGE Fortran)

add_dependencies(chameleon_starpu
  chameleon_include
  coreblas_include
  control_include
  runtime_starpu_include
)

# installation
# ------------
install(TARGETS chameleon_starpu
        DESTINATION lib)

###
### END CMakeLists.txt
###
