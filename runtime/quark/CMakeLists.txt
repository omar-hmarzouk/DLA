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

cmake_minimum_required(VERSION 2.8)

# Generate the magma headers for all possible precisions
# ------------------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
set(ZHDR
  include/quark_zblas.h
  include/quark_zcblas.h
)

precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "s;d;c;z;ds;zc"
                    TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
    include/morse_quark.h
    include/quark_blas.h
    include/core_blas_dag.h
    ${RUNTIME_HDRS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(runtime_quark_include ALL SOURCES ${RUNTIME_HDRS})

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
    control/runtime_async.c
    control/runtime_context.c
    control/runtime_control.c
    control/runtime_descriptor.c
    control/runtime_options.c
    control/runtime_profiling.c
    codelets/codelet_dataflush.c
    ${RUNTIME_COMMON_GENERATED}
    )

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
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

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}"
                    TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  )

# Add library
# -----------
add_library(chameleon_quark ${RUNTIME_SRCS})
set_property(TARGET chameleon_quark PROPERTY LINKER_LANGUAGE Fortran)

add_dependencies(chameleon_quark
  chameleon_include
  coreblas_include
  control_include
  runtime_quark_include
)

# installation
# ------------
install(TARGETS chameleon_quark
        DESTINATION lib)

###
### END CMakeLists.txt
###
