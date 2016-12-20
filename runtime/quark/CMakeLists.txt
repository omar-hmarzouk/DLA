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
    ${CODELETS_ZSRC}
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
set_property(TARGET chameleon_quark PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_quark coreblas ${QUARK_LIBRARIES_DEP})
if(CHAMELEON_USE_CUDA OR CHAMELEON_USE_MAGMA)
  target_link_libraries(chameleon_quark cudablas)
endif(CHAMELEON_USE_CUDA OR CHAMELEON_USE_MAGMA)

add_dependencies(chameleon_quark
  chameleon_include
  coreblas_include
  control_include
  runtime_quark_include
)
if (NOT CHAMELEON_SIMULATION)
    add_dependencies(chameleon_quark coreblas_include)
endif()

# installation
# ------------
install(TARGETS chameleon_quark
        DESTINATION lib)

###
### END CMakeLists.txt
###
