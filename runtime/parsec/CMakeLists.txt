###
#
# @copyright (c) 2009-2015 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2015 Inria. All rights reserved.
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

# check if magma_dgetrf_incpiv_gpu is accessible in libmagma and activate it in chameleon
if ( CBLAS_FOUND AND LAPACKE_FOUND AND LAPACK_FOUND AND CUDA_FOUND AND CUDA_CUBLAS_LIBRARIES AND MAGMA_FOUND )
    set(CMAKE_REQUIRED_LIBRARIES "${CBLAS_LIBRARIES};${LAPACKE_LIBRARIES};${LAPACK_SEQ_LIBRARIES};${CUDA_LIBRARIES};${CUDA_CUBLAS_LIBRARIES};${MAGMA_LIBRARIES};${COREBLAS_LIBRARIES}")
    unset(MAGMA_DGETRF_INCPIV_GPU_FOUND CACHE)
    check_function_exists(magma_dgetrf_incpiv_gpu MAGMA_DGETRF_INCPIV_GPU_FOUND)
    if ( MAGMA_DGETRF_INCPIV_GPU_FOUND )
        message(STATUS "Set HAVE_MAGMA_GETRF_INCPIV_GPU")
        set(HAVE_MAGMA_GETRF_INCPIV_GPU 1)
    endif()
endif()

# Generate headers for all possible precisions
# --------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
#set(ZHDR
#)
#
#precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
#                    PRECISIONS "s;d;c;z;ds;zc"
#                    TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
    include/morse_parsec.h
    ${RUNTIME_HDRS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(runtime_parsec_include ALL SOURCES ${RUNTIME_HDRS})

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
add_library(chameleon_parsec ${RUNTIME_SRCS})
set_property(TARGET chameleon_parsec PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_parsec PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_parsec coreblas ${PARSEC_LIBRARIES_DEP})
if(CHAMELEON_USE_CUDA)
  target_link_libraries(chameleon_parsec cudablas)
endif(CHAMELEON_USE_CUDA)

add_dependencies(chameleon_parsec
  chameleon_include
  coreblas_include
  control_include
  runtime_parsec_include
)

if (CHAMELEON_USE_CUDA AND NOT CHAMELEON_SIMULATION)
    add_dependencies(chameleon_starpu cudablas_include)
endif()

# installation
# ------------
install(TARGETS chameleon_parsec
        DESTINATION lib)

###
### END CMakeLists.txt
###
