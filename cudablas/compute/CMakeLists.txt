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
#  @author Florent Pruvost
#  @date 16-09-2015
#
###

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
set(CUDABLAS_SRCS_GENERATED "")
set(ZSRC
    cuda_zgemerge.c
    cuda_zgemm.c
    cuda_zhemm.c
    cuda_zher2k.c
    cuda_zherfb.c
    cuda_zherk.c
    cuda_zlarfb.c
    cuda_zparfb.c
    cuda_zsymm.c
    cuda_zsyr2k.c
    cuda_zsyrk.c
    cuda_ztpmqrt.c
    cuda_ztrmm.c
    cuda_ztrsm.c
    cuda_ztsmlq.c
    cuda_ztsmqr.c
    cuda_zttmqr.c
    cuda_zunmlqt.c
    cuda_zunmqrt.c
    )

if( CHAMELEON_USE_CUBLAS_V2 )
  set(ZSRC
    ${ZSRC}
    cuda_zgeadd.c
    )
endif( CHAMELEON_USE_CUBLAS_V2 )

# Former MAGMA files that are no longer supported
# if( CHAMELEON_USE_MAGMA )
#   set(ZSRC
#     ${ZSRC}
#     cuda_zgelqt.c
#     cuda_zgeqrt.c
#     cuda_zgessm.c
#     cuda_zgetrf.c
#     cuda_zlauum.c
#     cuda_zpotrf.c
#     cuda_zssssm.c
#     cuda_ztrtri.c
#     cuda_ztslqt.c
#     cuda_ztsqrt.c
#     cuda_ztstrf.c
#     )
# endif()

precisions_rules_py(
  CUDABLAS_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}")

set(CUDABLAS_SRCS
  ${CUDABLAS_SRCS_GENERATED}
  )

if (CHAMELEON_USE_CUBLAS_V2)
  set(CUDABLAS_SRCS
    ${CUDABLAS_SRCS}
    cudaglobal.c
    )
endif (CHAMELEON_USE_CUBLAS_V2)

# Compile step
# ------------
add_library(cudablas ${CUDABLAS_SRCS})
add_dependencies(cudablas cudablas_include)
set_property(TARGET cudablas PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries(cudablas coreblas ${CUDA_LIBRARIES})

# installation
# ------------
install(TARGETS cudablas
        DESTINATION lib)

###
### END CMakeLists.txt
###
