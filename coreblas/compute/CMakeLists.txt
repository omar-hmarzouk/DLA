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

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
set(COREBLAS_SRCS_GENERATED "")
set(ZSRC
    core_dzasum.c
    core_zgeadd.c
    core_zgelqt.c
    core_zgemm.c
    core_zgeqrt.c
    core_zgesplit.c
    core_zgessm.c
    core_zgessq.c
    core_zgetf2_nopiv.c
    core_zgetrf.c
    core_zgetrf_incpiv.c
    core_zgetrf_nopiv.c
    core_zhemm.c
    core_zher2k.c
    core_zherk.c
    core_zhessq.c    
    core_zlacpy.c
    core_zlag2c.c
    core_zlange.c
    core_zlanhe.c
    core_zlansy.c
    core_zlantr.c
    core_zlaset2.c
    core_zlaset.c
    core_zlauum.c
    core_zpamm.c
    core_zparfb.c
    core_zpemv.c
    core_zplghe.c
    core_zplgsy.c
    core_zplrnt.c
    core_zpotrf.c
    core_zssssm.c
    core_zsymm.c
    core_zsyr2k.c
    core_zsyrk.c
    core_zsyssq.c
    core_zsytf2_nopiv.c
    core_ztrasm.c
    core_ztrmm.c
    core_ztrsm.c
    core_ztrssq.c
    core_ztrtri.c
    core_ztslqt.c
    core_ztsmlq.c
    core_ztsmqr.c
    core_ztsqrt.c
    core_ztstrf.c
    core_zttlqt.c
    core_zttmlq.c
    core_zttmqr.c
    core_zttqrt.c
    core_zunmlq.c
    core_zunmqr.c
    )

precisions_rules_py(COREBLAS_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}")

set(COREBLAS_SRCS
    global.c
    ${COREBLAS_SRCS_GENERATED}
    )

# Compile step
# ------------
add_library(coreblas ${COREBLAS_SRCS})
add_dependencies(coreblas coreblas_include)
set_property(TARGET coreblas PROPERTY LINKER_LANGUAGE Fortran)

# # Dependencies step
# # -----------------
# foreach(_external_package "blas" "cblas" "lapack" "lapacke")
#     string(TOUPPER "${_external_package}" _NAMEVAR)
#     if(${_NAMEVAR}_EP)
#         add_dependencies(coreblas ${_external_package}_build)
#     endif()
# endforeach()

# installation
# ------------
install(TARGETS coreblas
        DESTINATION lib)

###
### END CMakeLists.txt
###
