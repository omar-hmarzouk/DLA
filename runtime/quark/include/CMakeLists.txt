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

# Generate the morse headers for all possible precisions
# ------------------------------------------------------
set(QUARKBLAS_HEADERS_GENERATED "")
set(ZHDR 
    quark_zblas.h
    quark_zcblas.h
)

precisions_rules_py(QUARKBLAS_HEADERS_GENERATED "${ZHDR}"
                    PRECISIONS "${CHAMELEON_PRECISION}"
                    DICTIONARY "MORSE")

# Define the list of headers
# --------------------------
set(QUARKBLAS_HEADERS
    core_blas_dag.h
    morse_quark.h
    quark_blas.h
    ${QUARKBLAS_HEADERS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(quarkblas_include ALL SOURCES ${QUARKBLAS_HEADERS})

# installation
# ------------
install(FILES ${MORSE_HEADERS} 
        DESTINATION include)

###
### END CMakeLists.txt
###
