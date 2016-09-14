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

# Define the list of sources
# --------------------------
set(SIMUCBLAS_SRCS
    simucblas.c
    )

# Compile step
# ------------
add_library(simucblas ${SIMUCBLAS_SRCS})
# this library is theoreticaly independent from coreblas but in order to avoid 
# the duplication of cblas headers in here we rely on the ones lying in 
# coreblas/include 
add_dependencies(simucblas coreblas_include)
set_property(TARGET simucblas PROPERTY LINKER_LANGUAGE Fortran)

# installation
# ------------
install(TARGETS simucblas
        DESTINATION lib)

###
### END CMakeLists.txt
###
