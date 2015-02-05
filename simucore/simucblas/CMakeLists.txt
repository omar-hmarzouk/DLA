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

# Define the list of headers
# --------------------------
set(SIMUCBLAS_HDRS
    cblas.h
    )

# include directories
# -------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Force generation of headers
# ---------------------------
add_custom_target(simucblas_include ALL SOURCES ${SIMUCBLAS_HDRS})

# Compile step
# ------------
add_library(simucblas ${SIMUCBLAS_SRCS})
add_dependencies(simucblas simucblas_include)
set_property(TARGET simucblas PROPERTY LINKER_LANGUAGE Fortran)

# installation
# ------------
install(FILES ${SIMUCBLAS_HDRS}
        DESTINATION include/chameleon/simucore/simucblas)
install(TARGETS simucblas
        DESTINATION lib)

###
### END CMakeLists.txt
###
