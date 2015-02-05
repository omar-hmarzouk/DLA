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
#  @date 02-10-2014
#
###

# Define the list of sources
# --------------------------
set(SIMULAPACKE_SRCS
    simulapacke.c
    )

# Define the list of headers
# --------------------------
set(SIMULAPACKE_HDRS
    lapacke.h
    lapacke_config.h
    lapacke_mangling.h
    )

# include directories
# -------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Force generation of headers
# ---------------------------
add_custom_target(simulapacke_include ALL SOURCES ${SIMULAPACKE_HDRS})

# Compile step
# ------------
add_library(simulapacke ${SIMULAPACKE_SRCS})
add_dependencies(simulapacke simulapacke_include)
set_property(TARGET simulapacke PROPERTY LINKER_LANGUAGE Fortran)

# installation
# ------------
install(FILES ${SIMULAPACKE_HDRS}
        DESTINATION include/chameleon/simucore/simulapacke)
install(TARGETS simulapacke
        DESTINATION lib)

###
### END CMakeLists.txt
###
