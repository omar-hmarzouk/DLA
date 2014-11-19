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
#  @version 1.0.0
#  @date 17-10-2014
#
###


add_subdirectory(basic_zposv)

if (CHAMELEON_PREC_D)
    add_subdirectory(lapack_to_morse)
else()
    message(WARNING "CHAMELEON_PREC_D is set to OFF so that lapack_to_morse"
    "tutorial cannot be built (use only double arithmetic precision).\n"
    "Please set CHAMELEON_PREC_D to ON if you want to build executables of"
    "this tutorial.")
endif()

###
### END CMakeLists.txt
###
