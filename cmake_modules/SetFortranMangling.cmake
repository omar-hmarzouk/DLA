###
#
# @copyright (c) 2017 Inria. All rights reserved.
#
###
#
#  @file SetFortranMangling.cmake
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
#  @author Florent Pruvost
#  @date 25-04-2017
#
# Detect Fortran mangling and define the proper API symbols
###

include(FortranCInterface)
## Ensure that the fortran compiler and c compiler specified are compatible
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/morse_mangling.h
                         MACRO_NAMESPACE "MORSE_"
                         SYMBOLS
                         MORSE_INIT
                         MORSE_FINALIZE
                         MORSE_ENABLE
                         MORSE_DISABLE
                         MORSE_SET
                         MORSE_GET
                         MORSE_DEALLOC_HANDLE
                         MORSE_VERSION
                         MORSE_DESC_CREATE
                         MORSE_DESC_DESTROY
                         MORSE_LAPACK_TO_TILE
                         MORSE_TILE_TO_LAPACK
                         )
