###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University 
#                          of Tennessee Research Foundation. 
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 IPB. All rights reserved. 
#
###
#
#  @file PrintOpts.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver. 
# 
#  @version 1.1.0
#  @author Florent Pruvost
#  @date 10-11-2014
#   
###
message("\nConfiguration of package `magmamorse':")
message("        BUILDNAME ...........: ${BUILDNAME}")
message("        SITE ................: ${SITE}")
message(" ")
message("        Compiler: C .........: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("                version .....: ${COMPILER_C_VERSION}")
if(CMAKE_CXX_COMPILER)
  message("        Compiler: C++ .......: ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
  message("                version .....: ${COMPILER_CXX_VERSION}")
endif()
if(CMAKE_Fortran_COMPILER)
  message("        Compiler: Fortran ...: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")
  message("                version .....: ${COMPILER_Fortran_VERSION}")
endif()
if(MAGMAMORSE_USE_MPI)
  message("        Compiler: MPI .......: ${MPI_C_COMPILER}")
  message("        compiler flags ......: ${MPI_C_COMPILE_FLAGS}")
endif()
message("        Linker: .............: ${CMAKE_LINKER}")
message(" ")
message("        Build type ..........: ${CMAKE_BUILD_TYPE}")
message("        Build shared ........: ${BUILD_SHARED_LIBS}")
message("        CFlags ..............: ${CMAKE_C_FLAGS}")
message("        CXXFlags ............: ${CMAKE_CXX_FLAGS}")
message("        LDFlags .............: ${CMAKE_C_LINK_FLAGS}")
message(" ")
message("        Implementation paradigm")
message("        CUDA ................: ${MAGMAMORSE_USE_CUDA}")
message("        MPI .................: ${MAGMAMORSE_USE_MPI}")
message(" ")
message("        Runtime specific")
message("        QUARK ...............: ${MAGMAMORSE_SCHED_QUARK}")
message("        StarPU ..............: ${MAGMAMORSE_SCHED_STARPU}")
message("        FxT .................: ${MAGMAMORSE_USE_FXT}")
message(" ")
message("        Kernels specific")
message("        BLAS ................: ${BLA_VENDOR}")
message("        MAGMA ...............: ${MAGMAMORSE_USE_MAGMA}")
message(" ")
message("        Simulation mode .....: ${MAGMAMORSE_SIMULATION}")
message(" ")
message("        Binaries to build")
message("        documentation ........: ${MAGMAMORSE_ENABLE_DOCS}")
message("        example ..............: ${MAGMAMORSE_ENABLE_EXAMPLE}")
message("        testing ..............: ${MAGMAMORSE_ENABLE_TESTING}")
message("        timing ...............: ${MAGMAMORSE_ENABLE_TIMING}")
message(" ")
message("        Magmamorse dependencies :")
foreach (_dep ${MAGMAMORSE_DEP})
    message("                                  ${_dep}")
endforeach ()
message(" ")
message("        INSTALL_PREFIX ......: ${CMAKE_INSTALL_PREFIX}")
