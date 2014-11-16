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
#  @file GenPkgConfig.cmake
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
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 10-11-2014
#   
###

###
#
# GENERATE_PKGCONFIG_FILE: generate a file .pc according to the options 
#
###
MACRO(GENERATE_PKGCONFIG_FILE _file) 
    
    # The link flags specific to this package and any required libraries 
    # that don't support PkgConfig
    set(MAGMAMORSE_PKGCONFIG_LIBS "")
    # The link flags for private libraries required by this package but not 
    # exposed to applications
    set(MAGMAMORSE_PKGCONFIG_LIBS_PRIVATE "")
    # A list of packages required by this package
    set(MAGMAMORSE_PKGCONFIG_REQUIRED "")
    # A list of private packages required by this package but not exposed to 
    # applications
    set(MAGMAMORSE_PKGCONFIG_REQUIRED_PRIVATE "")
    
    list(APPEND MAGMAMORSE_PKGCONFIG_LIBS -lmagmamorse)
    if(MAGMAMORSE_SCHED_STARPU)
        list(APPEND MAGMAMORSE_PKGCONFIG_LIBS -lmagmamorse_starpu)
        if ( MAGMAMORSE_USE_MPI )
            list(APPEND MAGMAMORSE_PKGCONFIG_REQUIRED 
            starpumpi-${MAGMAMORSE_STARPU_VERSION})
        else()
            list(APPEND MAGMAMORSE_PKGCONFIG_REQUIRED 
            starpu-${MAGMAMORSE_STARPU_VERSION})        
        endif()
    elseif(MAGMAMORSE_SCHED_QUARK)
        list(APPEND MAGMAMORSE_PKGCONFIG_LIBS -lmagmamorse_quark)
        list(APPEND MAGMAMORSE_PKGCONFIG_LIBS -lquark)
    endif()
    
    
    if(NOT MAGMAMORSE_SIMULATION)
    
        if(MAGMAMORSE_USE_CUDA)
            list(APPEND MAGMAMORSE_PKGCONFIG_LIBS ${CUDA_LIBRARIES})
        endif()
               
        if(MAGMAMORSE_USE_MAGMA)
            list(APPEND MAGMAMORSE_PKGCONFIG_REQUIRED magma)
        endif()
        
        list(APPEND MAGMAMORSE_PKGCONFIG_LIBS 
        -lcoreblas
        ${LAPACKE_LIBRARIES}
        ${CBLAS_LIBRARIES}    
        ${LAPACK_SEQ_LIBRARIES}
        ${BLAS_SEQ_LIBRARIES}
        ${EXTRA_LIBRARIES}
        )
        
        list(APPEND MAGMAMORSE_PKGCONFIG_REQUIRED hwloc)
    
    else(NOT MAGMAMORSE_SIMULATION)
        
        list(APPEND MAGMAMORSE_PKGCONFIG_LIBS 
        -lcoreblas
        -lsimulapacke    
        -lsimucblas
        ${EXTRA_LIBRARIES}        
        )
        
        list(APPEND MAGMAMORSE_PKGCONFIG_REQUIRED hwloc) 
               
    endif(NOT MAGMAMORSE_SIMULATION)
    
    list(REMOVE_DUPLICATES MAGMAMORSE_PKGCONFIG_LIBS)
    list(REMOVE_DUPLICATES MAGMAMORSE_PKGCONFIG_LIBS_PRIVATE)
    list(REMOVE_DUPLICATES MAGMAMORSE_PKGCONFIG_REQUIRED)
    list(REMOVE_DUPLICATES MAGMAMORSE_PKGCONFIG_REQUIRED_PRIVATE)

    
    # Define required package
    # -----------------------
    set(MAGMAMORSE_PKGCONFIG_LIBS_CPY "${MAGMAMORSE_PKGCONFIG_LIBS}")
    set(MAGMAMORSE_PKGCONFIG_LIBS "")
    foreach(_dep ${MAGMAMORSE_PKGCONFIG_LIBS_CPY})
        get_filename_component(dep_we ${_dep} NAME)
        STRING(REPLACE "lib"    "-l" dep_we "${dep_we}")
        STRING(REPLACE ".so"    ""   dep_we "${dep_we}")
        STRING(REPLACE ".a"     ""   dep_we "${dep_we}")
        STRING(REPLACE ".dylib" ""   dep_we "${dep_we}")
        STRING(REPLACE ".dll"   ""   dep_we "${dep_we}")
        list(APPEND MAGMAMORSE_PKGCONFIG_LIBS ${dep_we})
    endforeach()
    
    STRING(REPLACE ";" " " MAGMAMORSE_PKGCONFIG_LIBS "${MAGMAMORSE_PKGCONFIG_LIBS}")
    STRING(REPLACE ";" " " MAGMAMORSE_PKGCONFIG_LIBS_PRIVATE "${MAGMAMORSE_PKGCONFIG_LIBS_PRIVATE}")
    STRING(REPLACE ";" " " MAGMAMORSE_PKGCONFIG_REQUIRED "${MAGMAMORSE_PKGCONFIG_REQUIRED}")
    STRING(REPLACE ";" " " MAGMAMORSE_PKGCONFIG_REQUIRED_PRIVATE "${MAGMAMORSE_PKGCONFIG_REQUIRED_PRIVATE}")
      
    # Create .pc file
    # ---------------
    GET_FILENAME_COMPONENT(_output_file ${_file} NAME)
    STRING(REPLACE ".in" "" _output_file "${_output_file}")
    SET(_output_file "${CMAKE_BINARY_DIR}/${_output_file}")
    # TODO: add url of MORSE releases in .pc file
    CONFIGURE_FILE("${_file}" "${_output_file}" @ONLY)

    # installation
    # ------------
    INSTALL(FILES ${_output_file} DESTINATION lib/pkgconfig)

ENDMACRO(GENERATE_PKGCONFIG_FILE)

##
## @end file GenPkgConfig.cmake
##
