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
#  @version 0.9.0
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
    set(CHAMELEON_PKGCONFIG_LIBS "")
    # The link flags for private libraries required by this package but not
    # exposed to applications
    set(CHAMELEON_PKGCONFIG_LIBS_PRIVATE "")
    # A list of packages required by this package
    set(CHAMELEON_PKGCONFIG_REQUIRED "")
    # A list of private packages required by this package but not exposed to
    # applications
    set(CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE "")

    list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon)
    if(CHAMELEON_SCHED_STARPU)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_starpu)
        if ( CHAMELEON_USE_MPI )
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED
            starpumpi-${STARPU_VERSION_STRING})
        else()
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED
            starpu-${STARPU_VERSION_STRING})
        endif()
    elseif(CHAMELEON_SCHED_QUARK)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS -lchameleon_quark)
        list(APPEND CHAMELEON_PKGCONFIG_LIBS "-l${QUARK_quark_LIBRARY}")
    endif()


    if(NOT CHAMELEON_SIMULATION)

        if(CHAMELEON_USE_CUDA)
            list(APPEND CHAMELEON_PKGCONFIG_LIBS ${CUDA_LIBRARIES})
        endif()

        if(CHAMELEON_USE_MAGMA)
            list(APPEND CHAMELEON_PKGCONFIG_REQUIRED magma)
        endif()

        list(APPEND CHAMELEON_PKGCONFIG_LIBS
        -lcoreblas
        ${LAPACKE_LIBRARIES}
        ${CBLAS_LIBRARIES}
        ${EXTRA_LIBRARIES}
        )

        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED hwloc)

    else(NOT CHAMELEON_SIMULATION)

        list(APPEND CHAMELEON_PKGCONFIG_LIBS 
        -lcoreblas
        -lsimulapacke
        -lsimucblas
        ${EXTRA_LIBRARIES}
        )

        list(APPEND CHAMELEON_PKGCONFIG_REQUIRED hwloc)

    endif(NOT CHAMELEON_SIMULATION)

    # Define required package
    # -----------------------
    set(CHAMELEON_PKGCONFIG_LIBS_CPY "${CHAMELEON_PKGCONFIG_LIBS}")
    set(CHAMELEON_PKGCONFIG_LIBS "")
    foreach(_dep ${CHAMELEON_PKGCONFIG_LIBS_CPY})
        get_filename_component(dep_we ${_dep} NAME)
        STRING(REPLACE "lib"    "-l" dep_we "${dep_we}")
        STRING(REPLACE ".so"    ""   dep_we "${dep_we}")
        STRING(REPLACE ".a"     ""   dep_we "${dep_we}")
        STRING(REPLACE ".dylib" ""   dep_we "${dep_we}")
        STRING(REPLACE ".dll"   ""   dep_we "${dep_we}")
        list(APPEND CHAMELEON_PKGCONFIG_LIBS ${dep_we})
    endforeach()

    list(REMOVE_DUPLICATES CHAMELEON_PKGCONFIG_LIBS)
    list(REMOVE_DUPLICATES CHAMELEON_PKGCONFIG_LIBS_PRIVATE)
    list(REMOVE_DUPLICATES CHAMELEON_PKGCONFIG_REQUIRED)
    list(REMOVE_DUPLICATES CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE)

    STRING(REPLACE ";" " " CHAMELEON_PKGCONFIG_LIBS "${CHAMELEON_PKGCONFIG_LIBS}")
    STRING(REPLACE ";" " " CHAMELEON_PKGCONFIG_LIBS_PRIVATE "${CHAMELEON_PKGCONFIG_LIBS_PRIVATE}")
    STRING(REPLACE ";" " " CHAMELEON_PKGCONFIG_REQUIRED "${CHAMELEON_PKGCONFIG_REQUIRED}")
    STRING(REPLACE ";" " " CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE "${CHAMELEON_PKGCONFIG_REQUIRED_PRIVATE}")

    # Create .pc file
    # ---------------
    SET(_output_file "${CMAKE_BINARY_DIR}/chameleon.pc")
    # TODO: add url of MORSE releases in .pc file
    CONFIGURE_FILE("${_file}" "${_output_file}" @ONLY)

    # installation
    # ------------
    INSTALL(FILES ${_output_file} DESTINATION lib/pkgconfig)

ENDMACRO(GENERATE_PKGCONFIG_FILE)

##
## @end file GenPkgConfig.cmake
##
