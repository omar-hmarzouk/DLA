###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2015 Inria. All rights reserved.
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
#  @author Florent Pruvost
#  @date 14-09-2015
#
###

cmake_minimum_required(VERSION 2.8)

if (NOT EZTRACE_FOUND)
    find_package(EZTRACE)
endif()

if (EZTRACE_FOUND AND EZTRACE_DIR_FOUND)

    set(EZTRACE_eztrace_create_plugin_DIR "EZTRACE_eztrace_create_plugin_DIR-NOTFOUND")
    find_path(EZTRACE_eztrace_create_plugin_DIR
      NAMES eztrace_create_plugin
      HINTS ${EZTRACE_DIR_FOUND}/bin)

    if (EZTRACE_eztrace_create_plugin_DIR)

        set(EZTRACE_CREATE_PLUGIN "${EZTRACE_eztrace_create_plugin_DIR}/eztrace_create_plugin")
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output
            COMMAND ${EZTRACE_CREATE_PLUGIN}
            ARGS ${CMAKE_CURRENT_SOURCE_DIR}/lapacke_eztrace_module
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lapacke_eztrace_module
            )
        add_custom_target(
            eztrace-module-lapacke-dir ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output
            )
        set(LAPACKE_H_DIR ${CMAKE_SOURCE_DIR}/coreblas/include)
        set(LAPACKE_H_FILE_TO_COPY 
            "${LAPACKE_H_DIR}/lapacke_config.h"
            "${LAPACKE_H_DIR}/lapacke_mangling.h"
            "${LAPACKE_H_DIR}/lapacke.h")
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output/lapacke.h
            COMMAND cp
            ARGS ${LAPACKE_H_FILE_TO_COPY} ${CMAKE_CURRENT_BINARY_DIR}/output
            DEPENDS ${LAPACKE_H_FILE_TO_COPY} ${CMAKE_CURRENT_BINARY_DIR}/output
            )
        add_custom_target(
            eztrace-module-lapacke-headers ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output ${CMAKE_CURRENT_BINARY_DIR}/output/lapacke.h
            )
        add_custom_command(
            OUTPUT libeztrace-autostart-lapacke.so libeztrace-lapacke.so libeztrace-convert-lapacke.so
            COMMAND make
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output ${CMAKE_CURRENT_BINARY_DIR}/output/lapacke.h
            )
        add_custom_target(
            eztrace-module-lapacke-libs ALL
            DEPENDS libeztrace-autostart-lapacke.so libeztrace-lapacke.so libeztrace-convert-lapacke.so
            )
        # installation
        # ------------
        install(
            FILES
            ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-autostart-lapacke.so
            ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-lapacke.so
            ${CMAKE_CURRENT_BINARY_DIR}/output/libeztrace-convert-lapacke.so
            DESTINATION ${EZTRACE_LIBRARY_DIRS}
            )

    endif (EZTRACE_eztrace_create_plugin_DIR)

endif (EZTRACE_FOUND AND EZTRACE_DIR_FOUND)

###
### END CMakeLists.txt
###
