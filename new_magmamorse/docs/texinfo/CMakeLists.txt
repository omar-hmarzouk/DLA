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
#  @version 1.1.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 13-07-2012
#   
###

cmake_minimum_required(VERSION 2.8)

# Create file version.texi
# ------------------------
configure_file("version.texi.in"
               "version.texi"
               @ONLY)
configure_file("users_guide.texi.in"
               "users_guide.texi"
               @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/tile_lu.pdf
               ${CMAKE_CURRENT_BINARY_DIR}/tile_lu.pdf
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/tile_lu.jpg
               ${CMAKE_CURRENT_BINARY_DIR}/tile_lu.jpg
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/tile_layout.pdf
               ${CMAKE_CURRENT_BINARY_DIR}/tile_layout.pdf
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/tile_layout.jpg
               ${CMAKE_CURRENT_BINARY_DIR}/tile_layout.jpg
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/trace_qr.pdf
               ${CMAKE_CURRENT_BINARY_DIR}/trace_qr.pdf
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/trace_qr.jpg
               ${CMAKE_CURRENT_BINARY_DIR}/trace_qr.jpg
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/potri_async.png
               ${CMAKE_CURRENT_BINARY_DIR}/potri_async.png
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/morse_header.png
               ${CMAKE_CURRENT_BINARY_DIR}/morse_header.png
               COPYONLY)

# Looking which version we can compile
# ------------------------------------
FIND_PROGRAM(MAKEINFO_COMPILER makeinfo)
FIND_PROGRAM(TEXI2DVI_COMPILER texi2dvi)
FIND_PROGRAM(TEX_COMPILER tex)
FIND_PROGRAM(DOT_COMPILER dot)

# Looking for makeinfo
# --------------------
if(MAKEINFO_COMPILER)
    # Add target
    # ----------
    add_custom_command(OUTPUT  users_guide.info
                       COMMAND ${MAKEINFO_COMPILER}
                       ARGS    users_guide.texi
                       DEPENDS users_guide.texi.in
                      )
    add_custom_command(OUTPUT  users_guide.html
                       COMMAND ${MAKEINFO_COMPILER}
                       ARGS    --html
                               --no-split
                               --css-include=${CMAKE_CURRENT_SOURCE_DIR}/morse.css
                               users_guide.texi
                       DEPENDS users_guide.texi.in
                      )
    add_custom_target(doc-info-users_guide ALL DEPENDS users_guide.info)
    add_custom_target(doc-html-users_guide ALL DEPENDS users_guide.html)

    # Installation
    # ------------
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/users_guide.info
            DESTINATION share/info)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/users_guide.html
            DESTINATION share/html)

else(MAKEINFO_COMPILER)
    message(STATUS "Looking for makeinfo - not found")

endif(MAKEINFO_COMPILER)

# Looking for texi2dvi
# --------------------
if(TEXI2DVI_COMPILER AND TEX_COMPILER)
    # Add target
    # ----------
    add_custom_command(OUTPUT  users_guide.pdf
                       COMMAND ${TEXI2DVI_COMPILER}
                       ARGS    --pdf
                               --batch
                               users_guide.texi
                       DEPENDS users_guide.texi.in
                      )
    add_custom_target(doc-pdf-users_guide ALL DEPENDS users_guide.pdf)

    # Installation
    # ------------
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/users_guide.pdf
            DESTINATION share/pdf)

else()
    message(STATUS "Looking for texi2dvi - not found")

endif()

###
### END CMakeLists.txt
###
