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
#  @author Florent Pruvost
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

set(FIGURES 
    tile_lu.pdf
    tile_lu.jpg
    tile_layout.pdf
    tile_layout.jpg
    trace_qr.pdf
    trace_qr.jpg
    potri_async.png
    morse_header.png
    )
set(FIGURES_HTML 
    tile_lu.jpg
    tile_layout.jpg
    trace_qr.jpg
    potri_async.png
    morse_header.png
    )

foreach(_fig ${FIGURES})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figures/${_fig}
                   ${CMAKE_CURRENT_BINARY_DIR}/${_fig}
                   COPYONLY)
endforeach()

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
            DESTINATION share/chameleon/info)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/users_guide.html
            DESTINATION share/chameleon/html)
    foreach(_fig ${FIGURES_HTML})
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${_fig}
                DESTINATION share/chameleon/html)
    endforeach()

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
            DESTINATION share/chameleon/pdf)

else()
    message(STATUS "Looking for texi2dvi - not found")
endif()

###
### END CMakeLists.txt
###
