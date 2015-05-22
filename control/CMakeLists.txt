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


# Generate the morse headers for all possible precisions
# ------------------------------------------------------
set(CHAMELEON_HDRS_GENERATED "")
set(ZHDR
    compute_z.h
)

precisions_rules_py(CHAMELEON_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "s;d;c;z;ds;zc")

# Define the list of headers
# --------------------------
set(CHAMELEON_HDRS
    async.h
    auxiliary.h
    bulge.h
    common.h
    context.h
    descriptor.h
    gkkleader.h
    global.h
    morsewinthread.h
    tile.h
    workspace.h
    )

set(HDR_INSTALL ${CHAMELEON_HDRS})
foreach( hdr_file ${CHAMELEON_HDRS_GENERATED} )
    list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

list(APPEND CHAMELEON_HDRS ${CHAMELEON_HDRS_GENERATED})

# Force generation of headers
# ---------------------------
add_custom_target(control_include ALL SOURCES ${CHAMELEON_HDRS})


# installation
# ------------
install(FILES ${HDR_INSTALL}
        DESTINATION include/chameleon/control)

###
### END CMakeLists.txt
###
