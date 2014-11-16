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


# Generate the morse headers for all possible precisions
# ------------------------------------------------------
set(MAGMAMORSE_HDRS_GENERATED "")
set(ZHDR
    compute_z.h
)

precisions_rules_py(MAGMAMORSE_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "s;d;c;z;ds;zc")

# Define the list of headers
# --------------------------
set(MAGMAMORSE_HDRS
    async.h
    auxiliary.h
    common.h
    context.h
    descriptor.h
    morsewinthread.h
    workspace.h
    ${MAGMAMORSE_HDRS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(control_include ALL SOURCES ${MAGMAMORSE_HDRS})

set(HDR_INSTALL
    async.h
    auxiliary.h
    common.h
    context.h
    descriptor.h
    morsewinthread.h
    workspace.h
    )
foreach( hdr_file ${MAGMAMORSE_HDRS_GENERATED} )
    list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# installation
# ------------
install(FILES ${HDR_INSTALL}
        DESTINATION include)

###
### END CMakeLists.txt
###
