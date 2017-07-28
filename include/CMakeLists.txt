###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2017 Inria. All rights reserved.
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


# Generate the morse headers for all possible precisions
# ------------------------------------------------------
set(CHAMELEON_HDRS_GENERATED "")
set(ZHDR
    chameleon/morse_z.h
    chameleon/morse_zc.h
    chameleon/runtime_z.h
    chameleon/runtime_zc.h
)

precisions_rules_py(CHAMELEON_HDRS_GENERATED "${ZHDR}"
  TARGETDIR  "chameleon"
  PRECISIONS "z;c;d;s;zc;ds" )

configure_file("chameleon/chameleon_config.h.in"
               "chameleon/chameleon_config.h"
               @ONLY)

# Define the list of headers
# --------------------------
set(CHAMELEON_HDRS
    chameleon/chameleon_config.h
    chameleon/morse_constants.h
    chameleon/morse_kernels.h
    chameleon/morse_simulate.h
    chameleon/morse_struct.h
    chameleon/morse_types.h
    chameleon/morse_mangling.h
    chameleon/runtime.h
    ${CHAMELEON_HDRS_GENERATED}
    )

# set(flags_to_add "")
# foreach(_prec ${CHAMELEON_PRECISION})
#     set(flags_to_add "${flags_to_add} -DPRECISION_${_prec}")
# endforeach()
# set_source_files_properties(runtime.h PROPERTIES COMPILE_FLAGS "${flags_to_add}")

# Force generation of headers
# ---------------------------
add_custom_target(chameleon_include
  ALL SOURCES morse.h ${CHAMELEON_HDRS} )

# set(HDR_INSTALL
#     chameleon/chameleon_config.h
#     chameleon/morse_constants.h
#     chameleon/morse_kernels.h
#     chameleon/morse_simulate.h
#     chameleon/morse_struct.h
#     chameleon/morse_types.h
#     chameleon/morse_mangling.h
#     chameleon/runtime.h
#     ${CHAMELEON_HDRS_GENERATED}
#     )
# foreach( hdr_file ${CHAMELEON_HDRS_GENERATED} )
#     list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
# endforeach()

# installation
# ------------
install( FILES morse.h
         DESTINATION include )

install( FILES ${CHAMELEON_HDRS}
         DESTINATION include/chameleon )

###
### END CMakeLists.txt
###
