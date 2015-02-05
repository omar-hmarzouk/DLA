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

# Generate header files
# ---------------------
set(COREBLAS_HDRS_GENERATED "")
set(ZHDR
    coreblas_z.h
    coreblas_zc.h
)
precisions_rules_py(COREBLAS_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "s;d;c;z;zc;ds" )

# Define the list of headers
# --------------------------
set(COREBLAS_HDRS
    cblas.h
    coreblas.h
    lapacke.h
    lapacke_config.h
    lapacke_mangling.h
    ${COREBLAS_HDRS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(coreblas_include ALL SOURCES ${COREBLAS_HDRS})

set(HDR_INSTALL "cblas.h;coreblas.h;lapacke.h;lapacke_config.h;lapacke_mangling.h")
foreach( hdr_file ${COREBLAS_HDRS_GENERATED} )
    list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# installation
# ------------
install(FILES ${HDR_INSTALL}
        DESTINATION include/chameleon/coreblas/include)

###
### END CMakeLists.txt
###
