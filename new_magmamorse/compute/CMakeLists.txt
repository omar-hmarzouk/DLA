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


# Define the list of sources
# --------------------------
set(MAGMAMORSE_CONTROL
    ../control/async.c
    ../control/auxiliary.c
    ../control/context.c
    ../control/control.c
    ../control/descriptor.c
    ../control/workspace.c
    ../control/tile.c
    ../control/morse_f77.c
    ../control/morse_mf77.c
#    ../control/morsewinthread.c
   )

foreach(_prec ${MAGMAMORSE_PRECISION})
    set_source_files_properties(../control/tile.c PROPERTIES COMPILE_FLAGS "-DPRECISION_${_prec}")
endforeach()

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
set(MAGMAMORSE_SRCS_GENERATED "")
set(ZSRC
    ##################
    # BLAS 3
    ##################
    pzgemm.c
    pzhemm.c
    pzher2k.c
    pzherk.c
    pzsymm.c
    pzsyr2k.c
    pzsyrk.c
    pztrmm.c
    pztrsm.c
    pztrsmpl.c
    ###
    zgemm.c
    zhemm.c
    zher2k.c
    zherk.c
    zsymm.c
    zsyr2k.c
    zsyrk.c
    ztrmm.c
    ztrsm.c
    ztrsmpl.c
    ##################
    # LAPACK
    ##################
    pzgeadd.c
    pzgelqf.c
    pzgelqfrh.c
    pzgeqrf.c
    pzgeqrfrh.c
    pzgetrf_incpiv.c
    pzgetrf_nopiv.c
    pzlacpy.c
    pzlange.c
    pzlanhe.c
    pzlansy.c
    pzlantr.c
    pzlaset2.c
    pzlaset.c
    pzlauum.c
    pzplghe.c
    pzplgsy.c
    pzplrnt.c
    pzpotrf.c
    pzsytrf.c
    pztrtri.c
    pzpotrimm.c
    pzunglq.c
    pzunglqrh.c
    pzungqr.c
    pzungqrrh.c
    pzunmlq.c
    pzunmlqrh.c
    pzunmqr.c
    pzunmqrrh.c
    ###
    zgels.c
    zgelqf.c
    zgelqs.c
    zgeqrf.c
    zgeqrs.c
    #zgesv.c
    zgesv_incpiv.c
    zgesv_nopiv.c
    #zgetrf.c
    zgetrf_incpiv.c
    zgetrf_nopiv.c
    zgetrs_incpiv.c
    zgetrs_nopiv.c
    zlacpy.c
    zlange.c
    zlanhe.c
    zlansy.c
    zlantr.c
    zlaset.c
    zlauum.c
    zplghe.c
    zplgsy.c
    zplrnt.c
    zposv.c
    zsysv.c
    zpotrf.c
    zsytrf.c
    zpotri.c
    zpotrimm.c
    zpotrs.c
    zsytrs.c
    ztrtri.c
    zunglq.c
    zungqr.c
    zunmlq.c
    zunmqr.c
    ##################
    # MIXED PRECISION
    ##################
    pzlag2c.c
    ###
    #zcgels.c
    #zcgesv.c
    #zcposv.c
    #zcungesv.c
    ##################
    # OTHERS
    ##################
    #pzgebrd_ge2tb.c
    #pzgebrd_tb2bd.c
    #pzgetmi2.c
    #pzgetrf_reclap.c
    #pzgetrf_rectil.c
    #pzhbcpy_t2bl.c
    #pzhegst.c
    #pzherbt.c
    #pzhetrd_hb2st.c
    #pzhetrd_he2hb.c
    #pzlarft_blgtrd.c
    #pzlaswp.c
    #pzlaswpc.c
    #pztrsmrv.c
    #pzunmqr_blgtrd.c
    #########################
    #zgebrd.c
    #zgecfi.c
    #zgecfi2.c
    #zgesvd.c
    #zgetmi.c
    #zgetri.c
    #zgetrs.c
    #zheev.c
    #zheevd.c
    #zhegst.c
    #zhegv.c
    #zhegvd.c
    #zhetrd.c
    #zlaswp.c
    #zlaswpc.c
    #ztrsmrv.c
    ##################
    # CONTROL
    ##################
    #pzshift.c
    #pzpack.c
    pztile.c
    ztile.c
    )

precisions_rules_py(MAGMAMORSE_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${MAGMAMORSE_PRECISION}")

set(CONTROL_SRCS_GENERATED "")
set(ZSRC
    ../control/workspace_z.c
    ../control/morse_zf77.c
)

precisions_rules_py(CONTROL_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${MAGMAMORSE_PRECISION}"
                    TARGETDIR "control" )

set(MAGMAMORSE_SRCS
    ${MAGMAMORSE_CONTROL}
    ${MAGMAMORSE_SRCS_GENERATED}
    ${CONTROL_SRCS_GENERATED}
   )

# Generate the morse fortran sources for all possible precisions
# --------------------------------------------------------------
if(HAVE_ISO_C_BINDING)
    set(MAGMAMORSE_SRCS_F_GENERATED "")
    set(ZSRCF
        ../control/morse_zcf90.F90
        ../control/morse_zf90.F90
        ../control/morse_zf90_wrappers.F90
       )
    precisions_rules_py(MAGMAMORSE_SRCS_F_GENERATED "${ZSRCF}"
                        PRECISIONS "${MAGMAMORSE_PRECISION}"
                        TARGETDIR "control" )

    set(MAGMAMORSE_SRCSF
        ../control/morse_f90.f90
        ${MAGMAMORSE_SRCS_F_GENERATED}
       )
endif(HAVE_ISO_C_BINDING)

# Compile step
# ------------
add_library(magmamorse ${MAGMAMORSE_SRCS} ${MAGMAMORSE_SRCSF})

add_dependencies(magmamorse
  magmamorse_include
  coreblas_include
  control_include
  runtime_starpu_include
)

set_property(TARGET magmamorse PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET magmamorse PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")

# installation
# ------------
install(TARGETS magmamorse
        DESTINATION lib)

###
### END CMakeLists.txt
###
