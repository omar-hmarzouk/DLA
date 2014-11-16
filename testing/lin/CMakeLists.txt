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
#  @version 0.9.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 13-07-2012
#   
###

cmake_minimum_required(VERSION 2.8)

# Add include and link directories
# --------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# list of sources files
# ---------------------
set(ALINTST
    aladhd.f alaerh.f alaesm.f alahd.f  alareq.f
    alasum.f alasvm.f chkxer.f lsamen.f xlaenv.f xerbla.f
   )

set(SCLNTST slaord.f)

set(DZLNTST dlaord.f)

set(SLINTST 
    schkaa.f
    schkge.f serrge.f sdrvge.f
    serrvx.f
    sget02.f sget04.f sget06.f
    schkpo.f serrpo.f sdrvpo.f
    sposvx.f spotri.f sporfs.f
    spot01.f spot02.f spot03.f spot05.f
    sdrvls.f serrls.f
    schkqr.f serrqr.f
    sqrt01.f sqrt02.f sqrt03.f sqrt13.f sqrt14.f sqrt15.f sqrt16.f sqrt17.f
    schklq.f serrlq.f
    slqt01.f slqt02.f slqt03.f
    slarhs.f slatb4.f sgennd.f
    sgeequ.f slaqge.f
    spoequ.f slaqsy.f spocon.f slacn2.f slatrs.f slabad.f srscl.f
    slascl.f slaset.f slatms.f slartg.f slatm1.f slarnd.f slarot.f
    slaror.f slagsy.f slagge.f slaran.f
    slauum.f strtri.f
    slauu2.f strti2.f
   )

set(DLINTST
    dchkaa.f
    dchkge.f derrge.f ddrvge.f
    derrvx.f
    dget02.f dget04.f dget06.f
    dchkpo.f derrpo.f ddrvpo.f
    dposvx.f dpotri.f dporfs.f
    dpot01.f dpot02.f dpot03.f dpot05.f
    ddrvls.f derrls.f
    dchkqr.f derrqr.f
    dqrt01.f dqrt02.f dqrt03.f dqrt13.f dqrt14.f dqrt15.f dqrt16.f dqrt17.f
    dchklq.f derrlq.f
    dlqt01.f dlqt02.f dlqt03.f
    dlarhs.f dlatb4.f dgennd.f
    dgeequ.f dlaqge.f
    dpoequ.f dlaqsy.f dpocon.f dlacn2.f dlatrs.f dlabad.f drscl.f
    dlascl.f dlaset.f dlatms.f dlartg.f dlatm1.f dlarnd.f dlarot.f
    dlaror.f dlagsy.f dlagge.f dlaran.f
    dlauum.f dtrtri.f
    dlauu2.f dtrti2.f
   )

set(CLINTST
    cchkaa.f
    cchkge.f cerrge.f cdrvge.f
    cerrvx.f
    cget02.f cget04.f
    cchkpo.f cerrpo.f cdrvpo.f
    cposvx.f cpotri.f cporfs.f
    cpot01.f cpot02.f cpot03.f cpot05.f
    cdrvls.f cerrls.f
    cchkqr.f cerrqr.f
    cqrt01.f cqrt02.f cqrt03.f cqrt13.f cqrt14.f cqrt15.f cqrt16.f cqrt17.f
    cchklq.f cerrlq.f
    clqt01.f clqt02.f clqt03.f
    clarhs.f clatb4.f cgennd.f
    cgeequ.f claqge.f
    cpoequ.f claqsy.f cpocon.f clacn2.f clatrs.f csrscl.f
    clascl.f claset.f clatms.f clartg.f clatm1.f clarnd.f clarot.f clanhe.f claqhe.f
    claror.f clagsy.f clagge.f 
    claipd.f csbmv.f
    clauum.f ctrtri.f
    clauu2.f ctrti2.f
    cspmv.f csymv.f
    sget06.f slabad.f slaran.f slarnd.f slascl.f slatm1.f scsum1.f icmax1.f
   )
set(ZLINTST
    zchkaa.f
    zchkge.f zerrge.f zdrvge.f
    zerrvx.f
    zget02.f zget04.f
    zchkpo.f zerrpo.f zdrvpo.f
    zposvx.f zpotri.f zporfs.f
    zpot01.f zpot02.f zpot03.f zpot05.f
    zdrvls.f zerrls.f
    zchkqr.f zerrqr.f
    zqrt01.f zqrt02.f zqrt03.f zqrt13.f zqrt14.f zqrt15.f zqrt16.f zqrt17.f
    zchklq.f zerrlq.f
    zlqt01.f zlqt02.f zlqt03.f
    zlarhs.f zlatb4.f zgennd.f
    zgeequ.f zlaqge.f
    zpoequ.f zlaqsy.f zpocon.f zlacn2.f zlatrs.f zdrscl.f
    zlascl.f zlaset.f zlatms.f zlartg.f zlatm1.f zlarnd.f zlarot.f zlanhe.f zlaqhe.f
    zlaror.f zlagsy.f zlagge.f 
    zlaipd.f zsbmv.f
    zlauum.f ztrtri.f
    zlauu2.f ztrti2.f
    zspmv.f zsymv.f
    dget06.f dlabad.f dlaran.f dlarnd.f dlascl.f dlatm1.f dzsum1.f izmax1.f
   )


# Define what libraries we have to link with
# ------------------------------------------
set(libs_for_tests "chameleon")

# TODO: Check for name of following libraries (it certainly changed between morse and new_chameleon)
if(MORSE_SCHED_STARPU)
    list(APPEND libs_for_tests chameleon_starpu)
    list(APPEND libs_for_tests ${STARPU_LIBRARIES})
    list(APPEND libs_for_tests coreblas)
elseif(MORSE_SCHED_QUARK)
    list(APPEND libs_for_tests chameleon_starpu)
    list(APPEND libs_for_tests ${QUARK_LIBRARIES})
    list(APPEND libs_for_tests coreblas)
endif()

foreach(_dep MAGMA LAPACKE LAPACK CBLAS BLAS CUDA HWLOC MPI DL)
    if(HAVE_${_dep})
        list(APPEND libs_for_tests ${${_dep}_LIBRARY})
    endif()
endforeach()

list(APPEND libs_for_tests ${CMAKE_Fortran_FLAGS})
list(APPEND libs_for_tests ${CMAKE_Fortran_LDFLAGS})


# Define precisions to compile
# ----------------------------
if(BUILD_SINGLE)
    add_executable(morse_xlintsts ${ALINTST} ${SLINTST} ${SCLNTST})
    set_property(TARGET morse_xlintsts PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(morse_xlintsts ${libs_for_tests})
    install(TARGETS morse_xlintsts DESTINATION ${CMAKE_INSTALL_PREFIX}/testing/lin)
endif()

if(BUILD_DOUBLE)
    add_executable(morse_xlintstd ${ALINTST} ${DLINTST} ${DZLNTST})
    set_property(TARGET morse_xlintstd PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(morse_xlintstd ${libs_for_tests})
    install(TARGETS morse_xlintstd DESTINATION ${CMAKE_INSTALL_PREFIX}/testing/lin)
endif()

if(BUILD_COMPLEX)
    add_executable(morse_xlintstc ${ALINTST} ${CLINTST} ${SCLNTST})
    set_property(TARGET morse_xlintstc PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(morse_xlintstc ${libs_for_tests})
    install(TARGETS morse_xlintstc DESTINATION ${CMAKE_INSTALL_PREFIX}/testing/lin)
endif()

if(BUILD_COMPLEX16)
    add_executable(morse_xlintstz ${ALINTST} ${ZLINTST} ${DZLNTST})
    set_property(TARGET morse_xlintstz PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(morse_xlintstz ${libs_for_tests})
    install(TARGETS morse_xlintstz DESTINATION ${CMAKE_INSTALL_PREFIX}/testing/lin)
endif()

# Copy launcher
# -------------
add_custom_target(lin_launcher ALL
                  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/morse_lintesting.py
                                                   ${CMAKE_CURRENT_BINARY_DIR}/morse_lintesting.py)

# Add tests (C/CPP)
# -----------------
add_test(NAME    lapack_testing
         COMMAND ${PYTHON_EXECUTABLE}
                 ${CMAKE_CURRENT_BINARY_DIR}/morse_lintesting.py
        )

###
### END CMakeLists.txt
###