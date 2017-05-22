###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
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


# Check for subdirectories
# ------------------------
#add_subdirectory(lin)

# Generate morse auxiliary testing sources for all possible precisions
# --------------------------------------------------------------------
set(TESTING_HDRS_GENERATED "")
set(ZHDR
  testing_zauxiliary.h
 )

precisions_rules_py(TESTING_HDRS_GENERATED "${ZHDR}"
    PRECISIONS "s;d;c;z;ds;zc" )

add_custom_target(testing_include ALL SOURCES ${TESTING_HDRS_GENERATED})

# Generate morse testing sources for all possible precisions
# ----------------------------------------------------------
set(ZSRC
    testing_zauxiliary.c
    ##################
    # BLAS 3
    ##################
    testing_zgemm.c
    testing_zhemm.c
    testing_zher2k.c
    testing_zherk.c
    testing_zpemv.c
    testing_zsymm.c
    testing_zsyr2k.c
    testing_zsyrk.c
    testing_ztrmm.c
    testing_ztrsm.c
    ##################
    # LAPACK
    ##################
    testing_zgels.c
    testing_zgels_hqr.c
    testing_zgels_systolic.c
    #testing_zgesv.c
    testing_zgesv_incpiv.c
    #testing_zgetri.c
    testing_zlange.c
    testing_zposv.c
    testing_zpotri.c
    ##################
    # MIXED PRECISION
    ##################
    #testing_zcgels.c
    #testing_zcgesv.c
    #testing_zcposv.c
    #testing_zcungesv.c
    ##################
    # OTHERS
    ##################
    testing_zgeadd.c
    #testing_zgecfi.c
    testing_zgesvd.c
    #testing_zgetmi.c
    #testing_zheev.c
    testing_zheevd.c
    #testing_zhegst.c
    #testing_zhegv.c
    #testing_zhegvd.c
    testing_zgeqrf_qdwh.c
    )

# Add include and link directories
# --------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Define what libraries we have to link with
# ------------------------------------------
unset(libs_for_tests)
list(APPEND libs_for_tests chameleon)

if(CHAMELEON_SCHED_STARPU)
    list(APPEND libs_for_tests
      chameleon_starpu
      ${STARPU_LIBRARIES_DEP}
    )
    link_directories(${STARPU_LIBRARY_DIRS})
elseif(CHAMELEON_SCHED_PARSEC)
    list(APPEND libs_for_tests
        chameleon_parsec
        ${PARSEC_LIBRARIES_DEP}
      )
    link_directories(${PARSEC_LIBRARY_DIRS_DEP})
elseif(CHAMELEON_SCHED_QUARK)
    list(APPEND libs_for_tests
        chameleon_quark
        ${QUARK_LIBRARIES}
      )
    link_directories(${QUARK_LIBRARY_DIRS})
endif()

if(NOT CHAMELEON_SIMULATION)

    if(CHAMELEON_USE_CUDA)
        list(APPEND libs_for_tests
        cudablas)
    endif()
    if(CHAMELEON_USE_CUDA)
        link_directories(${CUDA_LIBRARY_DIRS})
    endif()

    list(APPEND libs_for_tests
    coreblas
    ${LAPACKE_LIBRARIES}
    ${TMG_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${CBLAS_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${HWLOC_LIBRARIES}
    ${LIBHQR_LIBRARIES}
    ${EXTRA_LIBRARIES}
    )

    link_directories(${LAPACKE_LIBRARY_DIRS})
    link_directories(${TMG_LIBRARY_DIRS})
    link_directories(${LAPACK_LIBRARY_DIRS})
    link_directories(${CBLAS_LIBRARY_DIRS})
    link_directories(${BLAS_LIBRARY_DIRS})

else()

    list(APPEND libs_for_tests ${EXTRA_LIBRARIES})

endif()

# message("LAPACKE_LIBRARY_DIRS: ${LAPACKE_LIBRARY_DIRS}")
# message("LAPACK_LIBRARY_DIRS: ${LAPACK_LIBRARY_DIRS}")
# message("CBLAS_LIBRARY_DIRS: ${CBLAS_LIBRARY_DIRS}")
# message("BLAS_LIBRARY_DIRS: ${BLAS_LIBRARY_DIRS}")
# message("HWLOC_LIBRARY_DIRS: ${HWLOC_LIBRARY_DIRS}")

# list(APPEND libs_for_tests ${CMAKE_Fortran_FLAGS}  )
# list(APPEND libs_for_tests ${CMAKE_Fortran_LDFLAGS})

# message(STATUS "libs testings: ${libs_for_tests}")

# Create one testing per precision with all files
# -----------------------------------------------
foreach(_precision ${CHAMELEON_PRECISION} )

    precisions_rules_py(${_precision}SRC_GENERATED "${ZSRC}"
        PRECISIONS "${_precision}" )

    add_executable(${_precision}testing ${${_precision}SRC_GENERATED})
    add_dependencies(${_precision}testing
        chameleon_include
        coreblas_include
        control_include
        testing_include
      )
    set_property(TARGET ${_precision}testing PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_precision}testing ${libs_for_tests})

    install(TARGETS ${_precision}testing
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/testing)

endforeach()

# Copy python scripts to use test drivers
# ---------------------------------------
if(CHAMELEON_SCHED_STARPU)
    # Copy launcher
    add_custom_target(testing_launcher_starpu ALL
                     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/chameleon_testing.py
                                                      ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py)
    # install file
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/testing)
endif()

if(CHAMELEON_SCHED_PARSEC)
    # Copy launcher
    add_custom_target(testing_launcher_parsec ALL
                     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/chameleon_testing.py
                                                      ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py)
    # install file
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/testing)
endif()

if(CHAMELEON_SCHED_QUARK)
    # Copy launcher
    add_custom_target(testing_launcher_quark ALL
                     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/chameleon_testing.py
                                                      ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py)
    # install file
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/chameleon_testing.py
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/chameleon/testing)
endif()

#-------- Tests ---------
include(CTestLists.cmake)

###
### END CMakeLists.txt
###
