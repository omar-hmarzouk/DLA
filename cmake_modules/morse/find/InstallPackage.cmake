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
# - Install a library
# This module installs a version of BLAS library from either:
#   a tarball downloaded (requires an internet access)
#   or from sources already present on the system
#
# We handle two different modes:
#
# - Install from a tarball downloaded on internet (we decide the distribution to ensure compatibility)
#   - activate it with BLAS_DOWNLOAD=ON
#   - will download and build NETLIB BLAS: http://netlib.org/blas/blas.tgz
#   - the resulting library is blas_LINUX.a, see their make.inc !
#
# - Install from sources somewhere on the system: 
#   - activate it with BLAS_SOURCE_DIR=path/to/sources, note you manage your build configuration through the make.inc 
#   - we check for blas_LINUX.a library by default, set BLAS_SOURCE_LIBNAME=your_blas_lib_name if different
#   - for now we handle the NETLIB blas only, see http://netlib.org/blas/blas.tgz, maybe more in the future
#
# Note: if BLAS_DOWNLOAD=ON then BLAS_SOURCE mode is not considered
#       turn it OFF to install sources using BLAS_SOURCE_DIR and optionnaly BLAS_SOURCE_LIBNAME

# List of libraries handle for now
set(MORSE_INSTALL_LIBRARY_LIST
    BLAS
    CBLAS
    CHAMELEON
    FFTW
    FXT
    HWLOC
    LAPACK
    LAPACKE
    MAGMA
    METIS
    MUMPS
    PARMETIS
    PASTIX
    PTSCOTCH
    QUARK
    SCALAPACK
    SCOTCH
    STARPU
    TMG
    )

# Set list if URLs to libraries
set(BLAS_TAR_URL      "http://netlib.org/blas/blas.tgz")
set(CBLAS_TAR_URL     "http://www.netlib.org/blas/blast-forum/cblas.tgz")
set(CHAMELEON_TAR_URL "")
set(FFTW_TAR_URL      "")
set(FXT_TAR_URL       "")
set(HWLOC_TAR_URL     "")
set(LAPACK_TAR_URL    "")
set(LAPACKE_TAR_URL   "")
set(MAGMA_TAR_URL     "")
set(METIS_TAR_URL     "")
set(MUMPS_TAR_URL     "")
set(PARMETIS_TAR_URL  "")
set(PASTIX_TAR_URL    "")
set(PTSCOTCH_TAR_URL  "")
set(QUARK_TAR_URL     "")
set(SCALAPACK_TAR_URL "")
set(SCOTCH_TAR_URL    "")
set(STARPU_TAR_URL    "")
set(TMG_TAR_URL       "")

# Set buildtool for each library
set(BLAS_BUILDTOOL      "make")
set(CBLAS_BUILDTOOL     "make")
set(CHAMELEON_BUILDTOOL "cmake")
set(FFTW_BUILDTOOL      "")
set(FXT_BUILDTOOL       "autotools")
set(HWLOC_BUILDTOOL     "autotools")
set(LAPACK_BUILDTOOL    "cmake")
set(LAPACKE_BUILDTOOL   "cmake")
set(MAGMA_BUILDTOOL     "make")
set(METIS_BUILDTOOL     "make")
set(MUMPS_BUILDTOOL     "make")
set(PARMETIS_BUILDTOOL  "make")
set(PASTIX_BUILDTOOL    "cmake")
set(PTSCOTCH_BUILDTOOL  "make")
set(QUARK_BUILDTOOL     "make")
set(SCALAPACK_BUILDTOOL "make")
set(SCOTCH_BUILDTOOL    "make")
set(STARPU_BUILDTOOL    "autotools")
set(TMG_BUILDTOOL       "cmake")

# Macro to install a library
macro(install_package _libname)

message(STATUS "Install package ${_libname}")

# _libname must be one we manage, check the list: LIBRARY_LIST
if (NOT MORSE_INSTALL_LIBRARY_LIST MATCHES "${_libname}")
    message(FATAL_ERROR "First argument for library name, ${_libname},"
    " is not managed in Morse for now. Check MORSE_INSTALL_LIBRARY_LIST list.")
endif()

message("ARGV: ${ARGV}")

# if download mode: call the corresponding macro
if (ARGV MATCHES "DOWNLOAD")

    # Download the tarball
    message(STATUS "Download: ${${_libname}_TAR_URL}")

    get_filename_component(tarball_name "${${_libname}_TAR_URL}" NAME)

    if (NOT EXISTS "${CMAKE_SOURCE_DIR}/externals/${tarball_name}")
        file(DOWNLOAD ${${_libname}_TAR_URL}
             ${CMAKE_SOURCE_DIR}/externals/${tarball_name}
             STATUS IS_GOT
             SHOW_PROGRESS
             TIMEOUT 30
        )
        if (EXISTS "${CMAKE_SOURCE_DIR}/externals/${tarball_name}")
            message(STATUS "${tarball_name} downloaded: ${CMAKE_SOURCE_DIR}/externals/${tarball_name}")
        else()
            message(WARNING "${${_libname}_TAR_URL} download has failed")
        endif()
    endif()

    # Untar
    if (NOT EXISTS "${CMAKE_SOURCE_DIR}/externals/${_libname}/")
        message(STATUS "Untar ${tarball_name}")
        execute_process(
            COMMAND tar xvf ${CMAKE_SOURCE_DIR}/externals/${tarball_name} ${_libname}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/externals
        )
        if (EXISTS "${CMAKE_SOURCE_DIR}/externals/${_libname}/")
            message(STATUS "${_libname} untared: ${CMAKE_SOURCE_DIR}/externals/${_libname}/")
        else()
            message(WARNING "${_libname} untar has failed")
        endif()
    endif()

endif()


endmacro(install_package)

# Macro to install a library from a downloaded tarball
#macro(install_library_from_download)
#endmacro()

