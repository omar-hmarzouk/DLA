/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file cudablas.h
 *
 *  MORSE cudablas headers
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 **/
#ifndef _CUDA_BLAS_H_
#define _CUDA_BLAS_H_

#include "include/chameleon_config.h"

#if !defined(CHAMELEON_USE_CUDA)
#error "This file should not be included"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <cuda.h>
#include <cuComplex.h>

#if defined(CHAMELEON_USE_CUBLAS_V2)

#include <cublas.h>
#include <cublas_v2.h>

#define CUBLAS_STREAM_PARAM cublasHandle_t handle
#define CUBLAS_STREAM_VALUE handle
#define CUBLAS_HANDLE handle,
#define CUBLAS_SADDR(_a_) (&(_a_))
#define CUBLAS_VALUE(_a_) (_a_)
#define CUBLAS_GET_STREAM                       \
    cudaStream_t stream;                        \
    cublasGetStream( handle, &stream )

#else

#include <cublas.h>
#define CUBLAS_STREAM_PARAM cudaStream_t stream
#define CUBLAS_STREAM_VALUE stream
#define CUBLAS_HANDLE
#define CUBLAS_SADDR(_a_) (_a_)
#define CUBLAS_VALUE(_a_) (*(_a_))
#define CUBLAS_GET_STREAM

#endif /* defined(CHAMELEON_USE_CUBLAS_V2) */

#if defined(CHAMELEON_USE_MAGMA)
#include <magma.h>
#endif

/** ****************************************************************************
 * MORSE types and constants
 **/
#include "morse_types.h"
#include "morse_struct.h"
#include "morse_constants.h"
//#include "control/auxiliary.h"
//#include "control/descriptor.h"
//#include "control/tile.h"

/** ****************************************************************************
 * CUDA BLAS headers
 **/
#include "cudablas/include/cudablas_z.h"
#include "cudablas/include/cudablas_d.h"
#include "cudablas/include/cudablas_c.h"
#include "cudablas/include/cudablas_s.h"

/** ****************************************************************************
 * Coreblas Error
 **/
#define cudablas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str)

/** ****************************************************************************
 *  LAPACK Constants
 **/
extern char *morse_lapack_constants[];
#define morse_lapack_const(morse_const) morse_lapack_constants[morse_const][0]

extern int morse_cublas_constants[];

#if defined(CHAMELEON_USE_CUBLAS_V2)
#define morse_cublas_const(morse_const) morse_cublas_constants[morse_const]
#else
#define morse_cublas_const(morse_const) morse_lapack_constants[morse_const][0]
#endif

#endif
