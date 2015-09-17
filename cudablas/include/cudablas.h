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

#include <stdio.h>
#include <math.h>
#include <string.h>

#if defined(CHAMELEON_USE_CUDA)
#include <cuda.h>
#include <cuComplex.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas_v2.h>
#else
#include <cublas.h>
#endif
#if defined(CHAMELEON_USE_MAGMA)
#include <magma.h>
#endif
/** ****************************************************************************
 * CUDA BLAS headers
 **/
#include "cudablas/include/cudablas_z.h"
#include "cudablas/include/cudablas_d.h"
#include "cudablas/include/cudablas_c.h"
#include "cudablas/include/cudablas_s.h"
#endif

/** ****************************************************************************
 * MORSE types and constants
 **/
#include "morse_types.h"
#include "morse_struct.h"
#include "morse_constants.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *  Global utilities
 **/
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef roundup
#define roundup(a, b) (b <= 0) ? (a) : (((a) + (b)-1) & ~((b)-1))
#endif

/** ****************************************************************************
 *  LAPACK Constants
 **/
extern char *morse_lapack_constants[];
#define morse_lapack_const(morse_const) morse_lapack_constants[morse_const][0]

#ifdef __cplusplus
}
#endif

#endif
