/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file coreblas.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @date 2010-11-15
 *
 **/
#ifndef _CORE_BLAS_H_
#define _CORE_BLAS_H_

#include <stdio.h>
#include <math.h>
#include <string.h>

/** ****************************************************************************
 *  CBLAS requires for scalar arguments to be passed
 *        by address rather than by value
 **/
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif
#include "coreblas/include/cblas.h"

/** ****************************************************************************
 * MORSE types and constants
 **/
#include "morse_types.h"
#include "morse_struct.h"
#include "morse_constants.h"
//#include "control/auxiliary.h"
//#include "control/descriptor.h"
//#include "control/tile.h"
//#include "control/bulge.h"

/** ****************************************************************************
 * CORE BLAS headers
 **/
#include "coreblas/include/coreblas_z.h"
#include "coreblas/include/coreblas_d.h"
#include "coreblas/include/coreblas_c.h"
#include "coreblas/include/coreblas_s.h"
#include "coreblas/include/coreblas_zc.h"
#include "coreblas/include/coreblas_ds.h"

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 * Coreblas Error
 **/
#define coreblas_error(k, str) fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str)

/** ****************************************************************************
 * CBlas enum
 **/
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#define CBLAS_UPLO      enum CBLAS_UPLO
#define CBLAS_DIAG      enum CBLAS_DIAG
#define CBLAS_SIDE      enum CBLAS_SIDE

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
