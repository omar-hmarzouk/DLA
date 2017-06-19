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
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
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

/** ****************************************************************************
 * CORE BLAS headers
 **/
#include "coreblas/include/coreblas_z.h"
#include "coreblas/include/coreblas_d.h"
#include "coreblas/include/coreblas_c.h"
#include "coreblas/include/coreblas_s.h"
#include "coreblas/include/coreblas_zc.h"
#include "coreblas/include/coreblas_ds.h"
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 * Coreblas Error
 **/
#define coreblas_error(k, str) do {                                     \
        fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str) ;  \
        assert(0);                                                      \
    } while(0)

/** ****************************************************************************
 * CBlas enum
 **/
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#define CBLAS_UPLO      enum CBLAS_UPLO
#define CBLAS_DIAG      enum CBLAS_DIAG
#define CBLAS_SIDE      enum CBLAS_SIDE

/** ****************************************************************************
 *  LAPACK Constants
 **/
extern char *morse_lapack_constants[];
#define morse_lapack_const(morse_const) morse_lapack_constants[morse_const][0]

void set_coreblas_gemm3m_enabled(int v) ;
int get_coreblas_gemm3m_enabled(void) ;

#ifdef __cplusplus
}
#endif

#endif
