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
 * @file quark_blas.h
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
#ifndef _QUARK_BLAS_H_
#define _QUARK_BLAS_H_

#include "coreblas/include/cblas.h"

#include "runtime/quark/include/quark_zblas.h"
#include "runtime/quark/include/quark_dblas.h"
#include "runtime/quark/include/quark_cblas.h"
#include "runtime/quark/include/quark_sblas.h"
#include "runtime/quark/include/quark_zcblas.h"
#include "runtime/quark/include/quark_dsblas.h"

void CORE_ztile_zero_quark(Quark *quark);
void CORE_dtile_zero_quark(Quark *quark);
void CORE_ctile_zero_quark(Quark *quark);
void CORE_stile_zero_quark(Quark *quark);

#endif
