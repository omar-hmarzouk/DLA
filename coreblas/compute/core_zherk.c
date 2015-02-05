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
 * @file core_zherk.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 **/
#include "coreblas/include/coreblas.h"

#undef REAL
#define COMPLEX
#ifdef COMPLEX
/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void CORE_zherk(MORSE_enum uplo, MORSE_enum trans,
                int N, int K,
                double alpha, const MORSE_Complex64_t *A, int LDA,
                double beta, MORSE_Complex64_t *C, int LDC)
{
    cblas_zherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        alpha, A, LDA,
        beta, C, LDC);
}


#endif
