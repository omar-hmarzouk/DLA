/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file cuda_ztrsm.c
 *
 *  MORSE cudablas kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c d s
 *
 **/
#include "cudablas.h"
#include "cudablas/cudablas_z.h"

int CUDA_ztrsm(MORSE_enum side, MORSE_enum uplo,
               MORSE_enum transa, MORSE_enum diag,
               int m, int n,
               cuDoubleComplex *alpha,
               const cuDoubleComplex *A, int lda,
               cuDoubleComplex *B, int ldb,
               CUBLAS_STREAM_PARAM)
{
    cublasZtrsm(CUBLAS_HANDLE
        morse_cublas_const(side), morse_cublas_const(uplo),
        morse_cublas_const(transa), morse_cublas_const(diag),
        m, n,
        CUBLAS_VALUE(alpha), A, lda,
        B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
