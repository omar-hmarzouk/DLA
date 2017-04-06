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
 * @file cuda_ztrmm.c
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
#include "cudablas/include/cudablas.h"
#include "cudablas/include/cudablas_z.h"

int CUDA_ztrmm(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        CUBLAS_STREAM_PARAM)
{
    cublasZtrmm(CUBLAS_HANDLE
        morse_lapack_const(side), morse_lapack_const(uplo),
        morse_lapack_const(transa), morse_lapack_const(diag),
        m, n,
        CUBLAS_VALUE(alpha), A, lda,
        B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}

