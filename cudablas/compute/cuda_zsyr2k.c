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
 * @file cuda_zsyr2k.c
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

int CUDA_zsyr2k(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUBLAS_STREAM_PARAM)
{
#if !defined(CHAMELEON_USE_CUBLAS_V2)
    cublasSetKernelStream( stream );
#endif

    cublasZsyr2k(CUBLAS_HANDLE
                 morse_lapack_const(uplo), morse_lapack_const(trans),
                 n, k,
                 CUBLAS_VALUE(alpha), A, lda,
                                      B, ldb,
                 CUBLAS_VALUE(beta),  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
