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
 * @file cuda_zherk.c
 *
 *  MORSE cudablas kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @author Florent Pruvost
 * @date 2015-09-17
 * @precisions normal z -> c
 *
 **/
#include "cudablas/include/cudablas.h"

int CUDA_zherk( MORSE_enum uplo, MORSE_enum trans,
                int n, int k,
                double *alpha,
                const cuDoubleComplex *A, int lda,
                double *beta,
                cuDoubleComplex *B, int ldb,
                CUBLAS_STREAM_PARAM)
{
#if !defined(CHAMELEON_USE_CUBLAS_V2)
    cublasSetKernelStream( stream );
#endif

    cublasZherk(
        morse_lapack_const(uplo), morse_lapack_const(trans),
        n, k,
        *alpha, A, lda,
        *beta,  B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
