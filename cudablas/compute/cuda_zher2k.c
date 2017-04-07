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
 * @file cuda_zher2k.c
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
#include "cudablas/include/cudablas_z.h"

int CUDA_zher2k(MORSE_enum uplo, MORSE_enum trans,
                int n, int k,
                cuDoubleComplex *alpha,
                const cuDoubleComplex *A, int lda,
                const cuDoubleComplex *B, int ldb,
                double *beta,
                cuDoubleComplex *C, int ldc,
                CUBLAS_STREAM_PARAM)
{
    cublasZher2k(CUBLAS_HANDLE
                 morse_cublas_const(uplo), morse_cublas_const(trans),
                 n, k,
                 CUBLAS_VALUE(alpha), A, lda,
                                      B, ldb,
                 CUBLAS_VALUE(beta),  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
