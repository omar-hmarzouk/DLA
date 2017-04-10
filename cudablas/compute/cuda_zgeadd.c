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
 * @file cuda_zgeadd.c
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

#if !defined(CHAMELEON_USE_CUBLAS_V2)
#error "This file requires cublas api v2 support"
#endif

int CUDA_zgeadd(MORSE_enum trans,
                int m, int n,
                const cuDoubleComplex *alpha,
                const cuDoubleComplex *A, int lda,
                const cuDoubleComplex *beta,
                cuDoubleComplex *B, int ldb,
                CUBLAS_STREAM_PARAM)
{
    cublasZgeam(CUBLAS_HANDLE
                morse_cublas_const(trans), morse_cublas_const(MorseNoTrans),
                m, n,
                CUBLAS_VALUE(alpha), A, lda,
                CUBLAS_VALUE(beta),  B, ldb,
                B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
