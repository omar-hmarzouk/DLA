/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
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

#if defined(CHAMELEON_USE_CUDA)
#if defined(CHAMELEON_USE_CUBLAS_V2)
int CUDA_zher2k_V2(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        double *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream)
{
    cublasHandle_t handle;
    cublasStatus_t stat;
    cublasFillMode_t cublasUplo;
    cublasOperation_t cublasTrans;

    stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    stat = cublasSetStream(handle, stream);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        printf ("cublasSetStream failed\n");
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    if (uplo == MorseUpper){
        cublasUplo = CUBLAS_FILL_MODE_UPPER;
    }else if(uplo == MorseLower){
        cublasUplo = CUBLAS_FILL_MODE_LOWER;
    }else if(uplo == MorseUpperLower){
        cublasUplo = 0;
    }else{
        fprintf(stderr, "Error in cl_zher2k_cuda_func: bad uplo parameter %d\n", uplo);
    }
    if (trans == MorseNoTrans){
        cublasTrans = CUBLAS_OP_N;
    }else if(trans == MorseTrans){
        cublasTrans = CUBLAS_OP_T;
    }else if(trans == MorseConjTrans){
        cublasTrans = CUBLAS_OP_C;
    }else{
        fprintf(stderr, "Error in cl_zher2k_cuda_func: bad trans parameter %d\n", trans);
    }

    stat = cublasZher2k( handle, cublasUplo, cublasTrans,
        n, k, (const cuDoubleComplex *) alpha, A, lda, B, ldb,
        (const double *) beta, C, ldc);
    if (stat != CUBLAS_STATUS_SUCCESS){
        printf ("cublasZher2k failed");
        cublasDestroy(handle);
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    cublasDestroy(handle);

    return MORSE_SUCCESS;
}
#else /* CHAMELEON_USE_CUBLAS_V2 */
int CUDA_zher2k(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        double *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream)
{
    cublasSetKernelStream( stream );

    cublasZher2k(
        morse_lapack_const(uplo), morse_lapack_const(trans),
        n, k,
        *alpha, A, lda,
               B, ldb,
        *beta,  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
#endif /* CHAMELEON_USE_CUBLAS_V2 */
#endif /* CHAMELEON_USE_CUDA */
