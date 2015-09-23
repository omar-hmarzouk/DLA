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
 * @file cuda_zgemm.c
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

#if defined(CHAMELEON_USE_CUDA)
#if defined(CHAMELEON_USE_CUBLAS_V2)
int CUDA_zgemm_V2(
        MORSE_enum transa, MORSE_enum transb,
        int m, int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream)
{
    cublasHandle_t handle;
    cublasStatus_t stat;
    cublasOperation_t cublasTransA;
    cublasOperation_t cublasTransB;

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

    if (transa == MorseNoTrans){
        cublasTransA = CUBLAS_OP_N;
    }else if(transa == MorseTrans){
        cublasTransA = CUBLAS_OP_T;
    }else if(transa == MorseConjTrans){
        cublasTransA = CUBLAS_OP_C;
    }else{
        fprintf(stderr, "Error in cl_zgemm_cuda_func: bad transA parameter %d\n", transa);
    }
    if (transb == MorseNoTrans){
        cublasTransB = CUBLAS_OP_N;
    }else if(transb == MorseTrans){
        cublasTransB = CUBLAS_OP_T;
    }else if(transb == MorseConjTrans){
        cublasTransB = CUBLAS_OP_C;
    }else{
        fprintf(stderr, "Error in cl_zgemm_cuda_func: bad transB parameter %d\n", transb);
    }

    stat = cublasZgemm(handle,
        cublasTransA, cublasTransB,
        m, n, k,
        (const cuDoubleComplex *) alpha, A, lda,
        B, ldb,
        (const cuDoubleComplex *) beta,  C, ldc);
    if (stat != CUBLAS_STATUS_SUCCESS){
        printf ("cublasZgemm failed");
        cublasDestroy(handle);
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    cublasDestroy(handle);

    return MORSE_SUCCESS;
}
#else /* CHAMELEON_USE_CUBLAS_V2 */
int CUDA_zgemm(
        int transa, int transb,
        int m, int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream
        )
{

    cublasSetKernelStream( stream );

    cublasZgemm(
        morse_lapack_const(transa), morse_lapack_const(transb),
        m, n, k,
        *alpha, A, lda,
                B, ldb,
        *beta,  C, ldc);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
#endif /* CHAMELEON_USE_CUBLAS_V2 */
#endif /* CHAMELEON_USE_CUDA */
