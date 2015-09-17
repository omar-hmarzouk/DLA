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

#if defined(CHAMELEON_USE_MAGMA)
#if defined(CHAMELEON_USE_CUBLAS_V2)
int CUDA_ztrmm_V2(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *C, int ldc,
        CUstream stream)
{
    cublasHandle_t handle;
    cublasStatus_t stat;
    cublasSideMode_t cublasSide;
    cublasFillMode_t cublasUplo;
    cublasOperation_t cublasTransA;
    cublasDiagType_t cublasDiag;

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

    if (side == MorseLeft){
        cublasSide = CUBLAS_SIDE_LEFT;
    }else if (side == MorseRight){
        cublasSide = CUBLAS_SIDE_RIGHT;
    }else{
        fprintf(stderr, "Error in cl_ztrmm_cuda_func: bad side parameter %d\n", side);
    }
    if (uplo == MorseUpper){
        cublasUplo = CUBLAS_FILL_MODE_UPPER;
    }else if(uplo == MorseLower){
        cublasUplo = CUBLAS_FILL_MODE_LOWER;
    }else if(uplo == MorseUpperLower){
        cublasUplo = 0;
    }else{
        fprintf(stderr, "Error in cl_ztrmm_cuda_func: bad uplo parameter %d\n", uplo);
    }
    if (transa == MorseNoTrans){
        cublasTransA = CUBLAS_OP_N;
    }else if(transa == MorseTrans){
        cublasTransa = CUBLAS_OP_T;
    }else if(transa == MorseConjTrans){
        cublasTransa = CUBLAS_OP_C;
    }else{
        fprintf(stderr, "Error in cl_ztrmm_cuda_func: bad transA parameter %d\n", transA);
    }
    if (diag == MorseNonUnit){
        cublasDiag = CUBLAS_DIAG_NON_UNIT;
    }else if(diag == MorseUnit){
        cublasDiag = CUBLAS_DIAG_UNIT;
    }else{
        fprintf(stderr, "Error in cl_ztrmm_cuda_func: bad diag parameter %d\n", diag);
    }

    stat = cublasZtrmm( handle,
        cublasSide, cublasUplo, cublasTransA, cublasDiag,
        m, n,
        (const cuDoubleComplex *) alpha, A, lda,
        B, ldb, C, ldc);
    if (stat != CUBLAS_STATUS_SUCCESS){
        printf ("cublasZtrmm failed");
        cublasDestroy(handle);
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    cublasDestroy(handle);

    return MORSE_SUCCESS;
}
#else /* CHAMELEON_USE_CUBLAS_V2 */
int CUDA_ztrmm(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        CUstream stream)
{
    cublasSetKernelStream( stream );

    cublasZtrmm(
        morse_lapack_const(side), morse_lapack_const(uplo),
        morse_lapack_const(transa), morse_lapack_const(diag),
        m, n,
        *alpha, A, lda,
        B, ldb);

    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

    return MORSE_SUCCESS;
}
#endif /* CHAMELEON_USE_CUBLAS_V2 */
#endif
