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
 * @file cuda_zmerge.c
 *
 *  MORSE cudablas kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 **/
#include "cudablas/include/cudablas.h"

#if defined(CHAMELEON_USE_MAGMA)
#if defined(CHAMELEON_USE_CUBLAS_V2)
int CUDA_zgemerge(
        magma_side_t side, magma_diag_t diag,
        magma_int_t M, magma_int_t N,
        magmaDoubleComplex *A, magma_int_t LDA,
        magmaDoubleComplex *B, magma_int_t LDB,
        CUstream stream)
{
    int i, j;
    magmaDoubleComplex *cola, *colb;
    cublasHandle_t handle;
    cublasStatus_t stat;

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

    if (M < 0) {
        return -1;
    }
    if (N < 0) {
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        return -7;
    }

    if (side == MagmaLeft){
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(handle, i+1, cola, 1, colb, 1);
            cudaMemcpyAsync(colb , cola,
                            (i+1)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }else{
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(handle, M-i, cola + i, 1, colb + i, 1);
            cudaMemcpyAsync(colb+i , cola+i,
                            (M-i)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }

    cublasDestroy(handle);

    return MORSE_SUCCESS;
}
#else /* CHAMELEON_USE_CUBLAS_V2 */
int CUDA_zgemerge(
        magma_side_t side, magma_diag_t diag,
        magma_int_t M, magma_int_t N,
        magmaDoubleComplex *A, magma_int_t LDA,
        magmaDoubleComplex *B, magma_int_t LDB,
        CUstream stream)
{
    int i, j;
    magmaDoubleComplex *cola, *colb;

    if (M < 0) {
        return -1;
    }
    if (N < 0) {
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        return -7;
    }

    if (side == MagmaLeft){
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(i+1, cola, 1, colb, 1);
            cudaMemcpyAsync(colb , cola,
                            (i+1)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }else{
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(M-i, cola + i, 1, colb + i, 1);
            cudaMemcpyAsync(colb+i , cola+i,
                            (M-i)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }

    return MORSE_SUCCESS;
}
#endif /* CHAMELEON_USE_CUBLAS_V2 */
#endif
