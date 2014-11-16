/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University 
 *                          of Tennessee Research Foundation. 
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 IPB. All rights reserved. 
 *
 **/

/**
 *
 * @file codelet_zgemm.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "morse_starpu.h"
#include "codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_zgemm(MORSE_option_t *options,
                      MORSE_enum transA, int transB,
                      int m, int n, int k, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                                                MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgemm;
    void (*callback)(void*) = options->profiling ? cl_zgemm_callback : NULL;

    if ( morse_desc_islocal( A, Am, An ) ||
         morse_desc_islocal( B, Bm, Bn ) ||
         morse_desc_islocal( C, Cm, Cn ) )
    {
        starpu_insert_task(
            codelet,
            STARPU_VALUE,    &transA,            sizeof(MORSE_enum),
            STARPU_VALUE,    &transB,            sizeof(MORSE_enum),
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_VALUE,    &k,                 sizeof(int),
            STARPU_VALUE,    &alpha,             sizeof(MORSE_Complex64_t),
            STARPU_R,         RTBLKADDR(A, MORSE_Complex64_t, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_R,         RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),
            STARPU_VALUE,    &ldb,               sizeof(int),
            STARPU_VALUE,    &beta,              sizeof(MORSE_Complex64_t),
            STARPU_RW,        RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn),
            STARPU_VALUE,    &ldc,               sizeof(int),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
            0);
    }
}



static void cl_zgemm_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum transA;
    MORSE_enum transB;
    int m;
    int n;
    int k;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *C;
    int ldc;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        m, n, k,
        CBLAS_SADDR(alpha), A, lda,
        B, ldb,
        CBLAS_SADDR(beta), C, ldc);
}

#ifdef CHAMELEON_USE_CUDA
static void cl_zgemm_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum transA;
    MORSE_enum transB;
    int m;
    int n;
    int k;
    cuDoubleComplex alpha;
    cuDoubleComplex *A;
    int lda;
    cuDoubleComplex *B;
    int ldb;
    cuDoubleComplex beta;
    cuDoubleComplex *C;
    int ldc;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);

    CUstream stream = starpu_cuda_get_local_stream();
    cublasSetKernelStream( stream );

    cublasZgemm(
        lapack_const(transA), lapack_const(transB),
        m, n, k,
        alpha, A, lda,
               B, ldb,
        beta,  C, ldc);
    assert( CUBLAS_STATUS_SUCCESS == cublasGetError() );

#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif

    return;
}
#endif

/*
 * Codelet definition
 */
CODELETS(zgemm, 3, cl_zgemm_cpu_func, cl_zgemm_cuda_func, STARPU_CUDA_ASYNC);