/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_zgelqt.c
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

#include "runtime/starpu/include/morse_starpu.h"
#include "runtime/starpu/include/runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgelqt - computes a LQ factorization of a complex M-by-N tile A: A = L * Q.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
 *  A(i,i+1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-min(M,N) lower trapezoidal tile L (L is
 *         lower triangular if M <= N); the elements above the diagonal,
 *         with the array TAU, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/

void MORSE_TASK_zgelqt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgelqt;
    void (*callback)(void*) = options->profiling ? cl_zgelqt_callback : NULL;
    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);

    if ( morse_desc_islocal( A, Am, An ) ||
         morse_desc_islocal( T, Tm, Tn ) )
    {
        starpu_insert_task(
            codelet,
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_VALUE,    &ib,                sizeof(int),
            STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_W,         RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),
            STARPU_VALUE,    &ldt,               sizeof(int),
            /* max( nb * (ib+1), ib * (ib+nb) ) */
            STARPU_SCRATCH,   options->ws_worker,
            /* /\* ib*n + 3*ib*ib + max(m,n) *\/ */
            STARPU_VALUE,    &h_work,            sizeof(MORSE_starpu_ws_t *),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
            0);
    }
}


static void cl_zgelqt_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *TAU, *WORK;

    A   = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    T   = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    TAU = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]); /* max(m,n) + ib*n */

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &ldt, &h_work);

    WORK = TAU + max( m, n );
    CORE_zgelqt(m, n, ib, A, lda, T, ldt, TAU, WORK);
}

#if defined(CHAMELEON_USE_MAGMA)
static void cl_zgelqt_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    cuDoubleComplex *h_A, *h_T, *h_D, *h_W, *h_TAU;
    cuDoubleComplex *d_A, *d_T, *d_D, *d_W;
    int lda, ldt;
    CUstream stream;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &ldt, &h_work);

    /* Gather pointer to data on device */
    d_A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    d_T = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    d_W = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]); /* m*ib + ib*ib*/
    d_D = d_W + m*ib;

    /* scratch data on host */
    /* ib*n + ib*ib + max(m,n) + ib*ib + ib*ib */
    h_A = (cuDoubleComplex*)RUNTIME_starpu_ws_getlocal(h_work);

    /* Gather pointer to scratch data on host */
    h_T   = h_A   + ib*n;
    h_TAU = h_T   + ib*ib;
    h_W   = h_TAU + max(m,n);
    h_D   = h_W   + ib*ib;

    stream = starpu_cuda_get_local_stream();
    cublasSetKernelStream( stream );

    CUDA_zgelqt(
            m, n, ib,
            d_A, lda, h_A, ib,
            d_T, ldt, h_T, ib,
            d_D, h_D, ib, h_TAU,
            h_W, d_W, stream);

    cudaThreadSynchronize();
}

#endif

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_MAGMA) || defined(CHAMELEON_SIMULATION)
CODELETS(zgelqt, 3, cl_zgelqt_cpu_func, cl_zgelqt_cuda_func, 0)
#else
CODELETS_CPU(zgelqt, 3, cl_zgelqt_cpu_func)
#endif
