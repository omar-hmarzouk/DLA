/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014, 2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_ztslqt.c
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
#undef REAL
#define COMPLEX

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_ztslqt computes a LQ factorization of a rectangular matrix
 *  formed by coupling side-by-side a complex M-by-M
 *  lower triangular tile A1 and a complex M-by-N tile A2:
 *
 *    | A1 A2 | = L * Q
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
 *  A2(i,1:n), and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *         The number of rows of the tile A1 and A2. M >= 0.
 *         The number of columns of the tile A1.
 *
 * @param[in] N
 *         The number of columns of the tile A2. N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A1
 *         On entry, the M-by-M tile A1.
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-M lower trapezoidal tile L;
 *         the elements above the diagonal are not referenced.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, all the elements with the array TAU, represent
 *         the unitary tile Q as a product of elementary reflectors
 *         (see Further Details).
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
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

void MORSE_TASK_ztslqt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_ztslqt;
    void (*callback)(void*) = options->profiling ? cl_ztslqt_callback : NULL;
    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);

    if ( morse_desc_islocal( A1, A1m, A1n ) ||
         morse_desc_islocal( A2, A2m, A2n ) ||
         morse_desc_islocal( T,  Tm,  Tn  ) )
    {
        starpu_insert_task(
            codelet,
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_VALUE,    &ib,                sizeof(int),
            STARPU_RW,        RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n),
            STARPU_VALUE,    &lda1,              sizeof(int),
            STARPU_RW,        RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n),
            STARPU_VALUE,    &lda2,              sizeof(int),
            STARPU_W,         RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn),
            STARPU_VALUE,    &ldt,               sizeof(int),
            /* max( nb * (ib+1), ib * (ib+nb) ) */
            STARPU_SCRATCH,   options->ws_worker,
            /* /\* 2 * ib * (nb+ib) + nb *\/ */
            STARPU_VALUE,    &h_work,            sizeof(MORSE_starpu_ws_t *),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "ztslqt",
#endif
            0);
    }
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztslqt_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A1;
    int lda1;
    MORSE_Complex64_t *A2;
    int lda2;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *TAU, *WORK;

    A1 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    A2 = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    T  = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    TAU= (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[3]); /* nb + ib*nb */

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda1, &lda2, &ldt, &h_work);

    WORK = TAU + max( m, n );
    CORE_ztslqt(m, n, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK);
}

#if defined(CHAMELEON_USE_MAGMA) && 0
static void cl_ztslqt_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    cuDoubleComplex *h_A2, *h_T, *h_D, *h_TAU, *h_W;
    cuDoubleComplex *d_A1, *d_A2, *d_T, *d_D, *d_W;
    int lda1, lda2, ldt;
    CUstream stream;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda1, &lda2, &ldt, &h_work);

    /* Gather pointer to data on device */
    d_A1 = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    d_A2 = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    d_T  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    d_W  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[3]); /* 2*ib*m + ib*ib */
    d_D  = d_W + 2*ib*m;

    /* scratch data on host */
    /* ib*n + ib*n + max(m,n) + ib*m + ib*ib */
    h_A2  = (cuDoubleComplex*)RUNTIME_starpu_ws_getlocal(h_work);
    h_T   = h_A2  + ib*n;
    h_TAU = h_T   + ib*n;
    h_W   = h_TAU + max(m,n);
    h_D   = h_W   + ib*m;

    stream = starpu_cuda_get_local_stream();
    CUDA_ztslqt(
            m, n, ib,
            d_A1, lda1, d_A2, lda2,
            h_A2, ib,
            d_T, ldt, h_T, ib,
            d_D, h_D, ib, h_TAU,
            h_W, d_W, stream);

    cudaThreadSynchronize();
}
#endif
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if (defined(CHAMELEON_USE_MAGMA)) && 0
CODELETS(ztslqt, 4, cl_ztslqt_cpu_func, cl_ztslqt_cuda_func, 0)
#else
CODELETS_CPU(ztslqt, 4, cl_ztslqt_cpu_func)
#endif
