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
 * @file codelet_zgessm.c
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
#include "coreblas/include/cblas.h"
#include "coreblas/include/lapacke.h"
#include "runtime/starpu/include/morse_starpu.h"
#include "runtime/starpu/include/runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgessm applies the factors L computed by CORE_zgetrf_incpiv to
 *  a complex M-by-N tile A.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] K
 *         The number of columns of the tile L. K >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in] IPIV
 *         The pivot indices array of size K as returned by
 *         CORE_zgetrf_incpiv.
 *
 * @param[in] L
 *         The M-by-K lower triangular tile.
 *
 * @param[in] LDL
 *         The leading dimension of the array L.  LDL >= max(1,M).
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, updated by the application of L.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *
 ******************************************************************************/

void MORSE_TASK_zgessm(MORSE_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       MORSE_desc_t *D, int Dm, int Dn, int ldd,
                       MORSE_desc_t *A, int Am, int An, int lda)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgessm;
    void (*callback)(void*) = options->profiling ? cl_zgessm_callback : NULL;

    if ( morse_desc_islocal( L, Lm, Ln ) ||
         morse_desc_islocal( D, Dm, Dn ) ||
         morse_desc_islocal( A, Am, An ) )
    {
        starpu_insert_task(
            codelet,
            STARPU_VALUE,     &m,                        sizeof(int),
            STARPU_VALUE,     &n,                        sizeof(int),
            STARPU_VALUE,     &k,                        sizeof(int),
            STARPU_VALUE,    &ib,                        sizeof(int),
            STARPU_VALUE,          &IPIV,                      sizeof(int*),
            STARPU_R,             RTBLKADDR(L, MORSE_Complex64_t, Lm, Ln),
            STARPU_VALUE,   &ldl,                        sizeof(int),
            STARPU_R,             RTBLKADDR(D, MORSE_Complex64_t, Dm, Dn),
            STARPU_VALUE,   &ldd,                        sizeof(int),
            STARPU_RW,            RTBLKADDR(A, MORSE_Complex64_t, Am, An),
            STARPU_VALUE,   &lda,                        sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
            0);
    }
}


static void cl_zgessm_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    MORSE_Complex64_t *L;
    int ldl;
    MORSE_Complex64_t *D;
    int ldd;
    MORSE_Complex64_t *A;
    int lda;

    L = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    D = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &k, &ib, &IPIV, &ldl, &ldd, &lda);
    CORE_zgessm(m, n, k, ib, IPIV, D, ldd, A, lda);
}

#if defined(CHAMELEON_USE_MAGMA) && defined(HAVE_MAGMA_GETRF_INCPIV_GPU)
static void cl_zgessm_cuda_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int k;
    int ib;
    int *IPIV;
    cuDoubleComplex *dL, *dD, *dA;
    int lddl, lddd, ldda;
    int info = 0;

    int ret;

    /*
     *  hwork => nb*nb
     */
    dL = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    dD = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    dA = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &k, &ib, &IPIV, &lddl, &lddd, &ldda);

    /* The kernel is just using the inverted part or nothing */
    if ( lddl >= 2*ib ) {
      dL += ib;
      ret = magma_zgessm_gpu( MagmaColMajor, m, n, k, ib,
			      IPIV, dL, lddl, dD, lddd, dA, ldda, &info );
    }
    else {
      ret = magma_zgessm_gpu( MagmaColMajor, m, n, k, ib,
			      IPIV, NULL, 1, dD, lddd, dA, ldda, &info );
    }

    if (ret != MAGMA_SUCCESS) {
        fprintf(stderr, "Error in MAGMA: %d\n", ret);
        exit(-1);
    }

    cudaThreadSynchronize();

    return;
}
#endif

/*
 * Codelet definition
 */
#if (defined(CHAMELEON_USE_MAGMA) && defined(HAVE_MAGMA_GETRF_INCPIV_GPU)) || defined(CHAMELEON_SIMULATION)
CODELETS(zgessm, 3, cl_zgessm_cpu_func, cl_zgessm_cuda_func, 0)
#else
CODELETS_CPU(zgessm, 3, cl_zgessm_cpu_func)
#endif
