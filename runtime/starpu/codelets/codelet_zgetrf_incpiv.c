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
 * @file codelet_zgetrf_incpiv.c
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
#include "coreblas/include/lapacke.h"
#include "runtime/starpu/include/morse_starpu.h"
#include "runtime/starpu/include/runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgetrf_incpiv computes an LU factorization of a general M-by-N tile A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *
 *    A = P * L * U
 *
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *  This is the right-looking Level 2.5 BLAS version of the algorithm.
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
 *         On entry, the M-by-N tile to be factored.
 *         On exit, the factors L and U from the factorization
 *         A = P*L*U; the unit diagonal elements of L are not stored.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] IPIV
 *         The pivot indices; for 1 <= i <= min(M,N), row i of the
 *         tile was interchanged with row IPIV(i).
 *
 * @param[out] INFO
 *         See returned value.
 *
 *******************************************************************************
 *
 * @return
 *         \retval MORSE_SUCCESS successful exit
 *         \retval <0 if INFO = -k, the k-th argument had an illegal value
 *         \retval >0 if INFO = k, U(k,k) is exactly zero. The factorization
 *              has been completed, but the factor U is exactly
 *              singular, and division by zero will occur if it is used
 *              to solve a system of equations.
 *
 ******************************************************************************/

void MORSE_TASK_zgetrf_incpiv(MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              MORSE_desc_t *A, int Am, int An, int lda,
                              MORSE_desc_t *L, int Lm, int Ln, int ldl,
                              int *IPIV,
                              MORSE_bool check_info, int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf_incpiv;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_incpiv_callback : NULL;

    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);

    if ( morse_desc_islocal( A, Am, An ) ||
         morse_desc_islocal( L, Lm, Ln ) )
    {
        starpu_insert_task(
            codelet,
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_VALUE,    &ib,                sizeof(int),
            STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_W,         RTBLKADDR(L, MORSE_Complex64_t, Lm, Ln),
            STARPU_VALUE,    &ldl,               sizeof(int),
            STARPU_VALUE,    &IPIV,              sizeof(int*),
            STARPU_VALUE,    &check_info,        sizeof(MORSE_bool),
            STARPU_VALUE,    &iinfo,             sizeof(int),
            STARPU_SCRATCH,   options->ws_worker,
            STARPU_VALUE,    &h_work,            sizeof(MORSE_starpu_ws_t *),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
            0);
    }
}


static void cl_zgetrf_incpiv_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_starpu_ws_t *h_work;
    int m;
    int n;
    int ib;
    MORSE_Complex64_t *A, *L;
    int lda, ldl;
    int *IPIV;
    MORSE_bool check_info;
    int iinfo;

    int info = 0;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    L = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &ldl, &IPIV, &check_info, &iinfo, &h_work);
    CORE_zgetrf_incpiv(m, n, ib, A, lda, IPIV, &info);

#if defined(CHAMELEON_USE_MAGMA)
    /*
     * L stores:
     *      L1     L2    L3     ...
     *      L1^-1  L2^-1 L3^-1  ...
     */
    /* Compute L-1 in lower rectangle of L */
    if ( ldl >= 2*ib )
    {
        int i, sb;

	L += ib;
        for (i=0; i<n; i+=ib) {
            sb = min( ib, n-i );
            CORE_zlacpy(MorseUpperLower, sb, sb, A+(i*lda+i), lda, L+(i*ldl), ldl );

            CORE_ztrtri( MorseLower, MorseUnit, sb, L+(i*ldl), ldl, &info );
            if (info != 0 ) {
                fprintf(stderr, "ERROR, trtri returned with info = %d\n", info);
            }
        }
    }
#endif
}


/*
 * Codelet GPU
 */
#if defined(CHAMELEON_USE_MAGMA) && defined(HAVE_MAGMA_GETRF_INCPIV_GPU)
static void cl_zgetrf_incpiv_cuda_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int ib;
    cuDoubleComplex *hA, *dA;
    cuDoubleComplex *hL, *dL;
    cuDoubleComplex *dwork;
    MORSE_starpu_ws_t *h_work;
    int lda, ldl;
    int *IPIV;
    MORSE_bool check_info;
    int iinfo;
    int info;

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &lda, &ldl, &IPIV, &check_info, &iinfo, &h_work);

    dA = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    dL = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    /*
     * hwork => at least (IB+NB)*IB contains all hA and hL
     * dwork => at least IB*NB
     */
    hA    = (cuDoubleComplex*)RUNTIME_starpu_ws_getlocal(h_work);
    dwork = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    hL = hA + lda*ib;

    /* Initialize L to 0 */
    memset(hL, 0, ib*ib*sizeof(cuDoubleComplex));

    if ( ldl >= 2*ib ) {
      /* Let's compute the inverses in the bottom part of L */
      dL += ib;
    } else {
      /* We prefer to stick with TRSM */
      dL = NULL;
      hL = NULL;
    }

    magma_zgetrf_incpiv_gpu( MagmaColMajor, m, n, ib,
                             hA, lda, dA, lda,
                             hL, ib,  dL, ldl,
                             IPIV,
                             dwork, lda,
                             &info );

    cudaThreadSynchronize();
}
#endif


/*
 * Codelet definition
 */
#if (defined(CHAMELEON_USE_MAGMA) && defined(HAVE_MAGMA_GETRF_INCPIV_GPU)) || defined(CHAMELEON_SIMULATION)
CODELETS(zgetrf_incpiv, 3, cl_zgetrf_incpiv_cpu_func, cl_zgetrf_incpiv_cuda_func, 0)
#else
CODELETS_CPU(zgetrf_incpiv, 3, cl_zgetrf_incpiv_cpu_func)
#endif
