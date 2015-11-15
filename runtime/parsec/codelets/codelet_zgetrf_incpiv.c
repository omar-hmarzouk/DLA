/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 * @file codelet_zgetrf_incpiv.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

/***************************************************************************//**
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
static int
CORE_zgetrf_incpiv_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *m;
    int *n;
    int *ib;
    MORSE_Complex64_t *A;
    int *lda;
    int *IPIV;
    MORSE_bool *check_info;
    int *iinfo;

    int info;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &m,
                          UNPACK_VALUE, &n,
                          UNPACK_VALUE, &ib,
                          UNPACK_DATA,  &A,
                          UNPACK_VALUE, &lda,
                          UNPACK_SCRATCH, &IPIV,
                          UNPACK_VALUE, &check_info,
                          UNPACK_VALUE, &iinfo
                          );



    CORE_zgetrf_incpiv(*m, *n, *ib, A, *lda, IPIV, &info);

    return 0;
}

void MORSE_TASK_zgetrf_incpiv(MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              MORSE_desc_t *A, int Am, int An, int lda,
                              MORSE_desc_t *L, int Lm, int Ln, int ldl,
                              int *IPIV,
                              MORSE_bool check_info, int iinfo)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_zgetrf_incpiv_parsec,        "getrf_inc",
                             sizeof(int),           &m,                                VALUE,
                             sizeof(int),           &n,                                VALUE,
                             sizeof(int),           &ib,                               VALUE,
                             PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | REGION_FULL,
                             sizeof(int),           &lda,                              VALUE,
                             sizeof(int)*nb,        IPIV,                              SCRATCH,
                             sizeof(int),           &check_info,                       VALUE,
                             sizeof(int),           &iinfo,                            VALUE,
                             0);
}
