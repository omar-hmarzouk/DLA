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
 * @file codelet_zgessm.c
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
static int
CORE_zgessm_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *m;
    int *n;
    int *k;
    int *ib;
    int *IPIV;
    dague_data_copy_t *gL;
    int *ldl;
    dague_data_copy_t *gD;
    int *ldd;
    dague_data_copy_t *gA;
    int *lda;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &m,
                          UNPACK_VALUE, &n,
                          UNPACK_VALUE, &k,
                          UNPACK_VALUE, &ib,
                          UNPACK_SCRATCH, &IPIV,
                          UNPACK_DATA,  &gL,
                          UNPACK_VALUE, &ldl,
                          UNPACK_DATA,  &gD,
                          UNPACK_VALUE, &ldd,
                          UNPACK_DATA,  &gA,
                          UNPACK_VALUE, &lda
                          );

    //MORSE_Complex64_t *L = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gL);
    MORSE_Complex64_t *D = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gD);
    MORSE_Complex64_t *A = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gA);

    CORE_zgessm(*m, *n, *k, *ib, IPIV, D, *ldd, A, *lda);

    return 0;
}

void MORSE_TASK_zgessm(MORSE_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       MORSE_desc_t *D, int Dm, int Dn, int ldd,
                       MORSE_desc_t *A, int Am, int An, int lda)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_zgessm_parsec,               "gessm",
                             sizeof(int),           &m,                                VALUE,
                             sizeof(int),           &n,                                VALUE,
                             sizeof(int),           &k,                                VALUE,
                             sizeof(int),           &ib,                               VALUE,
                             sizeof(int)*nb,        IPIV,                              SCRATCH,
                             PASSED_BY_REF,         RTBLKADDR( L, MORSE_Complex64_t, Lm, Ln ),     INPUT | REGION_FULL,
                             sizeof(int),           &ldl,                              VALUE,
                             PASSED_BY_REF,         RTBLKADDR( D, MORSE_Complex64_t, Dm, Dn ),     INPUT | REGION_FULL,
                             sizeof(int),           &ldd,                              VALUE,
                             PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | REGION_FULL,
                             sizeof(int),           &lda,                              VALUE,
                             0);
}
