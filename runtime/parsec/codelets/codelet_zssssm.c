/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 * @file codelet_zssssm.c
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
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static int
CORE_zssssm_parsec(parsec_execution_stream_t *context, parsec_task_t * this_task)
{
    int *m1;
    int *n1;
    int *m2;
    int *n2;
    int *k;
    int *ib;
    MORSE_Complex64_t *A1;
    int *lda1;
    MORSE_Complex64_t *A2;
    int *lda2;
    MORSE_Complex64_t *L1;
    int *ldl1;
    MORSE_Complex64_t *L2;
    int *ldl2;
    int *IPIV;

    parsec_dtd_unpack_args(
        this_task,
        UNPACK_VALUE,   &m1,
        UNPACK_VALUE,   &n1,
        UNPACK_VALUE,   &m2,
        UNPACK_VALUE,   &n2,
        UNPACK_VALUE,   &k,
        UNPACK_VALUE,   &ib,
        UNPACK_DATA,    &A1,
        UNPACK_VALUE,   &lda1,
        UNPACK_DATA,    &A2,
        UNPACK_VALUE,   &lda2,
        UNPACK_DATA,    &L1,
        UNPACK_VALUE,   &ldl1,
        UNPACK_DATA,    &L2,
        UNPACK_VALUE,   &ldl2,
        UNPACK_SCRATCH, &IPIV );

    CORE_zssssm(*m1, *n1, *m2, *n2, *k, *ib, A1, *lda1, A2, *lda2, L1, *ldl1, L2, *ldl2, IPIV);

    return 0;
}

void MORSE_TASK_zssssm(const MORSE_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       const MORSE_desc_t *L1, int L1m, int L1n, int ldl1,
                       const MORSE_desc_t *L2, int L2m, int L2n, int ldl2,
                       const int *IPIV)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zssssm_parsec, options->priority, "ssssm",
        sizeof(int),           &m1,                                VALUE,
        sizeof(int),           &n1,                                VALUE,
        sizeof(int),           &m2,                                VALUE,
        sizeof(int),           &n2,                                VALUE,
        sizeof(int),           &k,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A1, MORSE_Complex64_t, A1m, A1n ),     INOUT,
        sizeof(int),           &lda1,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( A2, MORSE_Complex64_t, A2m, A2n ),     INOUT,
        sizeof(int),           &lda2,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( L1, MORSE_Complex64_t, L1m, L1n ),     INPUT,
        sizeof(int),           &ldl1,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( L2, MORSE_Complex64_t, L2m, L2n ),     INPUT,
        sizeof(int),           &ldl2,                              VALUE,
        sizeof(int)*nb,        IPIV,                               SCRATCH,
        0);
}
