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
 * @file codelet_ztstrf.c
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
CORE_ztstrf_parsec(parsec_execution_stream_t *context, parsec_task_t *this_task)
{
    int *m;
    int *n;
    int *ib;
    int *nb;
    MORSE_Complex64_t *U;
    int *ldu;
    MORSE_Complex64_t *A;
    int *lda;
    MORSE_Complex64_t *L;
    int *ldl;
    int *IPIV;
    MORSE_Complex64_t *WORK;
    int *ldwork;
    MORSE_bool *check_info;
    int *iinfo;

    int info;

    parsec_dtd_unpack_args(
        this_task,
        UNPACK_VALUE,   &m,
        UNPACK_VALUE,   &n,
        UNPACK_VALUE,   &ib,
        UNPACK_VALUE,   &nb,
        UNPACK_DATA,    &U,
        UNPACK_VALUE,   &ldu,
        UNPACK_DATA,    &A,
        UNPACK_VALUE,   &lda,
        UNPACK_DATA,    &L,
        UNPACK_VALUE,   &ldl,
        UNPACK_SCRATCH, &IPIV,
        UNPACK_SCRATCH, &WORK,
        UNPACK_VALUE,   &ldwork,
        UNPACK_VALUE,   &check_info,
        UNPACK_VALUE,   &iinfo );

    CORE_ztstrf(*m, *n, *ib, *nb, U, *ldu, A, *lda, L, *ldl, IPIV, WORK, *ldwork, &info);

    return 0;
}

void MORSE_TASK_ztstrf(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *U, int Um, int Un, int ldu,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztstrf_parsec, options->priority, "tstrf",
        sizeof(int),           &m,                                VALUE,
        sizeof(int),           &n,                                VALUE,
        sizeof(int),           &ib,                               VALUE,
        sizeof(int),           &nb,                               VALUE,
        PASSED_BY_REF,         RTBLKADDR( U, MORSE_Complex64_t, Um, Un ),     INOUT | REGION_FULL,
        sizeof(int),           &ldu,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | REGION_FULL,
        sizeof(int),           &lda,                              VALUE,
        PASSED_BY_REF,         RTBLKADDR( L, MORSE_Complex64_t, Lm, Ln ),     OUTPUT | REGION_FULL,
        sizeof(int),           &ldl,                              VALUE,
        sizeof(int)*nb,        IPIV,                              SCRATCH,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,                 SCRATCH,
        sizeof(int),           &nb,                               VALUE,
        sizeof(int),           &check_info,                       VALUE,
        sizeof(int),           &iinfo,                            VALUE,
        0);
}
