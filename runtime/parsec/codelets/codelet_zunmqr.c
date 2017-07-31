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
 * @file codelet_zunmqr.c
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

static int
CORE_zunmqr_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    MORSE_enum *side;
    MORSE_enum *trans;
    int *m;
    int *n;
    int *k;
    int *ib;
    MORSE_Complex64_t *A;
    int *lda;
    MORSE_Complex64_t *T;
    int *ldt;
    MORSE_Complex64_t *C;
    int *ldc;
    MORSE_Complex64_t *WORK;
    int *ldwork;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE,   &side,
        UNPACK_VALUE,   &trans,
        UNPACK_VALUE,   &m,
        UNPACK_VALUE,   &n,
        UNPACK_VALUE,   &k,
        UNPACK_VALUE,   &ib,
        UNPACK_DATA,    &A,
        UNPACK_VALUE,   &lda,
        UNPACK_DATA,    &T,
        UNPACK_VALUE,   &ldt,
        UNPACK_DATA,    &C,
        UNPACK_VALUE,   &ldc,
        UNPACK_SCRATCH, &WORK,
        UNPACK_VALUE,   &ldwork );

    CORE_zunmqr(*side, *trans, *m, *n, *k, *ib,
                A, *lda, T, *ldt, C, *ldc, WORK, *ldwork);

    return 0;
}

void MORSE_TASK_zunmqr(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zunmqr_parsec, "unmqr",
        sizeof(MORSE_enum),    &side,                              VALUE,
        sizeof(MORSE_enum),    &trans,                             VALUE,
        sizeof(int),           &m,                                 VALUE,
        sizeof(int),           &n,                                 VALUE,
        sizeof(int),           &k,                                 VALUE,
        sizeof(int),           &ib,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT | REGION_L,
        sizeof(int),           &lda,                               VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ),     INPUT | REGION_FULL,
        sizeof(int),           &ldt,                               VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, MORSE_Complex64_t, Cm, Cn ),     INOUT | REGION_FULL,
        sizeof(int),           &ldc,                               VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,   NULL,                   SCRATCH,
        sizeof(int),           &nb,                                VALUE,
        0);
}
