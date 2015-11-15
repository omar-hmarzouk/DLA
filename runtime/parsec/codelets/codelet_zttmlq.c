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
 * @file codelet_zttmlq.c
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

static int
CORE_zttmlq_parsec(dague_execution_unit_t *context, dague_execution_context_t * this_task)
{
    MORSE_enum *side;
    MORSE_enum *trans;
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
    MORSE_Complex64_t *V;
    int *ldv;
    MORSE_Complex64_t *T;
    int *ldt;
    MORSE_Complex64_t *WORK;
    int *ldwork;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &side,
                          UNPACK_VALUE, &trans,
                          UNPACK_VALUE, &m1,
                          UNPACK_VALUE, &n1,
                          UNPACK_VALUE, &m2,
                          UNPACK_VALUE, &n2,
                          UNPACK_VALUE, &k,
                          UNPACK_VALUE, &ib,
                          UNPACK_DATA,  &A1,
                          UNPACK_VALUE, &lda1,
                          UNPACK_DATA,  &A2,
                          UNPACK_VALUE, &lda2,
                          UNPACK_DATA,  &V,
                          UNPACK_VALUE, &ldv,
                          UNPACK_DATA,  &T,
                          UNPACK_VALUE, &ldt,
                          UNPACK_SCRATCH, &WORK,
                          UNPACK_VALUE, &ldwork
                        );

    CORE_zttmlq(*side, *trans, *m1, *n1, *m2, *n2, *k, *ib, A1, *lda1,
                A2, *lda2, V, *ldv, T, *ldt, WORK, *ldwork);

    return 0;
}

void MORSE_TASK_zttmlq(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    int ldwork = side == MorseLeft ? ib : nb;

    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_zttmlq_parsec,        "ttmlq",
                            sizeof(MORSE_enum),     &side,                      VALUE,
                            sizeof(MORSE_enum),     &trans,                     VALUE,
                            sizeof(int),            &m1,                        VALUE,
                            sizeof(int),            &n1,                        VALUE,
                            sizeof(int),            &m2,                        VALUE,
                            sizeof(int),            &n2,                        VALUE,
                            sizeof(int),            &k,                         VALUE,
                            sizeof(int),            &ib,                        VALUE,
                            PASSED_BY_REF,          RTBLKADDR( A1, MORSE_Complex64_t, A1m, A1n ),    INOUT | REGION_FULL,
                            sizeof(int),            &lda1,                      VALUE,
                            PASSED_BY_REF,          RTBLKADDR( A2, MORSE_Complex64_t, A2m, A2n ),    INOUT | REGION_FULL,
                            sizeof(int),            &lda2,                      VALUE,
                            PASSED_BY_REF,          RTBLKADDR( V, MORSE_Complex64_t, Vm, Vn ),       INPUT | REGION_FULL,
                            sizeof(int),            &ldv,                       VALUE,
                            PASSED_BY_REF,          RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ),       INPUT | REGION_FULL,
                            sizeof(int),            &ldt,                       VALUE,
                            sizeof(MORSE_Complex64_t)*ib*nb,    NULL,           SCRATCH,
                            sizeof(int),            &ldwork,                    VALUE,
                             0);
}
