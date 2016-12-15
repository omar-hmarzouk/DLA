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
 * @file codelets_ztsmqr_hetra1.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

static inline int
CORE_ztsmqr_hetra1_parsec(dague_execution_unit_t    *context,
                          dague_execution_context_t *this_task)
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
                          UNPACK_VALUE,   &side,
                          UNPACK_VALUE,   &trans,
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
                          UNPACK_DATA,    &V,
                          UNPACK_VALUE,   &ldv,
                          UNPACK_DATA,    &T,
                          UNPACK_VALUE,   &ldt,
                          UNPACK_SCRATCH, &WORK,
                          UNPACK_VALUE,   &ldwork);

    CORE_ztsmqr_hetra1(*side, *trans, *m1, *n1, *m2, *n2, *k, *ib,
                       A1, *lda1, A2, *lda2,
                       V, *ldv, T, *ldt,
                       WORK, *ldwork);
}

void MORSE_TASK_ztsmqr_hetra1(const MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              const MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);
    int ldwork = side == MorseLeft ? ib : nb;

    dague_insert_task(
        DAGUE_dtd_handle, CORE_ztsmqr_hetra1_parsec, "tsmqr_hetra1",
        sizeof(MORSE_enum),              &side,   VALUE,
        sizeof(MORSE_enum),              &trans,  VALUE,
        sizeof(int),                     &m1,     VALUE,
        sizeof(int),                     &n1,     VALUE,
        sizeof(int),                     &m2,     VALUE,
        sizeof(int),                     &n2,     VALUE,
        sizeof(int),                     &k,      VALUE,
        sizeof(int),                     &ib,     VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(A1, MORSE_Complex64_t, A1m, A1n), INOUT | REGION_L | REGION_D,
        sizeof(int),                     &lda1,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(A2, MORSE_Complex64_t, A2m, A2n), INOUT | REGION_FULL,
        sizeof(int),                     &lda2,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,  RTBLKADDR(V,  MORSE_Complex64_t, Vm,  Vn),  INPUT | REGION_FULL,
        sizeof(int),                     &ldv,    VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,  RTBLKADDR(T,  MORSE_Complex64_t, Tm,  Tn),  INPUT | REGION_FULL,
        sizeof(int),                     &ldt,    VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,  NULL,   SCRATCH,
        sizeof(int),                     &ldwork, VALUE,
        0);
}
