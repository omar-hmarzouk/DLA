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
 * @file qwrapper_zherfb.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

static inline int
CORE_zherfb_parsec(dague_execution_unit_t    *context,
                   dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    int *n;
    int *k;
    int *ib;
    int *nb;
    MORSE_Complex64_t *A;
    int *lda;
    MORSE_Complex64_t *T;
    int *ldt;
    MORSE_Complex64_t *C;
    int *ldc;
    MORSE_Complex64_t *WORK;
    int *ldwork;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE,   &uplo,
                          UNPACK_VALUE,   &n,
                          UNPACK_VALUE,   &k,
                          UNPACK_VALUE,   &ib,
                          UNPACK_VALUE,   &nb,
                          UNPACK_DATA,    &A,
                          UNPACK_VALUE,   &lda,
                          UNPACK_DATA,    &T,
                          UNPACK_VALUE,   &ldt,
                          UNPACK_DATA,    &C,
                          UNPACK_VALUE,   &ldc,
                          UNPACK_SCRATCH, &WORK,
                          UNPACK_VALUE,   &ldwork);

    CORE_zherfb(*uplo, *n, *k, *ib, *nb,
                A, *lda, T, *ldt,
                C, *ldc, WORK, *ldwork);
}

void MORSE_TASK_zherfb(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zherfb_parsec, "herfb",
        sizeof(MORSE_enum),                &uplo, VALUE,
        sizeof(int),                       &n,    VALUE,
        sizeof(int),                       &k,    VALUE,
        sizeof(int),                       &ib,   VALUE,
        sizeof(int),                       &nb,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An), (uplo == MorseUpper) ? INOUT | REGION_U : INOUT | REGION_L,
        sizeof(int),                       &lda,  VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,    RTBLKADDR(T, MORSE_Complex64_t, Tm, Tn), INPUT | REGION_FULL,
        sizeof(int),                       &ldt,  VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn), (uplo == MorseUpper) ? INOUT | REGION_D | REGION_U : INOUT | REGION_D | REGION_L,
        sizeof(int),                       &ldc,  VALUE,
        sizeof(MORSE_Complex64_t)*2*nb*nb,  NULL, SCRATCH,
        sizeof(int),                       &nb,   VALUE,
        0);
}
