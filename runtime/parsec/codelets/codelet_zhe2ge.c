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
 * @file codelet_zhe2ge.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @precisions normal z -> c d s
 *
 **/

#include "runtime/parsec/include/morse_parsec.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
static inline int
CORE_zhe2ge_parsec(dague_execution_unit_t    *context,
                   dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    int *M;
    int *N;
    const MORSE_Complex64_t *A;
    int *LDA;
    MORSE_Complex64_t *B;
    int *LDB;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &uplo,
                          UNPACK_VALUE, &M,
                          UNPACK_VALUE, &N,
                          UNPACK_DATA,  &A,
                          UNPACK_VALUE, &LDA,
                          UNPACK_DATA,  &B,
                          UNPACK_VALUE, &LDB);

    CORE_zhe2ge(*uplo, *M, *N, A, *LDA, B, *LDB);
}


void MORSE_TASK_zhe2ge(const MORSE_option_t *options,
                       MORSE_enum uplo,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zhe2ge_parsec, "he2ge",
        sizeof(MORSE_enum),              &uplo,   VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(MORSE_Complex64_t)*mb*mb,  RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT  | REGION_FULL,
        sizeof(int),                     &lda,    VALUE,
        sizeof(MORSE_Complex64_t)*mb*mb,  RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn), OUTPUT | REGION_FULL,
        sizeof(int),                     &ldb,    VALUE,
        0);
}
