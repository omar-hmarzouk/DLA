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
 * @file qwrapper_zlatro.c
 *
 *  PLASMA core_blas quark wrapper
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlatro_parsec(dague_execution_unit_t    *context,
                   dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    MORSE_enum *trans;
    int *M;
    int *N;
    const MORSE_Complex64_t *A;
    int *LDA;
    MORSE_Complex64_t *B;
    int *LDB;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &trans,
        UNPACK_VALUE, &M,
        UNPACK_VALUE, &N,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &LDA,
        UNPACK_DATA,  &B,
        UNPACK_VALUE, &LDB);

    CORE_zlatro(*uplo, *trans, *M, *N,
                A, *LDA, B, *LDB);
}

/***************************************************************************//**
 *
 **/
void MORSE_TASK_zlatro(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int m, int n, int mb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zlatro_parsec, "latro",
        sizeof(MORSE_enum), &uplo,  VALUE,
        sizeof(MORSE_enum), &trans, VALUE,
        sizeof(int),        &m,     VALUE,
        sizeof(int),        &n,     VALUE,
        PASSED_BY_REF,       RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT  | REGION_FULL,
        sizeof(int),        &lda,   VALUE,
        PASSED_BY_REF,       RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn), OUTPUT | REGION_FULL,
        sizeof(int),        &ldb,   VALUE,
        0);
}
