/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014, 2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_zlascal.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
static inline int
CORE_zlascal_parsec(dague_execution_unit_t    *context,
                    dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    int *M;
    int *N;
    MORSE_Complex64_t *alpha;
    MORSE_Complex64_t *A;
    int *LDA;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &M,
        UNPACK_VALUE, &N,
        UNPACK_VALUE, &alpha,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &LDA);

    CORE_zlascal(*uplo, *M, *N, *alpha, A, *LDA);
}

void MORSE_TASK_zlascal(const MORSE_option_t *options,
                        MORSE_enum uplo,
                        int m, int n, int nb,
                        MORSE_Complex64_t alpha,
                        const MORSE_desc_t *A, int Am, int An, int lda)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zlascal_parsec, "lascal",
        sizeof(MORSE_enum),        &uplo,  VALUE,
        sizeof(int),               &m,     VALUE,
        sizeof(int),               &n,     VALUE,
        sizeof(MORSE_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF,              RTBLKADDR(A, MORSE_Complex64_t, Am, An), INOUT | REGION_FULL,
        sizeof(int),               &lda,   VALUE,
        0);
}


