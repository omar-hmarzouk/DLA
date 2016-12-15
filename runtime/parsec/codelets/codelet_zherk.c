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
 * @file codelet_zherk.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

static int
CORE_zherk_parsec(dague_execution_unit_t *context, dague_execution_context_t * this_task)
{
    MORSE_enum *uplo;
    MORSE_enum *trans;
    int *n;
    int *k;
    double *alpha;
    MORSE_Complex64_t *A;
    int *lda;
    double *beta;
    MORSE_Complex64_t *C;
    int *ldc;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &uplo,
                          UNPACK_VALUE, &trans,
                          UNPACK_VALUE, &n,
                          UNPACK_VALUE, &k,
                          UNPACK_VALUE, &alpha,
                          UNPACK_DATA,  &A,
                          UNPACK_VALUE, &lda,
                          UNPACK_VALUE, &beta,
                          UNPACK_DATA,  &C,
                          UNPACK_VALUE, &ldc
                        );


    CORE_zherk(*uplo, *trans, *n, *k,
               *alpha, A, *lda,
               *beta,  C, *ldc);

    return 0;
}

void MORSE_TASK_zherk(const MORSE_option_t *options,
                      MORSE_enum uplo, MORSE_enum trans,
                      int n, int k, int nb,
                      double alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      double beta, const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(DAGUE_dtd_handle,      CORE_zherk_parsec,                "herk",
                             sizeof(MORSE_enum),    &uplo,                             VALUE,
                             sizeof(MORSE_enum),    &trans,                            VALUE,
                             sizeof(int),           &n,                                VALUE,
                             sizeof(int),           &k,                                VALUE,
                             sizeof(double),        &alpha,                            VALUE,
                             PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT | REGION_FULL,
                             sizeof(int),           &lda,                              VALUE,
                             sizeof(double),        &beta,                             VALUE,
                             PASSED_BY_REF,         RTBLKADDR( C, MORSE_Complex64_t, Cm, Cn ),     INOUT | REGION_FULL,
                             sizeof(int),           &ldc,                              VALUE,
                             0);
}

