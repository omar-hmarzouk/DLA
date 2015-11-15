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
 * @file codelet_zasum.c
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
CORE_dzasum_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    MORSE_enum *storev;
    MORSE_enum *uplo;
    int *M;
    int *N;
    MORSE_Complex64_t *A;
    int *lda;
    double *work;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &storev,
                          UNPACK_VALUE, &uplo,
                          UNPACK_VALUE, &M,
                          UNPACK_VALUE, &N,
                          UNPACK_DATA,  &A,
                          UNPACK_VALUE, &lda,
                          UNPACK_DATA,  &work
                        );

    CORE_dzasum(*storev, *uplo, *M, *N, A, *lda, work);

    return 0;
}

void MORSE_TASK_dzasum(MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, int M, int N,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_dzasum_parsec,               "zasum",
                             sizeof(MORSE_enum),    &storev,                           VALUE,
                             sizeof(MORSE_enum),    &uplo,                             VALUE,
                             sizeof(int),           &M,                                VALUE,
                             sizeof(int),           &N,                                VALUE,
                             PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT | REGION_FULL,
                             sizeof(int),           &lda,                              VALUE,
                             PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ),     INOUT | REGION_FULL,
                             0);
}
