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
 * @file codelet_ztrasm.c
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
CORE_ztrasm_parsec(dague_execution_unit_t *context, dague_execution_context_t * this_task)
{
    MORSE_enum *storev;
    MORSE_enum *uplo;
    MORSE_enum *diag;
    int *M;
    int *N;
    dague_data_copy_t *gA;
    int *lda;
    dague_data_copy_t *gwork;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &storev,
                          UNPACK_VALUE, &uplo,
                          UNPACK_VALUE, &diag,
                          UNPACK_VALUE, &M,
                          UNPACK_VALUE, &N,
                          UNPACK_DATA,  &gA,
                          UNPACK_VALUE, &lda,
                          UNPACK_DATA,  &gwork
                        );


    void *A = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gA);
    void *work = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gwork);

    CORE_ztrasm(*storev, *uplo, *diag, *M, *N, A, *lda, work);

    return 0;
}

void MORSE_TASK_ztrasm(MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_ztrasm_parsec,    "trasm",
                            sizeof(MORSE_enum),     &storev,                VALUE,
                            sizeof(MORSE_enum),     &uplo,                  VALUE,
                            sizeof(MORSE_enum),     &diag,                  VALUE,
                            sizeof(int),            &M,                     VALUE,
                            sizeof(int),            &N,                     VALUE,
                            PASSED_BY_REF,          RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT | REGION_FULL,
                            sizeof(int),            &lda,                   VALUE,
                            PASSED_BY_REF,          RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ),     INOUT | REGION_FULL,
                             0);
}
