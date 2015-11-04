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
 * @file codelet_zpotrf.c
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

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
static int
CORE_zpotrf_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    int *tempkm, *ldak, *iinfo;
    dague_data_copy_t *data;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &uplo,
                          UNPACK_VALUE, &tempkm,
                          UNPACK_DATA,  &data,
                          UNPACK_VALUE, &ldak,
                          UNPACK_VALUE, &iinfo
                        );

    void *TT = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)data);

    CORE_zpotrf(*uplo, *tempkm, TT, *ldak, iinfo);

    return 0;
}

void MORSE_TASK_zpotrf(MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_zpotrf_parsec,               "potrf",
                             sizeof(MORSE_enum),    &uplo,                             VALUE,
                             sizeof(int),           &n,                                VALUE,
                             PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | REGION_FULL,
                             sizeof(int),           &lda,                              VALUE,
                             sizeof(int),           &iinfo,                            VALUE,
                             0);
}
