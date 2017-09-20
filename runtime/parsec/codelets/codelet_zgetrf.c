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
 * @file codelet_zgetrf.c
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
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static int
CORE_zgetrf_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *m;
    int *n;
    MORSE_Complex64_t *A;
    int *lda;
    int *IPIV;
    MORSE_bool *check_info;
    int *iinfo;
    int info;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &m,
        UNPACK_VALUE, &n,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &lda,
        UNPACK_SCRATCH, &IPIV,
        UNPACK_VALUE, &check_info,
        UNPACK_VALUE, &iinfo );

    CORE_zgetrf( *m, *n, A, *lda, IPIV, &info );

    return 0;
}

void MORSE_TASK_zgetrf(const MORSE_option_t *options,
                       int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zgetrf_parsec, "getrf",
        sizeof(int),        &m,                          VALUE,
        sizeof(int),        &n,                          VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | REGION_FULL,
        sizeof(int),        &lda,                        VALUE,
        sizeof(int)*nb,      IPIV,                        SCRATCH,
        sizeof(MORSE_bool), &check_info,                 VALUE,
        sizeof(int),        &iinfo,                      VALUE,
        0);
}
