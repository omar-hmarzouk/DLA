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
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
static int
CORE_zpotrf_parsec(parsec_execution_stream_t *context, parsec_task_t *this_task)
{
    MORSE_enum *uplo;
    int *tempkm, *ldak, *iinfo;
    MORSE_Complex64_t *A;

    parsec_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &tempkm,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &ldak,
        UNPACK_VALUE, &iinfo );

    CORE_zpotrf(*uplo, *tempkm, A, *ldak, iinfo);

    return 0;
}

void MORSE_TASK_zpotrf(const MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zpotrf_parsec, options->priority, "potrf",
        sizeof(MORSE_enum),    &uplo,                             VALUE,
        sizeof(int),           &n,                                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | AFFINITY,
        sizeof(int),           &lda,                              VALUE,
        sizeof(int),           &iinfo,                            VALUE,
        0);
}
