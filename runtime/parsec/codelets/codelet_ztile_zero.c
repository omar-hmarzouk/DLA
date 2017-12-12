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
 * @file codelet_ztile_zero.c
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
CORE_ztile_zero_parsec(parsec_execution_stream_t *context, parsec_task_t *this_task)
{
    int *X1;
    int *X2;
    int *Y1;
    int *Y2;
    MORSE_Complex64_t *A;
    int *lda;
    int x, y;

    parsec_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &X1,
        UNPACK_VALUE, &X2,
        UNPACK_VALUE, &Y1,
        UNPACK_VALUE, &Y2,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &lda );

    for (x = *X1; x < *X2; x++)
        for (y = *Y1; y < *Y2; y++)
            A[*lda*x+y] = 0.0;

    return 0;
}

void MORSE_TASK_ztile_zero(const const MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           const MORSE_desc_t *A, int Am, int An, int lda)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztile_zero_parsec, options->priority, "tile zero",
        sizeof(int),       &X1,                       VALUE,
        sizeof(int),       &X2,                       VALUE,
        sizeof(int),       &Y1,                       VALUE,
        sizeof(int),       &Y2,                       VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     OUTPUT | REGION_FULL,
        sizeof(int),       &lda,                      VALUE,
        0);
}
