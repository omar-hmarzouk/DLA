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
 * @file codelet_ztrsm.c
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
CORE_ztrsm_parsec(parsec_execution_stream_t *context, parsec_task_t *this_task)
{
    MORSE_enum *side, *uplo, *trans, *diag;
    int  *tempmm, *nb, *ldak, *ldam;
    MORSE_Complex64_t *alpha;
    MORSE_Complex64_t *T;
    MORSE_Complex64_t *C;

    parsec_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &side,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &trans,
        UNPACK_VALUE, &diag,
        UNPACK_VALUE, &tempmm,
        UNPACK_VALUE, &nb,
        UNPACK_VALUE, &alpha,
        UNPACK_DATA,  &T,
        UNPACK_VALUE, &ldak,
        UNPACK_DATA,  &C,
        UNPACK_VALUE, &ldam );

    CORE_ztrsm(*side, *uplo, *trans, *diag,
           *tempmm, *nb, *alpha, T, *ldak,
           C, *ldam);

    return 0;
}

void MORSE_TASK_ztrsm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ztrsm_parsec, options->priority, "Trsm",
        sizeof(MORSE_enum),    &side,                     VALUE,
        sizeof(MORSE_enum),    &uplo,                     VALUE,
        sizeof(MORSE_enum),    &transA,                   VALUE,
        sizeof(MORSE_enum),    &diag,                     VALUE,
        sizeof(int),           &m,                        VALUE,
        sizeof(int),           &n,                        VALUE,
        sizeof(MORSE_Complex64_t),           &alpha,      VALUE,
        PASSED_BY_REF,     RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT,
        sizeof(int),           &lda,                      VALUE,
        PASSED_BY_REF,     RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ),     INOUT | AFFINITY,
        sizeof(int),           &ldb,                      VALUE,
        0);
}
