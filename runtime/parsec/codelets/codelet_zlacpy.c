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
 * @file codelet_zlacpy.c
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
static inline int
CORE_zlacpyx_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_enum *uplo;
    int *M;
    int *N;
    int *displA;
    MORSE_Complex64_t *A;
    int *LDA;
    int *displB;
    MORSE_Complex64_t *B;
    int *LDB;

    parsec_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &M,
        UNPACK_VALUE, &N,
        UNPACK_VALUE, &displA,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &LDA,
        UNPACK_VALUE, &displB,
        UNPACK_DATA,  &B,
        UNPACK_VALUE, &LDB );

    CORE_zlacpy(*uplo, *M, *N, A + (*displA), *LDA, B + (*displB), *LDB);

    (void)context;
    return 0;
}

void MORSE_TASK_zlacpyx( const MORSE_option_t *options,
                         MORSE_enum uplo, int m, int n, int nb,
                         int displA, const MORSE_desc_t *A, int Am, int An, int lda,
                         int displB, const MORSE_desc_t *B, int Bm, int Bn, int ldb )
{

    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlacpyx_parsec, options->priority, "lacpy",
        sizeof(MORSE_enum),    &uplo,                      VALUE,
        sizeof(int),           &m,                         VALUE,
        sizeof(int),           &n,                         VALUE,
        sizeof(int),           &displA,                    VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),   INPUT,
        sizeof(int),           &lda,                       VALUE,
        sizeof(int),           &displB,                    VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ),   OUTPUT,
        sizeof(int),           &ldb,                       VALUE,
        0);
}

void MORSE_TASK_zlacpy(const MORSE_option_t *options,
                       MORSE_enum uplo, int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    MORSE_TASK_zlacpyx( options, uplo, m, n, nb,
                        0, A, Am, An, lda,
                        0, B, Bm, Bn, ldb );
}
