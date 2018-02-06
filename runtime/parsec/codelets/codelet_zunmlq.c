/**
 *
 * @file codelet_zunmlq.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 **/
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zunmlq_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_enum side;
    MORSE_enum trans;
    int m;
    int n;
    int k;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *C;
    int ldc;
    MORSE_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m, &n, &k, &ib, &A, &lda, &T, &ldt, &C, &ldc, &WORK, &ldwork );

    CORE_zunmlq( side, trans, m, n, k, ib,
                 A, lda, T, ldt, C, ldc, WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zunmlq(const MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int k, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       const MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zunmlq_parsec, options->priority, "unmlq",
        sizeof(MORSE_enum),                 &side,              VALUE,
        sizeof(MORSE_enum),                 &trans,             VALUE,
        sizeof(int),                        &m,                 VALUE,
        sizeof(int),                        &n,                 VALUE,
        sizeof(int),                        &k,                 VALUE,
        sizeof(int),                        &ib,                VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT,
        sizeof(int),                        &lda,               VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ),     INPUT,
        sizeof(int),                        &ldt,               VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, MORSE_Complex64_t, Cm, Cn ),     INOUT,
        sizeof(int),                        &ldc,               VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,               SCRATCH,
        sizeof(int),                        &nb,                VALUE,
        PARSEC_DTD_ARG_END );
}
