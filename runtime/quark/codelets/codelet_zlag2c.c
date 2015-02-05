/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_zlag2c.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#include "coreblas/include/lapacke.h"
#include "runtime/quark/include/morse_quark.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
void MORSE_TASK_zlag2c(MORSE_option_t *options,
                       int m, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAG2C;
    QUARK_Insert_Task(opt->quark, CORE_zlag2c_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(MORSE_Complex32_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex32_t, Bm, Bn),                 OUTPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(MORSE_sequence_t*),           &(options->sequence),  VALUE,
        sizeof(MORSE_request_t*),            &(options->request),   VALUE,
        0);
}



void CORE_zlag2c_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex32_t *B;
    int ldb;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    int info;

    quark_unpack_args_8(quark, m, n, A, lda, B, ldb, sequence, request);
    info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
    if (sequence->status == MORSE_SUCCESS && info != 0)
        RUNTIME_sequence_flush(quark, sequence, request, info);
}

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_clag2z(MORSE_option_t *options,
                       int m, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    QUARK_Insert_Task(opt->quark, CORE_clag2z_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(MORSE_Complex32_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex32_t, Am, An),             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}


void CORE_clag2z_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex32_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;

    quark_unpack_args_6(quark, m, n, A, lda, B, ldb);
    LAPACKE_clag2z_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

