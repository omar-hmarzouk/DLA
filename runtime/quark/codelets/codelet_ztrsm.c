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
 * @file codelet_ztrsm.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/quark/include/morse_quark.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_ztrsm(MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TRSM;
    QUARK_Insert_Task(opt->quark, CORE_ztrsm_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &side,      VALUE,
        sizeof(MORSE_enum),                &uplo,      VALUE,
        sizeof(MORSE_enum),                &transA,    VALUE,
        sizeof(MORSE_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(MORSE_Complex64_t),         &alpha,     VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}


void CORE_ztrsm_quark(Quark *quark)
{
    MORSE_enum side;
    MORSE_enum uplo;
    MORSE_enum transA;
    MORSE_enum diag;
    int m;
    int n;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    cblas_ztrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m, n,
        CBLAS_SADDR(alpha), A, lda,
        B, ldb);
}
