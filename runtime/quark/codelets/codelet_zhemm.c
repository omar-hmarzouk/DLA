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
 * @file codelet_zhemm.c
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
 * @precisions normal z -> c
 *
 **/
#include "runtime/quark/include/morse_quark.h"

#undef REAL
#define COMPLEX
#ifdef COMPLEX
/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_zhemm(MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_HEMM;
    QUARK_Insert_Task(opt->quark, CORE_zhemm_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &side,    VALUE,
        sizeof(MORSE_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(MORSE_Complex64_t),         &alpha,   VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(MORSE_Complex64_t),         &beta,    VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn),               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}


void CORE_zhemm_quark(Quark *quark)
{
    MORSE_enum side;
    MORSE_enum uplo;
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    cblas_zhemm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}
#endif
