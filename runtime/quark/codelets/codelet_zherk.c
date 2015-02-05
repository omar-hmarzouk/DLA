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
 * @file codelet_zherk.c
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

void MORSE_TASK_zherk(MORSE_option_t *options,
                      MORSE_enum uplo, MORSE_enum trans,
                      int n, int k, int nb,
                      double alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      double beta, MORSE_desc_t *C, int Cm, int Cn, int ldc)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_HERK;
    QUARK_Insert_Task(opt->quark, CORE_zherk_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &uplo,      VALUE,
        sizeof(MORSE_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),                     &alpha,     VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double),                     &beta,      VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(C, MORSE_Complex64_t, Cm, Cn),                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}


void CORE_zherk_quark(Quark *quark)
{
    MORSE_enum uplo;
    MORSE_enum trans;
    int n;
    int k;
    double alpha;
    MORSE_Complex64_t *A;
    int lda;
    double beta;
    MORSE_Complex64_t *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_zherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc);
}
#endif
