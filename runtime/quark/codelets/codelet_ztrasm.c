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
 * @file codelet_ztrasm.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
 *  University of Bordeaux, Bordeaux INP
 *
 * @version 2.6.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/quark/include/morse_quark.h"

void MORSE_TASK_ztrasm(MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    int szeW = uplo == MorseRowwise ? M : N ;
    QUARK_Insert_Task(opt->quark, CORE_ztrasm_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &storev,    VALUE,
        sizeof(MORSE_enum),              &uplo,      VALUE,
        sizeof(MORSE_enum),              &diag,      VALUE,
        sizeof(int),                     &M,         VALUE,
        sizeof(int),                     &N,         VALUE,
        sizeof(MORSE_Complex64_t)*lda*N, RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INPUT,
        sizeof(int),                     &lda,       VALUE,
        sizeof(double)*szeW,             RTBLKADDR(B, double, Bm, Bn), INOUT,
        0);
}


void CORE_ztrasm_quark(Quark *quark)
{
    MORSE_enum storev;
    MORSE_enum uplo;
    MORSE_enum diag;
    int M;
    int N;
    MORSE_Complex64_t *A;
    int lda;
    double *work;

    quark_unpack_args_8(quark, storev, uplo, diag, M, N, A, lda, work);
    CORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
}
