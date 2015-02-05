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
 * @file codelet_zgessq.c
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

void MORSE_TASK_zgessq( MORSE_option_t *options,
                        int m, int n,
                        MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(opt->quark, CORE_zgessq_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                     &m,    VALUE,
        sizeof(int),                     &n,    VALUE,
        sizeof(MORSE_Complex64_t)*lda*n, RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &lda,  VALUE,
        sizeof(double)*2,                RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn), INOUT,
        0);
}


void CORE_zgessq_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    quark_unpack_args_5( quark, m, n, A, lda, SCALESUMSQ );
    CORE_zgessq( m, n, A, lda, &SCALESUMSQ[0], &SCALESUMSQ[1]);
}
