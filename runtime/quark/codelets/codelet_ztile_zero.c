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
 * @file codelet_ztile_zero.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/quark/include/morse_quark.h"

/*****************************************************************************
 *
 **/
void MORSE_TASK_ztile_zero(MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           MORSE_desc_t *A, int Am, int An, int lda)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(opt->quark, CORE_ztile_zero_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                       &X1,                                       VALUE,
        sizeof(int),                       &X2,                                       VALUE,
        sizeof(int),                       &Y1,                                       VALUE,
        sizeof(int),                       &Y2,                                       VALUE,
        sizeof(MORSE_Complex64_t)*A->bsiz,  RTBLKADDR(A, MORSE_Complex64_t, Am, An),  OUTPUT | LOCALITY,
        sizeof(int),                       &lda,                                      VALUE,
        0);
}

/*****************************************************************************
 *
 **/
void CORE_ztile_zero_quark(Quark *quark)
{
    int X1;
    int X2;
    int Y1;
    int Y2;
    MORSE_Complex64_t *A;
    int lda;

    int x, y;

    quark_unpack_args_6(quark, X1, X2, Y1, Y2, A, lda);

    for (x = X1; x < X2; x++)
        for (y = Y1; y < Y2; y++)
            A[lda*x+y] = 0.0;

}
