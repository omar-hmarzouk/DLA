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
 * @file codelet_zaxpy.c
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

void MORSE_TASK_zaxpy(MORSE_option_t *options,
                      int M, MORSE_Complex64_t *alpha,
                      MORSE_desc_t *A, int Am, int An, int incA,
                      MORSE_desc_t *B, int Bm, int Bn, int incB)
{
	quark_option_t *opt = (quark_option_t*)(options->schedopt);
	DAG_CORE_AXPY;
    QUARK_Insert_Task(opt->quark, CORE_zaxpy_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &M,         VALUE,
        sizeof(MORSE_Complex64_t),          alpha,      VALUE,
        sizeof(MORSE_Complex64_t)*M,        RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT,
        sizeof(int),                        &incA,      VALUE,
        sizeof(MORSE_Complex64_t)*M,        RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn), INOUT,
        sizeof(int),                        &incB,      VALUE,
        0);
}


void CORE_zaxpy_quark(Quark *quark)
{
    int M;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int incA;
    MORSE_Complex64_t *B;
    int incB;

    quark_unpack_args_6(quark, M, alpha, A, incA, B, incB);
    cblas_zaxpy(M, CBLAS_SADDR(alpha), A, incA, B, incB);
}
