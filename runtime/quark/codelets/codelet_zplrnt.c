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
 * @file codelet_zplrnt.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/quark/include/morse_quark.h"



void MORSE_TASK_zplrnt( MORSE_option_t *options,
                        int m, int n, MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_PLRNT;
    QUARK_Insert_Task(opt->quark, CORE_zplrnt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(MORSE_Complex64_t)*lda*n, RTBLKADDR(A, MORSE_Complex64_t, Am, An),         OUTPUT,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}


void CORE_zplrnt_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_8( quark, m, n, A, lda, bigM, m0, n0, seed );
    CORE_zplrnt( m, n, A, lda, bigM, m0, n0, seed );
}

