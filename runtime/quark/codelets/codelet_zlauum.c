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
 * @file codelet_zlauum.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "coreblas/include/lapacke.h"
#include "runtime/quark/include/morse_quark.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_zlauum(MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LAUUM;
    QUARK_Insert_Task(opt->quark, CORE_zlauum_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &uplo,  VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),             INOUT,
        sizeof(int),                        &lda,   VALUE,
        0);
}


void CORE_zlauum_quark(Quark *quark)
{
    MORSE_enum uplo;
    int N;
    MORSE_Complex64_t *A;
    int LDA;

    quark_unpack_args_4(quark, uplo, N, A, LDA);
    LAPACKE_zlauum_work(LAPACK_COL_MAJOR, morse_lapack_const(uplo), N, A, LDA);
}
