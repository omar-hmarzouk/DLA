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
 * @file pzgeadd.c
 *
 *  MORSE auxiliary routines
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
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
/***************************************************************************//**
 *
 **/

/***************************************************************************//**
 *
 **/
void morse_pzgeadd(MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int X, Y;
    int m, n;
    int ldam, ldbm;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < A->mt; m++) {
        X = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(B, m);

        for (n = 0; n < A->nt; n++) {
            Y = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            MORSE_TASK_zgeadd(
                &options,
                X, Y,
                alpha, A(m, n), ldam,
                       B(m, n), ldbm);
        }
    }
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
}
