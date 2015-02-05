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
 * @file pzplghe.c
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
 * @precisions normal z -> c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void morse_pzplghe( double bump, MORSE_desc_t *A, unsigned long long int seed,
                           MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n;
    int ldam;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

            MORSE_TASK_zplghe( 
                &options,
                bump, tempmm, tempnn, A(m, n), ldam,
                A->m, m*A->mb, n*A->nb, seed );
        
        }
    }
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
}
