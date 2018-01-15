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
 * @file pztile.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "control/common.h"

#define A(m, n) A, m, n
#define B(m, n) &B, m, n

/*******************************************************************************
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void morse_pzlapack_to_tile(MORSE_Complex64_t *Af77, int ldaf77, MORSE_desc_t *A,
                            MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t B;
    int m, n;
    int ldam;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    B = morse_desc_init(
        MorseComplexDouble, A->mb, A->nb, A->bsiz,
        ldaf77, A->n, 0, 0, A->m, A->n, 1, 1);

    B.get_blkaddr = morse_getaddr_cm;
    B.get_blkldd  = morse_getblkldd_cm;
    B.mat = Af77;
    B.styp = MorseCM;

    RUNTIME_desc_create( &B );

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        for (n = 0; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            MORSE_TASK_zlacpy(
                &options,
                MorseUpperLower,
                tempmm, tempnn, A->mb,
                B(m, n), ldaf77,
                A(m, n), ldam);
        }
    }

    RUNTIME_desc_flush( &B, sequence );
    RUNTIME_sequence_wait( morse, sequence );
    RUNTIME_options_finalize( &options, morse );
    RUNTIME_desc_destroy( &B );
}

/*******************************************************************************
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void morse_pztile_to_lapack(MORSE_desc_t *A, MORSE_Complex64_t *Af77, int ldaf77,
                            MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t B;
    int m, n;
    int ldam;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    B = morse_desc_init(
        MorseComplexDouble, A->mb, A->nb, A->bsiz,
        ldaf77, A->n, 0, 0, A->m, A->n, 1, 1);

    B.get_blkaddr = morse_getaddr_cm;
    B.get_blkldd  = morse_getblkldd_cm;
    B.mat  = Af77;
    B.styp = MorseCM;

    RUNTIME_desc_create( &B );

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        for (n = 0; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            MORSE_TASK_zlacpy(
                &options,
                MorseUpperLower,
                tempmm, tempnn, A->mb,
                A(m, n), ldam,
                B(m, n), ldaf77);
        }
    }

    RUNTIME_desc_flush( &B, sequence );
    RUNTIME_sequence_wait( morse, sequence );
    RUNTIME_options_finalize( &options, morse );
    RUNTIME_desc_destroy( &B );
}
