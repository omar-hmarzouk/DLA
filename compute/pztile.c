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

#define A(m, n) dA, m, n
#define B(m, n) &dB, m, n

/*******************************************************************************
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void morse_pzlapack_to_tile(MORSE_Complex64_t *Af77, int lda, MORSE_desc_t *dA,
                            MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t dB;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    dB = morse_desc_init(
        MorseComplexDouble, dA->mb, dA->nb, dA->bsiz,
        lda, dA->n, dA->i, dA->j, dA->m, dA->n, 1, 1);

    dB.get_blkaddr = morse_getaddr_cm;
    dB.get_blkldd  = morse_getblkldd_cm;
    dB.mat = Af77;
    dB.styp = MorseCM;

    RUNTIME_desc_create( &dB );

    for (m = 0; m < dA->mt; m++)
    {
        ldt = BLKLDD(dA, m);
        for (n = 0; n < dA->nt; n++)
        {
            X1 = n == 0 ? dA->j%dA->nb : 0;
            Y1 = m == 0 ? dA->i%dA->mb : 0;
            X2 = n == dA->nt-1 ? (dA->j+dA->n-1)%dA->nb+1 : dA->nb;
            Y2 = m == dA->mt-1 ? (dA->i+dA->m-1)%dA->mb+1 : dA->mb;

            MORSE_TASK_zlacpy(
                &options,
                MorseUpperLower,
                (Y2-Y1), (X2-X1), dA->mb,
                B(m, n), lda,
                A(m, n), ldt);
        }
    }

    RUNTIME_barrier( morse );
    RUNTIME_options_finalize( &options, morse );
    MORSE_TASK_dataflush_all();
    RUNTIME_desc_getoncpu( &dB );
    RUNTIME_desc_destroy( &dB );
}

/*******************************************************************************
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void morse_pztile_to_lapack(MORSE_desc_t *dA, MORSE_Complex64_t *Af77, int lda,
                            MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    MORSE_desc_t dB;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    dB = morse_desc_init(
        MorseComplexDouble, dA->mb, dA->nb, dA->bsiz,
        lda, dA->n, dA->i, dA->j, dA->m, dA->n, 1, 1);

    dB.get_blkaddr = morse_getaddr_cm;
    dB.get_blkldd  = morse_getblkldd_cm;
    dB.mat  = Af77;
    dB.styp = MorseCM;

    RUNTIME_desc_create( &dB );

    for (m = 0; m < dA->mt; m++)
    {
        Y1 = m == 0 ? dA->i%dA->mb : 0;
        Y2 = m == dA->mt-1 ? (dA->i+dA->m-1)%dA->mb+1 : dA->mb;
        ldt = BLKLDD(dA, m);

        for (n = 0; n < dA->nt; n++)
        {
            X1 = n == 0 ? dA->j%dA->nb : 0;
            X2 = n == dA->nt-1 ? (dA->j+dA->n-1)%dA->nb+1 : dA->nb;

            MORSE_TASK_zlacpy(
                &options,
                MorseUpperLower,
                (Y2-Y1), (X2-X1), dA->mb,
                A(m, n), ldt,
                B(m, n), lda);
        }
    }

    RUNTIME_barrier( morse );
    RUNTIME_options_finalize( &options, morse );
    MORSE_TASK_dataflush_all();
    RUNTIME_desc_getoncpu( &dB );
    RUNTIME_desc_destroy( &dB );
}


/*******************************************************************************
 *  Zeroes a submatrix in tile layout - dynamic scheduling
 **/
void morse_pztile_zero(MORSE_desc_t *dA, MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < dA->mt; m++)
    {
        Y1 = m == 0 ? dA->i%dA->mb : 0;
        Y2 = m == dA->mt-1 ? (dA->i+dA->m-1)%dA->mb+1 : dA->mb;
        ldt = BLKLDD(dA, m);
        for (n = 0; n < dA->nt; n++)
        {
            X1 = n == 0 ? dA->j%dA->nb : 0;
            X2 = n == dA->nt-1 ? (dA->j+dA->n-1)%dA->nb+1 : dA->nb;

            MORSE_TASK_ztile_zero(
                &options,
                X1, X2, Y1, Y2,
                A(m, n), ldt);
        }
    }

    RUNTIME_barrier( morse );
    RUNTIME_options_finalize( &options, morse );
    MORSE_TASK_dataflush_all();
}
