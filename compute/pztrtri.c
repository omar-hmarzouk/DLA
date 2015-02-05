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
 * @file pztrtri.c
 *
 *  MORSE auxiliary routines
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
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
/***************************************************************************//**
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 **/
void morse_pztrtri(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldam, ldan;
    int tempkn, tempmm, tempnn;

    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);
    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (n = 0; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            ldan = BLKLDD(A, n);
            for (m = n+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, diag,
                    tempmm, tempnn, A->mb,
                    mzone, A(n, n), ldan,
                           A(m, n), ldam);
            }
            for (m = n+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                for (k = 0; k < n; k++) {
                    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        tempmm, tempkn, tempnn, A->mb,
                        zone, A(m, n), ldam,
                              A(n, k), ldan,
                        zone, A(m, k), ldam);
                }
            }
            for (m = 0; m < n; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, diag,
                    tempnn, tempmm, A->mb,
                    zone, A(n, n), ldan,
                          A(n, m), ldan);
            }
            MORSE_TASK_ztrtri(
                &options,
                uplo, diag,
                tempnn, A->mb,
                A(n, n), ldan, A->nb*n);
        }
    }
    /*
     *  MorseUpper
     */
    else {
        for (m = 0; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            for (n = m+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, diag,
                    tempmm, tempnn, A->mb,
                    mzone, A(m, m), ldam,
                           A(m, n), ldam);
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                ldan = BLKLDD(A, n);
                for (k = m+1; k < A->nt; k++) {
                    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        tempnn, tempkn, tempmm, A->mb,
                        zone, A(n, m), ldan,
                              A(m, k), ldam,
                        zone, A(n, k), ldan);
                }
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, diag,
                    tempnn, tempmm, A->mb,
                    zone, A(m, m), ldam,
                          A(n, m), ldan);
            }
            MORSE_TASK_ztrtri(
                &options,
                uplo, diag,
                tempmm, A->mb,
                A(m, m), ldam, A->mb*m);
        }
    }
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
}
