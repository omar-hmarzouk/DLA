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
 * @file pzpotrf.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
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
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void morse_pzpotrf(MORSE_enum uplo, MORSE_desc_t *A,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldak, ldam, ldan;
    int tempkm, tempmm, tempnn;
    size_t ws_host   = 0;

    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

#ifdef CHAMELEON_USE_MAGMA
    if (0) /* Disable the workspace as long as it is is not used*/
    {
        int nb = magma_get_zpotrf_nb(A->nb);
        ws_host = sizeof(MORSE_Complex64_t)*nb*nb;
    }
#endif
    RUNTIME_options_ws_alloc( &options, 0, ws_host );

    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);

            MORSE_TASK_zpotrf(
                &options,
                MorseLower, tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, MorseLower, MorseConjTrans, MorseNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                ldan = BLKLDD(A, n);
                MORSE_TASK_zherk(
                    &options,
                    MorseLower, MorseNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k), ldan,
                     1.0, A(n, n), ldan);

                for (m = n+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k), ldam,
                               A(n, k), ldan,
                        zone,  A(m, n), ldam);
                }
                MORSE_TASK_dataflush( &options, A(n, k) );
            }
        }
    }
    /*
     *  MorseUpper
     */
    else {
        for (k = 0; k < A->nt; k++) {
            tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            MORSE_TASK_zpotrf(
                &options,
                MorseUpper,
                tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, MorseUpper, MorseConjTrans, MorseNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                MORSE_TASK_zherk(
                    &options,
                    MorseUpper, MorseConjTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam);

                for (n = m+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m), ldak,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam);
                }
                MORSE_TASK_dataflush( &options, A(k, m) );
            }
        }
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}
