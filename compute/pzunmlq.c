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
 * @file pzunmlq.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Azzam Haidar
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
#define T(m,n) T,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(k) DIAG, k, 0
#else
#define DIAG(k) A, k, k
#endif

/***************************************************************************//**
 *  Parallel application of Q using tile V - LQ factorization - dynamic scheduling
 **/
void morse_pzunmlq(MORSE_enum side, MORSE_enum trans,
                   MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int ldak, ldbk, ldbm;
    int tempmm, tempnn, tempkn, tempkm, tempkmin;
    int ib, minMT, minM;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    if (A->m > A->n) {
        minM  = A->n;
        minMT = A->nt;
    } else {
        minM  = A->m;
        minMT = A->mt;
    }

    /*
     * zunmlq = A->nb * ib
     * ztsmlq = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_MAGMA)
    /* Worker space
     *
     * zunmlq = A->nb * ib
     * ztsmlq = 2 * A->nb * ib
     */
    ws_worker = max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* necessary to avoid dependencies between tasks regarding the diag tile */
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, minMT*A->mb, A->nb, 0, 0, minMT*A->mb, A->nb);

    if (side == MorseLeft ) {
        if (trans == MorseNoTrans) {
            /*
             *  MorseLeft / MorseNoTrans
             */
            for (k = 0; k < minMT; k++) {
                tempkm   = k == B->mt -1 ? B->m -k*B->mb : B->mb;
                tempkmin = k == minMT-1 ? minM-k*A->nb : A->nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkmin, tempkm, A->nb,
                    A(k, k), ldak,
                    DIAG(k), A->mb );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkmin, tempkm,
                    0., 1.,
                    DIAG(k), A->mb );
#endif
                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    MORSE_TASK_zunmlq(
                        &options,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T->nb,
                        DIAG(k), A->mb,
                        T(k, k), T->mb,
                        B(k, n), ldbk);
                }
                for (m = k+1; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztsmlq(
                            &options,
                            side, trans,
                            B->mb, tempnn, tempmm, tempnn, tempkmin, ib, T->nb,
                            B(k, n), ldbk,
                            B(m, n), ldbm,
                            A(k, m), ldak,
                            T(k, m), T->mb);
                    }
                }
            }
        }
        else {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = minMT-1; k >= 0; k--) {
                tempkm   = k == B->mt -1 ? B->m -k*B->mb : B->mb;
                tempkmin = k == minMT-1 ? minM-k*A->nb : A->nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (m = B->mt-1; m > k; m--) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B->nt; n++) {
                        tempnn   = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztsmlq(
                            &options,
                            side, trans,
                            B->mb, tempnn, tempmm, tempnn, tempkmin, ib, T->nb,
                            B(k, n), ldbk,
                            B(m, n), ldbm,
                            A(k, m), ldak,
                            T(k, m), T->mb);
                    }
                }
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkmin, tempkm, A->nb,
                    A(k, k), ldak,
                    DIAG(k), A->mb );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkmin, tempkm,
                    0., 1.,
                    DIAG(k), A->mb );
#endif
                for (n = 0; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    MORSE_TASK_zunmlq(
                        &options,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T->nb,
                        DIAG(k), A->mb,
                        T(k, k), T->mb,
                        B(k, n), ldbk);
                }
            }
        }
    }
    else {
        if (trans == MorseNoTrans) {
            /*
             *  MorseRight / MorseNoTrans
             */
            for (k = minMT-1; k >= 0; k--) {
                tempkn   = k == B->nt -1 ? B->n -k*B->nb : B->nb;
                tempkmin = k == minMT-1 ? minM-k*A->nb : A->nb;
                ldak = BLKLDD(A, k);
                for (n = B->nt-1; n > k; n--) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztsmlq(
                            &options,
                            side, trans,
                            tempmm, B->nb, tempmm, tempnn, tempkmin, ib, T->nb,
                            B(m, k), ldbm,
                            B(m, n), ldbm,
                            A(k, n), ldak,
                            T(k, n), T->mb);
                    }
                }
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkmin, tempkn, A->nb,
                    A(k, k), ldak,
                    DIAG(k), A->mb );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkmin, tempkn,
                    0., 1.,
                    DIAG(k), A->mb );
#endif
                for (m = 0; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    MORSE_TASK_zunmlq(
                        &options,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T->nb,
                        DIAG(k), A->mb,
                        T(k, k), T->mb,
                        B(m, k), ldbm);
                }
            }
        }
        else {
            /*
             *  MorseRight / MorseConjTrans
             */
            for (k = 0; k < minMT; k++) {
                tempkn   = k == B->nt -1 ? B->n -k*B->nb : B->nb;
                tempkmin = k == minMT-1 ? minM-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkmin, tempkn, A->nb,
                    A(k, k), ldak,
                    DIAG(k), A->mb );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkmin, tempkn,
                    0., 1.,
                    DIAG(k), A->mb );
#endif
                for (m = 0; m < B->mt; m++) {
                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    MORSE_TASK_zunmlq(
                        &options,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T->nb,
                        DIAG(k), A->mb,
                        T(k, k), T->mb,
                        B(m, k), ldbm);
                }
                for (n = k+1; n < B->nt; n++) {
                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_ztsmlq(
                            &options,
                            side, trans,
                            tempmm, B->nb, tempmm, tempnn, tempkmin, ib, T->nb,
                            B(m, k), ldbm,
                            B(m, n), ldbm,
                            A(k, n), ldak,
                            T(k, n), T->mb);
                    }
                }
            }
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();

    morse_desc_mat_free(DIAG);
    free(DIAG);
}
