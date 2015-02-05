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
 * @file pzunmqrrh.c
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
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  (m),  (n)
#define B(m,n) B,  (m),  (n)
#define T(m,n) T,  (m),  (n)
#define T2(m,n) T,  (m),  ((n)+A->nt)
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(m,n) DIAG, ((m)/BS), 0
#else
#define DIAG(m,n) A, (m), (n)
#endif

/***************************************************************************//**
 *  Parallel application of Q using tile V - QR factorization (reduction
 *  Householder) - dynamic scheduling
 **/
void morse_pzunmqrrh(MORSE_enum side, MORSE_enum trans,
                     MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, int BS,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int K, M, RD, lastRD;
    int ldaM, ldam, ldan, ldaMRD;
    int ldbM, ldbm, ldbMRD;
    int tempMm, tempkn, tempnn, tempmm, tempMRDm, tempkmin;
    int ib;
    int nblk;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zunmqr = A->nb * ib
     * ztsmqr = A->nb * ib
     * zttmqr = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_MAGMA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* necessary to avoid dependencies between tasks regarding the diag tile */
    nblk = ( A->mt + BS -1 ) / BS;
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, nblk * A->mb, A->nb, 0, 0, nblk * A->mb, A->nb);

    K = min(A->mt, A->nt);
    if (side == MorseLeft ) {
        if (trans == MorseConjTrans) {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = 0; k < K; k++) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    ldbM = BLKLDD(B, M);
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempMm, tempkmin, A->nb,
                        A(M, k), ldaM,
                        DIAG(M, k), ldaM );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempMm, tempkmin,
                        0., 1.,
                        DIAG(M, k), ldaM );
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempMm, tempnn,
                            tempkmin, ib, T->nb,
                            DIAG(M, k), ldaM,
                            T(M, k), T->mb,
                            B(M, n), ldbM);
                    }
                    for (m = M+1; m < min(M+BS, A->mt); m++) {
                        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                        ldbm = BLKLDD(B, m);
                        ldam = BLKLDD(A, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                A->nb, tempnn, tempmm, tempnn,
                                tempkn, ib, T->nb,
                                B(M, n), ldbM,
                                B(m, n), ldbm,
                                A(m, k), ldam,
                                T(m, k), T->mb);
                        }
                    }
                }
                for (RD = BS; RD < A->mt-k; RD *= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldbM   = BLKLDD(B, M   );
                        ldbMRD = BLKLDD(B, M+RD);
                        ldaMRD = BLKLDD(A, M+RD);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                A->nb, tempnn, tempMRDm, tempnn,
                                tempkn, ib, T->nb,
                                B (M,    n), ldbM,
                                B (M+RD, n), ldbMRD,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb);
                        }
                    }
                }
            }
        } else {
            /*
             *  MorseLeft / MorseNoTrans
             */
            for (k = K-1; k >= 0; k--) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                lastRD = 0;
                for (RD = BS; RD < A->mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldbM   = BLKLDD(B, M   );
                        ldbMRD = BLKLDD(B, M+RD);
                        ldaMRD = BLKLDD(A, M+RD);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                A->nb, tempnn, tempMRDm, tempnn,
                                tempkn, ib, T->nb,
                                B (M,    n), ldbM,
                                B (M+RD, n), ldbMRD,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb);
                        }
                    }
                }
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    ldbM = BLKLDD(B, M);
                    for (m = min(M+BS, A->mt)-1; m > M; m--) {
                        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                        ldbm = BLKLDD(B, m);
                        ldam = BLKLDD(A, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                A->nb, tempnn, tempmm, tempnn,
                                tempkn, ib, T->nb,
                                B(M, n), ldbM,
                                B(m, n), ldbm,
                                A(m, k), ldam,
                                T(m, k), T->mb);
                        }
                    }
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempMm, tempkmin, A->nb,
                        A(M, k), ldaM,
                        DIAG(M, k), ldaM );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempMm, tempkmin,
                        0., 1.,
                        DIAG(M, k), ldaM );
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempMm, tempnn,
                            tempkmin, ib, T->nb,
                            DIAG(M, k), ldaM,
                            T(M, k), T->mb,
                            B(M, n), ldbM);
                    }
                }
            }
        }
    } else {
        if (trans == MorseConjTrans) {
            /*
             *  MorseRight / MorseConjTrans
             */
            for (k = K-1; k >= 0; k--) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                lastRD = 0;
                for (RD = BS; RD < A->mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldaMRD = BLKLDD(A, M+RD);
                        for (m = 0; m < B->mt; m++) {
                            ldbm   = BLKLDD(B, m);
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempMRDm,
                                tempkn, ib, T->nb,
                                B (m, M), ldbm,
                                B (m, M+RD), ldbm,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb);
                        }
                    }
                }
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    ldbM = BLKLDD(B, M);
                    for (n = min(M+BS, A->mt)-1; n > M; n--) {
                        ldan = BLKLDD(A, n);
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        for (m = 0; m < B->mt; m++) {
                            ldbm = BLKLDD(B, m);
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                tempmm, tempMm, tempmm, tempnn,
                                tempkn, ib, T->nb,
                                B(m, M), ldbm,
                                B(m, n), ldbm,
                                A(n, k), ldan,
                                T(n, k), T->mb);
                        }
                    }
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempMm, tempkmin, A->nb,
                        A(M, k), ldaM,
                        DIAG(M, k), ldaM );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempMm, tempkmin,
                        0., 1.,
                        DIAG(M, k), ldaM );
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempMm,
                            tempkmin, ib, T->nb,
                            DIAG(M, k), ldaM,
                            T(M, k), T->mb,
                            B(m, M), ldbm);
                    }
                }
            }
        } else {
            /*
             *  MorseRight / MorseNoTrans
             */
            for (k = 0; k < K; k++) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                for (M = k; M < A->mt; M += BS) {
                    tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempMm, tempkmin, A->nb,
                        A(M, k), ldaM,
                        DIAG(M, k), ldaM );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempMm, tempkmin,
                        0., 1.,
                        DIAG(M, k), ldaM );
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempMm,
                            tempkmin, ib, T->nb,
                            DIAG(M, k), ldaM,
                            T(M, k), T->mb,
                            B(m, M), ldbm);
                    }
                    for (n = M+1; n < min(M+BS,  A->mt); n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        ldan = BLKLDD(A, n);
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                tempmm, tempMm, tempmm, tempnn,
                                tempkn, ib, T->nb,
                                B(m, M), ldbm,
                                B(m, n), ldbm,
                                A(n, k), ldan,
                                T(n, k), T->mb);
                        }
                    }
                }
                for (RD = BS; RD < A->mt-k; RD *= 2) {
                    for (M = k; M+RD < A->mt; M += 2*RD) {
                        tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                        ldaMRD = BLKLDD(A, M+RD);
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm   = BLKLDD(B, m);
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempMRDm,
                                tempkn, ib, T->nb,
                                B (m, M   ), ldbm,
                                B (m, M+RD), ldbm,
                                A (M+RD, k), ldaMRD,
                                T2(M+RD, k), T->mb);
                        }
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
