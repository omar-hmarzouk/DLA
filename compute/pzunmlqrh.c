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
 * @file pzunmlqrh.c
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
#define T2(m,n) T,  (m),  (n)+A->nt
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(m,n) DIAG, ((n)/BS), 0
#else
#define DIAG(m,n) A, (m), (n)
#endif

/***************************************************************************//**
 *  Parallel application of Q using tile V - LQ factorization (reduction
 *  Householder) - dynamic scheduling
 **/
void morse_pzunmlqrh(MORSE_enum side, MORSE_enum trans,
                     MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, int BS,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int K, N, RD, lastRD;
    int ldaN, ldak;
    int ldbN, ldbm, ldbNRD;
    int tempNn, tempkm, tempnn, tempmm, tempNRDn, tempkmin;
    int ib;
    int nblk;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zunmlq = A->nb * ib
     * ztsmlq = A->nb * ib
     * zttmlq = A->nb * ib
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
    nblk = ( A->nt + BS -1 ) / BS;
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, nblk * A->mb, A->nb, 0, 0, nblk * A->mb, A->nb);

    K = min(A->mt, A->nt);
    if (side == MorseLeft ) {
        if (trans == MorseNoTrans) {
            /*
             *  MorseLeft / MorseNoTrans
             */
            for (k = 0; k < K; k++) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                for (N = k; N < A->nt; N += BS) {
                    tempNn   = N == A->nt-1 ? A->n-N*A->nb : A->nb;
                    tempkmin = min(tempkm,tempNn);
                    ldaN = BLKLDD(A, N);
                    ldbN = BLKLDD(B, N);
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempNn, A->nb,
                        A(k, N), ldak,
                        DIAG(k, N), ldak );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempNn,
                        0., 1.,
                        DIAG(k, N), ldak );
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempNn, tempnn,
                            tempkmin, ib, T->nb,
                            DIAG(k, N), ldak,
                            T(k, N), T->mb,
                            B(N, n), ldbN);
                    }
                    for (m = N+1; m < min(N+BS, A->nt); m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                B->nb, tempnn, tempmm, tempnn,
                                tempkm, ib, T->nb,
                                B(N, n), ldbN,
                                B(m, n), ldbm,
                                A(k, m), ldak,
                                T(k, m), T->mb);
                        }
                    }
                }
                for (RD = BS; RD < A->nt-k; RD *= 2) {
                    for (N = k; N+RD < A->nt; N += 2*RD) {
                        tempNRDn = N+RD == A->nt-1 ? A->n-(N+RD)*A->nb : A->nb;
                        ldbN   = BLKLDD(B, N   );
                        ldbNRD = BLKLDD(B, N+RD);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempNRDn, tempnn,
                                tempkm, ib, T->nb,
                                B (N,    n), ldbN,
                                B (N+RD, n), ldbNRD,
                                A (k, N+RD), ldak,
                                T2(k, N+RD), T->mb);
                        }
                    }
                }
            }
        } else {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = K-1; k >= 0; k--) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                lastRD = 0;
                for (RD = BS; RD < A->nt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (N = k; N+RD < A->nt; N += 2*RD) {
                        tempNRDn = N+RD == A->nt-1 ? A->n-(N+RD)*A->nb : A->nb;
                        ldbN   = BLKLDD(B, N   );
                        ldbNRD = BLKLDD(B, N+RD);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                B->nb, tempnn, tempNRDn, tempnn,
                                tempkm, ib, T->nb,
                                B (N,    n), ldbN,
                                B (N+RD, n), ldbNRD,
                                A (k, N+RD), ldak,
                                T2(k, N+RD), T->mb);
                        }
                    }
                }
                for (N = k; N < A->nt; N += BS) {
                    tempNn   = N == A->nt-1 ? A->n-N*A->nb : A->nb;
                    tempkmin = min(tempkm,tempNn);
                    ldaN = BLKLDD(A, N);
                    ldbN = BLKLDD(B, N);
                    for (m = min(N+BS, A->nt)-1; m > N; m--) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn,
                                tempkm, ib, T->nb,
                                B(N, n), ldbN,
                                B(m, n), ldbm,
                                A(k, m), ldak,
                                T(k, m), T->mb);
                        }
                    }
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempNn, A->nb,
                        A(k, N), ldak,
                        DIAG(k, N), ldak );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempNn,
                        0., 1.,
                        DIAG(k, N), ldak );
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempNn, tempnn,
                            tempkmin, ib, T->nb,
                            DIAG(k, N), ldak,
                            T(k, N), T->mb,
                            B(N, n), ldbN);
                    }
                }
            }

        }
    } else {
        if (trans == MorseNoTrans) {
            /*
             *  MorseRight / MorseNoTrans
             */
            for (k = K-1; k >= 0; k--) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                lastRD = 0;
                for (RD = BS; RD < A->nt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (N = k; N+RD < A->nt; N += 2*RD) {
                        tempNRDn = N+RD == A->nt-1 ? A->n-(N+RD)*A->nb : A->nb;
                        for (m = 0; m < B->mt; m++) {
                            ldbm   = BLKLDD(B, m);
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempNRDn,
                                tempkm, ib, T->nb,
                                B (m, N   ), ldbm,
                                B (m, N+RD), ldbm,
                                A (k, N+RD), ldak,
                                T2(k, N+RD), T->mb);
                        }
                    }
                }
                for (N = k; N < A->nt; N += BS) {
                    tempNn   = N == A->nt-1 ? A->n-N*A->nb : A->nb;
                    tempkmin = min(tempkm,tempNn);
                    for (n = min(N+BS, A->nt)-1; n > N; n--) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn,
                                tempkm, ib, T->nb,
                                B(m, N), ldbm,
                                B(m, n), ldbm,
                                A(k, n), ldak,
                                T(k, n), T->mb);
                        }
                    }
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempNn, A->nb,
                        A(k, N), ldak,
                        DIAG(k, N), ldak );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempNn,
                        0., 1.,
                        DIAG(k, N), ldak );
#endif
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, tempNn,
                            tempkmin, ib, T->nb,
                            DIAG(k, N), ldak,
                            T(k, N), T->mb,
                            B(m, N), ldbm);
                    }
                }
            }
        } else {
            /*
             *  MorseRight / MorseConjTrans
             */
            for (k = 0; k < K; k++) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                for (N = k; N < A->nt; N += BS) {
                    tempNn = N == A->nt-1 ? A->n-N*A->nb : A->nb;
                    tempkmin = min(tempkm,tempNn);
                    ldaN = BLKLDD(A, N);
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempNn, A->nb,
                        A(k, N), ldaN,
                        DIAG(k, N), ldaN );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempNn,
                        0., 1.,
                        DIAG(k, N), ldaN );
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, tempNn,
                            tempkmin, ib, T->nb,
                            DIAG(k, N), ldaN,
                            T(k, N), T->mb,
                            B(m, N), ldbm);
                    }
                    for (n = N+1; n < min(N+BS, A->nt); n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                tempmm, tempNn, tempmm, tempnn,
                                tempkm, ib, T->nb,
                                B(m, N), ldbm,
                                B(m, n), ldbm,
                                A(k, n), ldak,
                                T(k, n), T->mb);
                        }
                    }
                }
                for (RD = BS; RD < A->nt-k; RD *= 2) {
                    for (N = k; N+RD < A->nt; N += 2*RD) {
                        tempNRDn = N+RD == A->nt-1 ? A->n-(N+RD)*A->nb : A->nb;
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm   = BLKLDD(B, m);
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempNRDn,
                                tempkm, ib, T->nb,
                                B (m, N   ), ldbm,
                                B (m, N+RD), ldbm,
                                A (k, N+RD), ldak,
                                T2(k, N+RD), T->mb);
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
