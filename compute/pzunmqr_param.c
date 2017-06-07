/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file pzunmqr_param.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2017-05-17
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define TS(m,n) TS,  m,  n
#define TT(m,n) TT,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define D(m,n)   D,  m,  n
#else
#define D(m,n)   A,  m,  n
#endif

/**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 */
void morse_pzunmqr_param(const libhqr_tree_t *qrtree,
                         MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int ldam, ldan, ldbm, ldbp;
    int tempnn, tempkmin, tempmm, tempkn;
    int ib, K;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    K = chameleon_min(A->mt, A->nt);

    /*
     * zunmqr = A->nb * ib
     * ztsmqr = A->nb * ib
     * zttmqr = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of tiles */

    tiles = (int*)malloc((qrtree->mt)*sizeof(int));
    memset( tiles, 0, (qrtree->mt)*sizeof(int) );

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    if (side == MorseLeft ) {
        if (trans == MorseConjTrans) {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkn   = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempmm, tempkmin, A->nb,
                        A(m, k), ldam,
                        D(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempmm, tempkmin,
                        0., 1.,
                        D(m, k), ldam );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( m, k), ldam,
                            TS(m, k), TS->mb,
                            B( m, n), ldbm);
                    }
                }
                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = k; i < B->mt-1; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);
                    if(qrtree->gettype(qrtree, k, m) == 0){
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkn, ib, TS->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( m, k), ldam,
                                TS(m, k), TS->mb);
                        }
                    }
                    else {
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkn, ib, TT->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( m, k), ldam,
                                TT(m, k), TT->mb);
                        }
                    }
                }
                RUNTIME_iteration_pop(morse);
            }
        } else {
            /*
             *  MorseLeft / MorseNoTrans
             */
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = B->mt-2; i >= k; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);

                    /* TT or TS */

                    if(qrtree->gettype(qrtree, k, m) == 0){
                        for (n = k; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkn, ib, TS->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( m, k), ldam,
                                TS(m, k), TS->mb);
                        }
                    }
                    else {
                        for (n = k; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkn, ib, TT->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( m, k), ldam,
                                TT(m, k), TT->mb);
                        }
                    }
                }
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempmm, tempkmin, A->nb,
                        A(m, k), ldam,
                        D(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempmm, tempkmin,
                        0., 1.,
                        D(m, k), ldam );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( m, k), ldam,
                            TS(m, k), TS->mb,
                            B( m, n), ldbm);
                    }
                }
                RUNTIME_iteration_pop(morse);
            }
        }
    } else {
        if (trans == MorseConjTrans) {
            /*
             *  MorseRight / MorseConjTrans
             */
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == A->nt-1 ? A->n - k*A->nb : A->nb;

                /* Setting the order of tiles */
                libhqr_treewalk(qrtree, k, tiles);

                for (i = B->nt-2; i >= k; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    ldbp = BLKLDD(B, p);

                    /* TS or TT */
                    if(qrtree->gettype(qrtree, k, n) == 0){
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkn, ib, TS->nb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( n, k), ldan,
                                TS(n, k), TS->mb);
                        }
                    }
                    else{
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkn, ib, TT->nb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( n, k), ldan,
                                TT(n, k), TT->mb);
                        }
                    }
                }
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    n = qrtree->getm(qrtree, k, i);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    tempkmin = chameleon_min(tempnn, tempkn);
                    ldan = BLKLDD(A, n);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempnn, tempkmin, A->nb,
                        A(n, k), ldan,
                        D(n, k), ldan );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempnn, tempkmin,
                        0., 1.,
                        D(n, k), ldan );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( n, k), ldan,
                            TS(n, k), TS->mb,
                            B( m, n), ldbm);
                    }
                }

                RUNTIME_iteration_pop(morse);
            }
        } else {
            /*
             *  MorseRight / MorseNoTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;

                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    n = qrtree->getm(qrtree, k, i);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    tempkmin = chameleon_min(tempnn, tempkn);
                    ldan = BLKLDD(A, n);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseLower, tempnn, tempkmin, A->nb,
                        A(n, k), ldan,
                        D(n, k), ldan );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseUpper, tempnn, tempkmin,
                        0., 1.,
                        D(n, k), ldan );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmqr(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( n, k), ldan,
                            TS(n, k), TS->mb,
                            B( m, n), ldbm);
                    }
                }
                /* Setting the order of tiles */
                libhqr_treewalk(qrtree, k, tiles);

                for (i = k; i < B->nt-1; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldan = BLKLDD(A, n);
                    ldbp = BLKLDD(B, p);
                    if(qrtree->gettype(qrtree, k, n) == 0){
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_ztsmqr(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkn, ib, TS->nb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( n, k), ldan,
                                TS(n, k), TS->mb);
                        }
                    }
                    else {
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_zttmqr(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkn, ib, TT->nb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( n, k), ldan,
                                TT(n, k), TT->mb);
                        }
                    }
                }

                RUNTIME_iteration_pop(morse);
            }
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();

    (void)D;
}
