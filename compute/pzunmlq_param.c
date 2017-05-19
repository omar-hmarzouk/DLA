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
 * @file pzunmlq_param.c
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
 *  Parallel application of Q using tile V - LQ factorization - dynamic scheduling
 */
void morse_pzunmlq_param(const libhqr_tree_t *qrtree,
                         MORSE_enum side, MORSE_enum trans,
                         MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *TS, MORSE_desc_t *TT,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *D = NULL;

    int k, m, n, i, p;
    int ldan, ldam, ldbm, ldbn, ldak, ldbp;
    int tempnn, temppn, tempkmin, tempmm, tempkn, tempkm;
    int ib, K;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    K = chameleon_min(A->mt, A->nt);

    /*
     * zunmlq = A->nb * ib
     * ztsmlq = A->nb * ib
     * zttmlq = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq = A->nb * ib
     * ztsmlq = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of tiles */
    tiles = (int*)malloc((qrtree->nt)*sizeof(int));
    memset( tiles, 0, (qrtree->nt)*sizeof(int) );

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* necessary to avoid dependencies between tasks regarding the diag tile */
#if defined(CHAMELEON_COPY_DIAG)
    D = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc_diag(*D, A->mb, A->nb, K*A->mb, A->nb, 0, 0, K*A->mb, A->nb, A->p, A->q);
#endif

    if (side == MorseLeft ) {
        if (trans == MorseNoTrans) {
            /*
             *  MorseLeft / MorseNoTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);
                    ldbp = BLKLDD(B, p);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            temppn, tempnn, tempkmin, ib, TS->nb,
                            D( k, p), ldak,
                            TS(k, p), TS->mb,
                            B( p, n), ldbp);
                    }
                }

                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = k; i < A->nt-1; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbp = BLKLDD(B, p);
                    ldbm = BLKLDD(B, m);

                    /* TT or TS */
                    if(qrtree->gettype(qrtree, k, m) == 0){
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkm, ib, TS->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( k, m), ldak,
                                TS(k, m), TS->mb);
                        }
                    }
                    else {
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;

                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkm, ib, TT->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( k, m), ldak,
                                TT(k, m), TS->mb);
                        }
                    }
                }
                RUNTIME_iteration_pop(morse);
            }
        } else {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = A->nt-2; i >= k; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbp = BLKLDD(B, p);
                    ldbm = BLKLDD(B, m);

                    /* TT or TS */
                    if(qrtree->gettype(qrtree, k, m) == 0){
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkm, ib, TS->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( k, m), ldak,
                                TS(k, m), TS->mb);
                        }
                    }
                    else {
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkm, ib, TT->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( k, m), ldak,
                                TT(k, m), TT->mb);
                        }
                    }
                }
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);
                    ldbp = BLKLDD(B, p);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmim, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            temppn, tempnn, tempkmin, ib, TS->nb,
                            D( k, p), ldak,
                            TS(k, p), TS->mb,
                            B( p, n), ldbp);
                    }
                }
                RUNTIME_iteration_pop(morse);
            }
        }
    } else {
        if (trans == MorseNoTrans) {
            /*
             *  MorseRight / MorseNoTrans
             */
            for (k = K-1; k >= 0; k--) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = A->nt-2; i >= k; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldbp = BLKLDD(B, p);

                    /* TT or TS */

                    if(qrtree->gettype(qrtree, k, n) == 0){
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkm, ib, TS->nb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TS(k, n), TS->mb);
                        }
                    }
                    else {
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkm, ib, TT->nb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TT(k, n), TT->mb);
                        }
                    }
                }
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->mt-1 ? A->m-p*A->mb : A->mb;
                    tempkmin = chameleon_min(tempkm, temppn);
                    ldbp = BLKLDD(B, p);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmim, temppn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, temppn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, temppn, tempkmin, ib, TS->nb,
                            D( k, p), ldak,
                            TS(k, p), TS->mb,
                            B( m, p), ldbm);
                    }
                }
                RUNTIME_iteration_pop(morse);
            }
        } else {
            /*
             *  MorseRight / MorseConjTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
                    tempkmin = chameleon_min(tempkm, temppn);
                    ldbp = BLKLDD(B, p);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempkpn, A->nb,
                        A(k, p), ldak,
                        D(k, p), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempkpn,
                        0., 1.,
                        D(k, p), ldak );
#endif
#endif
                    for (m = 0; m < B->mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, temppn, tempkmin, ib, TS->nb,
                            D( k, p), ldak,
                            TS(k, p), TS->mb,
                            B( m, p), ldbm);
                    }
                }
                /* Setting the order of tiles */
                libhqr_treewalk(qrtree, k, tiles);

                for (i = k; i < A->nt-1; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                    ldbp = BLKLDD(B, p);

                    if(qrtree->gettype(qrtree, k, n) == 0){
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_ztsmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkm, ib, TS->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TS(k, n), TS->mb);
                        }
                    }
                    else {
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                tempmm, B->nb, tempmm, tempnn, tempkm, ib, TT->nb,
                                B( p, n), ldbp,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TT(k, n), TT->mb);
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

#if defined(CHAMELEON_COPY_DIAG)
    MORSE_Sequence_Wait(sequence);
    morse_desc_mat_free(D);
    free(D);
#endif
    (void)D;
}
