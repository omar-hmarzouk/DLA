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
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Raphael Boucherie
 * @date 2010-11-15
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
    int ldam, ldan, ldbm, ldak, ldbp;
    int tempnn, tempkmin, tempmm, tempkn, tempkm;
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
    morse_zdesc_alloc_diag(*D, A->mb, A->nb, K*A->nb, A->nb, 0, 0, K*A->nb, A->nb, A->p, A->q);
#endif

    if (side == MorseLeft ) {
        if (trans == MorseConjTrans) {
            /*
             *  MorseLeft / MorseConjTrans
             */
            for (k = 0; k < K; k++) {
                RUNTIME_iteration_push(morse, k);

                tempkm   = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkm);
                    ldbm = BLKLDD(B, m);
#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempmm, A->nb,
                        A(k, m), ldak,
                        D(k, m), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempmm, tempkmin,
                        0., 1.,
                        D(k, m), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( k, m), ldak,
                            TS(k, m), TS->mb,
                            B( m, n), ldbm);
                    }
                }
                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = k; i < B->mt-1; i++) {
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
                                B->nb, tempnn, tempmm, tempnn, tempkm, ib, TS->nb,
                                B( m, p), ldbp,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TS(k, n), TS->mb);
                        }
                    }
                    else {
                        for (m = 0; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            ldbm = BLKLDD(B, m);
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkn, ib, TT->mb,
                                B( m, p), ldbm,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TT(k, n), TT->mb);
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

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                /* Setting the order of the tiles*/
                libhqr_treewalk(qrtree, k, tiles);

                for (i = B->mt-2; i >= k; i--) {
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
                                B->nb, tempnn, tempmm, tempnn, tempkm, ib, TS->nb,
                                B( m, p), ldbp,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TS(k, n), TS->mb);
                        }
                    }
                    else {
                        for (m = k; m < B->mt; m++) {
                            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                            MORSE_TASK_zttmlq(
                                &options,
                                side, trans,
                                B->mb, tempnn, tempmm, tempnn, tempkm, ib, TT->nb,
                                B( m, p), ldbp,
                                B( m, n), ldbm,
                                A( k, n), ldak,
                                TT(k, n), TT->mb);
                        }
                    }
                }
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    tempkmin = chameleon_min(tempmm, tempkn);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmim, tempmm, A->nb,
                        A(k, m), ldak,
                        D(k, m), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempmm,
                        0., 1.,
                        D(k, m), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( k, m), ldak,
                            TS(k, m), TS->mb,
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

                tempkm = k == A->mt-1 ? A->m - k*A->mb : A->mb;
                ldak = BLKLDD(A, k);

                /* Setting the order of tiles */
                libhqr_treewalk(qrtree, k, tiles);

                for (i = B->nt-2; i >= k; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);

                    /* TS or TT */
                    if(qrtree->gettype(qrtree, k, m) == 0){
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
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
                    else{
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
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
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    tempkmin = chameleon_min(tempmm, tempkm);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempkmm, A->nb,
                        A(k, m), ldak,
                        D(k, m), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempkmm,
                        0., 1.,
                        D(k, m), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        ldbm = BLKLDD(B, m);
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( k, n), ldak,
                            TS(k, n), TS->mb,
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

                tempkm = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                ldak = BLKLDD(A, k);
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    m = qrtree->getm(qrtree, k, i);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    tempkmin = chameleon_min(tempmm, tempkm);
                    ldbm = BLKLDD(B, m);

#if defined(CHAMELEON_COPY_DIAG)
                    MORSE_TASK_zlacpy(
                        &options,
                        MorseUpper, tempkmin, tempkmm, A->nb,
                        A(k, m), ldak,
                        D(k, m), ldak );
#if defined(CHAMELEON_USE_CUDA)
                    MORSE_TASK_zlaset(
                        &options,
                        MorseLower, tempkmin, tempkmm,
                        0., 1.,
                        D(k, m), ldak );
#endif
#endif
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_zunmlq(
                            &options,
                            side, trans,
                            tempmm, tempnn, tempkmin, ib, TS->nb,
                            D( k, m), ldan,
                            TS(k, m), TS->mb,
                            B( m, n), ldbm);
                    }
                }
                /* Setting the order of tiles */
                libhqr_treewalk(qrtree, k, tiles);

                for (i = k; i < B->mt-1; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                    ldbm = BLKLDD(B, m);
                    ldbp = BLKLDD(B, p);
                    if(qrtree->gettype(qrtree, k, m) == 0){
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
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
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
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
