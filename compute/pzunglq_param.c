/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file pzunglq_pram.c
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

#define A(m,n) A,  (m),  (n)
#define Q(m,n) Q,  (m),  (n)
#define TS(m,n) TS,  (m),  (n)
#define TT(m,n) TT,  (m),  (n)
#if defined(CHAMELEON_COPY_DIAG)
#define D(m,n) D, (m), (n)
#else
#define D(m,n) A, (m), (n)
#endif

/**
 *  Parallel construction of Q using tile V - dynamic scheduling
 */
void morse_pzunglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                         MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int K;
    int ldak, ldqm;
    int tempkm, tempkmin, temppn, tempnn, tempmm;
    int ib;
    int *tiles;

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

    K = chameleon_min(A->mt, A->nt);

    for (k = K-1; k >= 0; k--) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        ldak = BLKLDD(A, k);

        /* Setting the order of the tiles*/
        libhqr_treewalk(qrtree, k, tiles);

        for (i = A->nt-2; i >= k; i--) {
            n = tiles[i];
            p = qrtree->currpiv(qrtree, k, n);

            tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

            /* TT or TS */

            if(qrtree->gettype(qrtree, k, n) == 0){
                for (m = k; m < Q->mt; m++) {
                    tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                    ldqm = BLKLDD(Q, m);
                    MORSE_TASK_ztsmlq(
                        &options,
                        MorseRight, MorseNoTrans,
                        tempmm, Q->nb, tempmm, tempnn, tempkm, ib, TS->nb,
                        Q( m, p), ldqm,
                        Q( m, n), ldqm,
                        A( k, n), ldak,
                        TS(k, n), TS->mb);
                }
            }
            else {
                for (m = k; m < Q->mt; m++) {
                    tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                    ldqm = BLKLDD(Q, m);
                    MORSE_TASK_zttmlq(
                        &options,
                        MorseRight, MorseNoTrans,
                        tempmm, Q->nb, tempmm, tempnn, tempkm, ib, TT->nb,
                        Q( m, p), ldqm,
                        Q( m, n), ldqm,
                        A( k, n), ldak,
                        TT(k, n), TT->mb);
                }
            }
        }
        for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
            p = qrtree->getm(qrtree, k, i);

            temppn = p == A->nt-1 ? A->n-p*A->nb : A->nb;
            tempkmin = chameleon_min(tempkm, temppn);

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
            for (m = k; m < Q->mt; m++) {
                tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
                ldqm = BLKLDD(Q, m);
                MORSE_TASK_zunmlq(
                    &options,
                    MorseRight, MorseNoTrans,
                    tempmm, temppn, tempkmin, ib, TS->nb,
                    D( k, p), ldak,
                    TS(k, p), TS->mb,
                    Q( m, p), ldqm);
            }
        }
        RUNTIME_iteration_pop(morse);
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
    (void)D;
}
