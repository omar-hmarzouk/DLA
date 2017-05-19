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
 * @file pzunglq_pram.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2011-05-24
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  (m),  (n)
#define Q(m,n) Q,  (m),  (n)
#define TS(m,n) TS,  (m),  (n)
#define TT(m,n) TT,  (m),  (n)
#if defined(CHAMELEON_COPY_DIAG)
#define D(m,n) D, ((n)/BS), 0
#else
#define D(m,n) A, (m), (n)
#endif

/**
 *  Parallel construction of Q using tile V - dynamic scheduling
 */
void morse_pzunglq_param(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *Q,
                        MORSE_desc_t *TS, MORSE_desc_t *TT,
                        MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *D = NULL;

    int k, m, n, i, p;
    int K;
    int ldak, ldqp, ldqm;
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

#if defined(CHAMELEON_COPY_DIAG)
    {
        /* necessary to avoid dependencies between tasks regarding the diag tile */
        int nblk = ( A->nt + BS -1 ) / BS;
        D = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
        morse_zdesc_alloc_diag(*DIAG, A->mb, A->nb, nblk * A->mb, A->nb, 0, 0, nblk * A->mb, A->nb, A->p, A->q);
    }
#endif

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
            ldqp = BLKLDD(Q, p);

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

            temppn = p == A->mt-1 ? A->m-p*A->mb : A->mb;
            tempkmin = chameleon_min(tempkm, temppn);
            ldqp = BLKLDD(Q, p);

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
            for (m = k; m < Q->mt; m++) {
                tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
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

#if defined(CHAMELEON_COPY_DIAG)
    MORSE_Sequence_Wait(sequence);
    morse_desc_mat_free(D);
    free(D);
#endif
    (void)D;
}
