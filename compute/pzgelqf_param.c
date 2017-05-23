/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University of
 *                          Tennessee Research Foundation.  All rights reserved.
 * @copyright (c) 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                          Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file pzgelqf_param.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Raphael Boucherie
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"
#include "libhqr.h"

#define A(m,n)  A,  (m), (n)
#define TS(m,n) TS, (m), (n)
#define TT(m,n) TT, (m), (n)
#if defined(CHAMELEON_COPY_DIAG)
#define D(m,n)  D,  (m), (n)
#else
#define D(m,n)  A,  (m), (n)
#endif

/**
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 */
void morse_pzgelqf_param( const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *D = NULL;

    int k, m, n, i, p;
    int K;
    int ldak, ldam, ldap;
    int tempkmin, tempkm, tempnn, tempmm;
    int ib;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zgelqt = A->nb * (ib+1)
     * zunmlq = A->nb * ib
     * ztslqt = A->nb * (ib+1)
     * zttlqt = A->nb * (ib+1)
     * ztsmlq = A->nb * ib
     * zttmlq = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq = A->nb * ib
     * ztsmlq = 2 * A->nb * ib
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
        D = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
        morse_zdesc_alloc(*D, A->mb, A->nb, A->m, A->n, 0, 0, A->m, A->n, );
    }
#endif

    K = chameleon_min(A->mt, A->nt);

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(morse, k);
        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

        /* The number of geqrt to apply */
        for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
            n = qrtree->getm(qrtree, k, i);
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            tempkmin = chameleon_min(tempnn, tempkm);
            ldak = BLKLDD(A, k);

            MORSE_TASK_zgelqt(
                &options,
                tempkm, tempnn, ib, TS->nb,
                A( k, n), ldak,
                TS(k, n), TS->mb);
            if ( k < (A->nt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseUpper, tempkm, tempnn, A->nb,
                    A(k, n), ldak,
                    D(k, n), ldak );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseLower, tempkm, tempnn,
                    0., 1.,
                    D(k, n), ldak );
#endif
#endif
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_zunmlq(
                    &options,
                    MorseRight, MorseConjTrans,
                    tempmm, tempnn, tempkmin, ib, TS->nb,
                    D( k, n), ldak,
                    TS(k, n), TS->mb,
                    A( m, n), ldam);
            }
        }

        /* Setting the order of the tiles */
        libhqr_treewalk(qrtree, k, tiles);

        for (i = k; i < A->nt-1; i++) {
            n = tiles[i];
            p = qrtree->currpiv(qrtree, k, n);

            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            ldap = BLKLDD(A, p);
            /* Tiles killed is a TS */
            if(qrtree->gettype(qrtree, k, n) == 0){
                MORSE_TASK_ztslqt(
                    &options,
                    tempkm, tempnn, ib, TS->nb,
                    A( k, p), ldak,
                    A( k, n), ldak,
                    TS(k, n), TS->mb);

                for (m = k+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_ztsmlq(
                        &options,
                        MorseRight, MorseConjTrans,
                        tempmm, A->nb, tempmm, tempnn, tempkm, ib, TS->nb,
                        A( m, p), ldam,
                        A( m, n), ldam,
                        A( k, n), ldak,
                        TS(k, n), TS->mb);
                }
            }

            /* Tiles killed is a TT */
            else {
                MORSE_TASK_zttlqt(
                    &options,
                    tempkm, tempnn, ib, TT->nb,
                    A( k, p), ldak,
                    A( k, n), ldak,
                    TT(k, n), TT->mb);

                for (m = k+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zttmlq(
                        &options,
                        MorseRight, MorseConjTrans,
                        tempmm, A->nb, tempmm, tempnn, tempkm, ib, TT->nb,
                        A( m, p), ldam,
                        A( m, n), ldam,
                        A( k, n), ldak,
                        TT(k, n), TT->mb);
                }
            }
        }
        RUNTIME_iteration_pop(morse);
    }

    free(tiles);
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
