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
 * @file pzgeqrf_param.c
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
#include <stdlib.h>
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
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 */
void morse_pzgeqrf_param( const libhqr_tree_t *qrtree, MORSE_desc_t *A,
                          MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *D,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int K;
    int ldap, ldam;
    int tempkmin, tempkn, tempnn, tempmm;
    int ib;
    int *tiles;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zgeqrt = A->nb * (ib+1)
     * zunmqr = A->nb * ib
     * ztsqrt = A->nb * (ib+1)
     * zttqrt = A->nb * (ib+1)
     * ztsmqr = A->nb * ib
     * zttmqr = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    /* Initialisation of tiles */

    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);

    /* The number of the factorization */
    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(morse, k);
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        /* The number of geqrt to apply */
        for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
            m = qrtree->getm(qrtree, k, i);
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            tempkmin = chameleon_min(tempmm, tempkn);
            ldam = BLKLDD(A, m);

            MORSE_TASK_zgeqrt(
                &options,
                tempmm, tempkn, ib, TS->nb,
                A( m, k), ldam,
                TS(m, k), TS->mb);
            if ( k < (A->nt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempmm, A->nb, A->nb,
                    A(m, k), ldam,
                    D(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempmm, A->nb,
                    0., 1.,
                    D(m, k), ldam );
#endif
#endif
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zunmqr(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, tempkmin, ib, TS->nb,
                    D( m, k), ldam,
                    TS(m, k), TS->mb,
                    A( m, n), ldam);
            }
        }

        /* Setting the order of the tiles */
        libhqr_walk_stepk( qrtree, k, tiles + (k+1) );

        for (i = k+1; i < A->mt; i++) {
            m = tiles[i];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldap = BLKLDD(A, p);
            ldam = BLKLDD(A, m);

            /* Tiles killed is a TS */
            if(qrtree->gettype(qrtree, k, m) == 0){
                MORSE_TASK_ztsqrt(
                    &options,
                    tempmm, tempkn, ib, TS->nb,
                    A( p, k), ldap,
                    A( m, k), ldam,
                    TS(m, k), TS->mb);

                for (n = k+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    MORSE_TASK_ztsmqr(
                        &options,
                        MorseLeft, MorseConjTrans,
                        A->nb, tempnn, tempmm, tempnn, A->nb, ib, TS->nb,
                        A( p, n), ldap,
                        A( m, n), ldam,
                        A( m, k), ldam,
                        TS(m, k), TS->mb);
                }
            }

            /* Tiles killed is a TT */
            else {
                MORSE_TASK_zttqrt(
                    &options,
                    tempmm, tempkn, ib, TT->nb,
                    A( p, k), ldap,
                    A( m, k), ldam,
                    TT(m, k), TT->mb);

                for (n = k+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    MORSE_TASK_zttmqr(
                        &options,
                        MorseLeft, MorseConjTrans,
                        A->mb, tempnn, tempmm, tempnn, A->nb, ib, TT->nb,
                        A( p, n), ldap,
                        A( m, n), ldam,
                        A( m, k), ldam,
                        TT(m, k), TT->mb);
                }
            }
        }
        RUNTIME_iteration_pop(morse);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    (void)D;
}
