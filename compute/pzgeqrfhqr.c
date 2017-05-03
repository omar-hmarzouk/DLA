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
 * @file pzgeqrfhqr.c
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

#define A(m,n) A,  (m),  (n)
#define TS(m,n) TS,  (m),  (n)
#define TT(m,n) TT,  (m), (n)
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(m,n) DIAG, (m), (n)
#else
#define DIAG(m,n) A,  (m),  (n)
#endif

/***************************************************************************//**
                                                                              *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
                                                                              **/
void morse_pzgeqrfhqr( const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT,
                      MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n, i, j, p;
    int K, M, RD;
    int ldap, ldam, ldaMRD;
    int tempkmin, tempkn, tempMm, tempnn, tempmm, tempMRDm;
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

    tiles = (int*)malloc((qrtree->mt)*sizeof(int));
    memset( tiles, 0, (qrtree->mt)*sizeof(int) );

#if defined(CHAMELEON_USE_MAGMA)
    /* Worker space
     *
     * zgeqrt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     * ztsqrt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     */
    ws_worker = chameleon_max( ws_worker, ib * (ib + A->nb) );

    /* Host space
     *
     * zgeqrt = ib * (A->nb+3*ib) + A->nb )
     * ztsqrt = 2 * ib * (A->nb+ib) + A->nb
     */
    ws_host = chameleon_max( ws_host, ib * (A->mb + 3 * ib) + A->mb );
    ws_host = chameleon_max( ws_host,  2 * ib * (A->nb + ib) + A->nb );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

#if defined(CHAMELEON_COPY_DIAG)
    {
        /* necessary to avoid dependencies between tasks regarding the diag tile */
        int nblk = ( A->mt + BS -1 ) / BS;
        DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
        morse_zdesc_alloc_diag(*DIAG, A->mb, A->nb, nblk * A->mb, A->nb, 0, 0, nblk * A->mb, A->nb, A->p, A->q);
    }
#endif

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
                A(m, k), ldam,
                TS(m, k), TS->mb);
            if ( k < (A->nt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
                MORSE_TASK_zlacpy(
                    &options,
                    MorseLower, tempmm, A->nb, A->nb,
                    A(m, k), ldam,
                    DIAG(m, k), ldam );
#if defined(CHAMELEON_USE_CUDA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempmm, A->nb,
                    0., 1.,
                    DIAG(m, k), ldam );
#endif
#endif
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zunmqr(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, tempkmin, ib, TS->nb,
                    DIAG(m, k), ldam,
                    TS(m, k), TS->mb,
                    A(m, n), ldam);
            }
        }

        /* Setting the order of the tiles */
        libhqr_treewalk(qrtree, k, tiles);

        for (j = k; j < A->mt-1; j++) {
            m = tiles[j];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            ldap = BLKLDD(A, p);

            /* Tiles killed is a TS */
            if(qrtree->gettype(qrtree, k, m) == 0){
                MORSE_TASK_ztsqrt(
                    &options,
                    tempmm, tempkn, ib, TS->nb,
                    A(p, k), ldap,
                    A(m, k), ldam,
                    TS(m, k), TS->mb);

                for (n = k+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    MORSE_TASK_ztsmqr(
                        &options,
                        MorseLeft, MorseConjTrans,
                        A->nb, tempnn, tempmm, tempnn, A->nb, ib, TS->nb,
                        A(p, n), ldap,
                        A(m, n), ldam,
                        A(m, k), ldam,
                        TS(m, k), TS->mb);
                }
            }

            /* Tiles killed is a TT */
            else {
                MORSE_TASK_zttqrt(
                    &options,
                    tempmm, tempkn, ib, TT->nb,
                    A(p, k), ldap,
                    A(m, k), ldam,
                    TT(m, k), TT->mb);

                for (n = k+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    MORSE_TASK_zttmqr(
                        &options,
                        MorseLeft, MorseConjTrans,
                        A->mb, tempnn, tempmm, tempnn, A->nb, ib, TT->nb,
                        A(p, n), ldap,
                        A(m, n), ldam,
                        A(m, k), ldam,
                        TT(m, k), TT->mb);
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
    morse_desc_mat_free(DIAG);
    free(DIAG);
#endif
    (void)DIAG;
}
