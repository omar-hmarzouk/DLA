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
 * @file pzungqr.c
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
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
#define Q(m,n) Q,  m,  n
#define T(m,n) T,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(k) DIAG, k, 0
#else
#define DIAG(k) A, k, k
#endif

/***************************************************************************//**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 **/
void morse_pzungqr(MORSE_desc_t *A, MORSE_desc_t *Q, MORSE_desc_t *T,
                   MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int ldak, ldqk, ldam, ldqm;
    int tempmm, tempnn, tempkmin, tempkm;
    int tempAkm, tempAkn;
    int ib;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zunmqr = A->nb * ib
     * ztsmqr = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
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
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, min(A->m, A->n), A->nb, 0, 0, min(A->m, A->n), A->nb);

    for (k = min(A->mt, A->nt)-1; k >= 0; k--) {
        tempAkm  = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempAkn  = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        tempkmin = min( tempAkn, tempAkm );
        tempkm   = k == Q->mt-1 ? Q->m-k*Q->mb : Q->mb;
        ldak = BLKLDD(A, k);
        ldqk = BLKLDD(Q, k);
        for (m = Q->mt - 1; m > k; m--) {
            tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
            ldam = BLKLDD(A, m);
            ldqm = BLKLDD(Q, m);
            for (n = 0; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;
                MORSE_TASK_ztsmqr(
                    &options,
                    MorseLeft, MorseNoTrans,
                    Q->mb, tempnn, tempmm, tempnn, tempAkn, ib, T->nb,
                    Q(k, n), ldqk,
                    Q(m, n), ldqm,
                    A(m, k), ldam,
                    T(m, k), T->mb);
            }
        }
#if defined(CHAMELEON_COPY_DIAG)
        MORSE_TASK_zlacpy(
            &options,
            MorseLower, tempkm, tempkmin, A->nb,
            A(k, k), ldak,
            DIAG(k), ldak );
#endif
#if defined(CHAMELEON_USE_MAGMA)
        MORSE_TASK_zlaset(
            &options,
            MorseUpper, tempkm, tempkmin,
            0., 1.,
            DIAG(k), ldak );
#endif
        for (n = 0; n < Q->nt; n++) {
            tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;
            MORSE_TASK_zunmqr(
                &options,
                MorseLeft, MorseNoTrans,
                tempkm, tempnn, tempkmin, ib, T->nb,
                DIAG(k), ldak,
                T(k, k), T->mb,
                Q(k, n), ldqk);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();

    morse_desc_mat_free(DIAG);
    free(DIAG);
}
