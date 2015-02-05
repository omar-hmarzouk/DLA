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
 * @file pzgelqfrh.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
//ALLOC_WS :  A->nb + ib*T->nb
//WS_ADD :  A->nb + ib*T->nb
#include "control/common.h"

#define A(m,n) A,  (m),  (n)
#define T(m,n) T,  (m),  (n)
#define T2(m,n) T,  (m),  (n)+A->nt
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(m,n) DIAG, ((n)/BS), 0
#else
#define DIAG(m,n) A,  (m),  (n)
#endif

/***************************************************************************//**
 *  Parallel tile LQ factorization (reduction Householder) - dynamic scheduling
 **/
void morse_pzgelqfrh(MORSE_desc_t *A, MORSE_desc_t *T, int BS,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int N, RD;
    int ldak, ldam;
    int tempkmin, tempkm, tempNn, tempnn, tempmm, tempNRDn;
    int ib;
    int nblk;

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
#if defined(CHAMELEON_USE_MAGMA)
    /* Worker space
     *
     * zgelqt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     * zunmlq = A->nb * ib
     * ztslqt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     * ztsmlq = 2 * A->nb * ib
     */
    ws_worker = max( ws_worker, ib * (ib + A->nb) );
    ws_worker = max( ws_worker, ib * A->nb * 2 );

    /* Host space
     *
     * zgelqt =     ib * A->nb + 3 * ib * ib + A->nb
     * ztslqt = 3 * ib * A->nb +     ib * ib + A->nb
     */
    ws_host = max( ws_host,     ib * A->nb + 3 * ib * ib + A->nb );
    ws_host = max( ws_host, 3 * ib * A->nb +     ib * ib + A->nb );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* necessary to avoid dependencies between tasks regarding the diag tile */
    nblk = ( A->nt + BS -1 ) / BS;
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, nblk * A->mb, A->nb, 0, 0, nblk * A->mb, A->nb);

    for (k = 0; k < min(A->mt, A->nt); k++) {
        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        ldak = BLKLDD(A, k);
        for (N = k; N < A->nt; N += BS) {
            tempNn = N == A->nt-1 ? A->n-N*A->nb : A->nb;
            tempkmin = min(tempkm, tempNn);
            MORSE_TASK_zgelqt(
                &options,
                tempkm, tempNn, ib, T->nb,
                A(k, N), ldak,
                T(k, N), T->mb);
        if ( k < (A->mt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
            MORSE_TASK_zlacpy(
                &options,
                MorseUpper, tempkm, tempNn, A->nb,
                A(k, N), ldak,
                DIAG(k, N), ldak );
#endif
#if defined(CHAMELEON_USE_MAGMA)
            MORSE_TASK_zlaset(
                &options,
                MorseLower, tempkm, tempNn,
                0., 1.,
                DIAG(k, N), ldak );
#endif
        }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_zunmlq(
                    &options,
                    MorseRight, MorseConjTrans,
                    tempmm, tempNn, tempkmin, ib, T->nb,
                    DIAG(k, N), ldak,
                    T(k, N), T->mb,
                    A(m, N), ldam);
            }
            for (n = N+1; n < min(N+BS, A->nt); n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztslqt(
                    &options,
                    tempkm, tempnn, ib, T->nb,
                    A(k, N), ldak,
                    A(k, n), ldak,
                    T(k, n), T->mb);

                for (m = k+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_ztsmlq(
                        &options,
                        MorseRight, MorseConjTrans,
                        tempmm, A->nb, tempmm, tempnn, tempkm, ib, T->nb,
                        A(m, N), ldam,
                        A(m, n), ldam,
                        A(k, n), ldak,
                        T(k, n), T->mb);
                }
            }
        }
        for (RD = BS; RD < A->nt-k; RD *= 2) {
            for (N = k; N+RD < A->nt; N += 2*RD) {
                tempNRDn = N+RD == A->nt-1 ? A->n-(N+RD)*A->nb : A->nb;
                MORSE_TASK_zttlqt(
                    &options,
                    tempkm, tempNRDn, ib, T->nb,
                    A (k, N   ), ldak,
                    A (k, N+RD), ldak,
                    T2(k, N+RD), T->mb);

                for (m = k+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                    ldam   = BLKLDD(A, m );
                    MORSE_TASK_zttmlq(
                        &options,
                        MorseRight, MorseConjTrans,
                        tempmm, A->nb, tempmm, tempNRDn, tempkm, ib, T->nb,
                        A (m, N   ), ldam,
                        A (m, N+RD), ldam,
                        A (k, N+RD), ldak,
                        T2(k, N+RD), T->mb);
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
