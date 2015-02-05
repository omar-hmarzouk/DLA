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
 * @file pzgeqrfrh.c
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
#include "control/common.h"

#define A(m,n) A,  (m),  (n)
#define T(m,n) T,  (m),  (n)
#define T2(m,n) T,  (m), ((n)+A->nt)
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(m,n) DIAG, ((m)/BS), 0
#else
#define DIAG(m,n) A,  (m),  (n)
#endif

/***************************************************************************//**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 **/
void morse_pzgeqrfrh(MORSE_desc_t *A, MORSE_desc_t *T, int BS,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int K, M, RD;
    int ldaM, ldam, ldaMRD;
    int tempkmin, tempkn, tempMm, tempnn, tempmm, tempMRDm;
    int ib;
    int nblk;

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
#if defined(CHAMELEON_USE_MAGMA)
    /* Worker space
     *
     * zgeqrt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     * zunmqr = A->nb * ib
     * ztsqrt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = max( ws_worker, ib * (ib + A->nb) );
    ws_worker = max( ws_worker, ib * A->nb * 2 );

    /* Host space
     *
     * zgeqrt = ib * (A->nb+3*ib) + A->nb )
     * ztsqrt = 2 * ib * (A->nb+ib) + A->nb
     */
    ws_host = max( ws_host, ib * (A->mb + 3 * ib) + A->mb );
    ws_host = max( ws_host,  2 * ib * (A->nb + ib) + A->nb );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* necessary to avoid dependencies between tasks regarding the diag tile */
    nblk = ( A->mt + BS -1 ) / BS;
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, nblk * A->mb, A->nb, 0, 0, nblk * A->mb, A->nb);

    K = min(A->mt, A->nt);
    for (k = 0; k < K; k++) {
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        for (M = k; M < A->mt; M += BS) {
            tempMm = M == A->mt-1 ? A->m-M*A->mb : A->mb;
            tempkmin = min(tempMm, tempkn);
            ldaM = BLKLDD(A, M);
            MORSE_TASK_zgeqrt(
                &options,
                tempMm, tempkn, ib, T->nb,
                A(M, k), ldaM,
                T(M, k), T->mb);
            if ( k < (A->nt-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
            MORSE_TASK_zlacpy(
                &options,
                MorseLower, tempMm, A->nb, A->nb,
                A(M, k), ldaM,
                DIAG(M, k), ldaM );
#endif
#if defined(CHAMELEON_USE_MAGMA)
                MORSE_TASK_zlaset(
                    &options,
                    MorseUpper, tempMm, A->nb,
                    0., 1.,
                    DIAG(M, k), ldaM );
#endif
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zunmqr(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempMm, tempnn, tempkmin, ib, T->nb,
                    DIAG(M, k), ldaM,
                    T(M, k), T->mb,
                    A(M, n), ldaM);
            }
            for (m = M+1; m < min(M+BS, A->mt); m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztsqrt(
                    &options,
                    tempmm, tempkn, ib, T->nb,
                    A(M, k), ldaM,
                    A(m, k), ldam,
                    T(m, k), T->mb);

                for (n = k+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    MORSE_TASK_ztsmqr(
                        &options,
                        MorseLeft, MorseConjTrans,
                        A->nb, tempnn, tempmm, tempnn, A->nb, ib, T->nb,
                        A(M, n), ldaM,
                        A(m, n), ldam,
                        A(m, k), ldam,
                        T(m, k), T->mb);
                }
            }
        }
        for (RD = BS; RD < A->mt-k; RD *= 2) {
            for (M = k; M+RD < A->mt; M += 2*RD) {
                tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                ldaM   = BLKLDD(A, M   );
                ldaMRD = BLKLDD(A, M+RD);
                MORSE_TASK_zttqrt(
                    &options,
                    tempMRDm, tempkn, ib, T->nb,
                    A (M   , k), ldaM,
                    A (M+RD, k), ldaMRD,
                    T2(M+RD, k), T->mb);

                for (n = k+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                    MORSE_TASK_zttmqr(
                        &options,
                        MorseLeft, MorseConjTrans,
                        A->nb, tempnn, tempMRDm, tempnn, A->nb, ib, T->nb,
                        A (M,    n), ldaM,
                        A (M+RD, n), ldaMRD,
                        A (M+RD, k), ldaMRD,
                        T2(M+RD, k), T->mb);
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
