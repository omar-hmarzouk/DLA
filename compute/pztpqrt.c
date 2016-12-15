/**
 *
 * @copyright (c) 2009-2016 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                          Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file pztpqrt.c
 *
 *  MORSE computational routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define T(m,n) T,  m,  n
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(k) DIAG, k, 0
#else
#define DIAG(k) A, k, k
#endif

/***************************************************************************//**
 *  Parallel tile QR factorization - dynamic scheduling
 **/
void morse_pztpqrt( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T,
                    MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;
    MORSE_desc_t *DIAG = NULL;

    int k, m, n;
    int ldak, ldbm;
    int tempkm, tempkn, tempnn, tempmm, templm;
    int ib;

    /* Dimension of the first column */
    int maxm  = B->m - L;
    int maxmt = (maxm % B->mb == 0) ? (maxm / B->mb) : (maxm / B->mb + 1);

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * zgeqrt = A->nb * (ib+1)
     * zunmqr = A->nb * ib
     * ztsqrt = A->nb * (ib+1)
     * ztsmqr = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr = A->nb * ib
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = max( ws_worker, ib * A->nb * 2 );
#endif

#if defined(CHAMELEON_USE_MAGMA)
    /* Worker space
     *
     * zgeqrt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     * ztsqrt = max( A->nb * (ib+1), ib * (ib + A->nb) )
     */
    ws_worker = max( ws_worker, ib * (ib + A->nb) );

    /* Host space
     *
     * zgeqrt = ib * (A->mb+3*ib) + A->mb )
     * ztsqrt = 2 * ib * (A->nb+ib) + A->nb
     */
    ws_host = max( ws_host, ib * (A->mb + 3 * ib) + A->mb );
    ws_host = max( ws_host,  2 * ib * (A->nb + ib) + A->nb );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

#if defined(CHAMELEON_COPY_DIAG)
    /* necessary to avoid dependencies between tsqrt and unmqr tasks regarding the diag tile */
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc_diag(*DIAG, A->mb, A->nb, min(A->m, A->n), A->nb, 0, 0, min(A->m, A->n), A->nb, A->p, A->q);
#endif

    for (k = 0; k < A->nt; k++) {
        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);

        for (m = 0; m < maxmt; m++) {
            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
            templm = m == maxmt-1 ? tempmm       : 0;
            ldbm = BLKLDD(B, m);
            MORSE_TASK_ztpqrt(
                &options,
                tempmm, tempkn, templm, ib, T->nb,
                A(k, k), ldak,
                B(m, k), ldbm,
                T(m, k), T->mb );

            for (n = k+1; n < B->nt; n++) {
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, tempkm, templm, ib, T->nb,
                    B(m, k), ldbm,
                    T(m, k), T->mb,
                    A(k, n), ldak,
                    B(m, n), ldbm );
            }
        }

        maxmt = min( B->mt, maxmt+1 );
    }
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
