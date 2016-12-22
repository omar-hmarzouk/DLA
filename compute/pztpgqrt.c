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
 * @file pztpgqrt.c
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

#define V(m,n) V,  m,  n
#define T(m,n) T,  m,  n
#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n

/***************************************************************************//**
 *  Parallel tile QR factorization - dynamic scheduling
 **/
void morse_pztpgqrt( int L, MORSE_desc_t *V, MORSE_desc_t *T, MORSE_desc_t *A, MORSE_desc_t *B,
                     MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int ldak, ldvm, ldbm;
    int tempkn, tempnn, tempmm, templm;
    int ib;

    /* Dimension of the first column */
    int maxm  = B->m - L;
    int maxmt = (maxm % B->mb == 0) ? (maxm / B->mb) : (maxm / B->mb + 1);
    int maxmtk;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    /*
     * ztsmqr = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * ztsmqr = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(MORSE_Complex64_t);
    ws_host   *= sizeof(MORSE_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = V->nt-1; k >= 0; k--) {
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);

        maxmtk = chameleon_min( B->mt, maxmt+k ) - 1;
        for (m = maxmtk; m > -1; m--) {
            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
            templm = m == maxmtk  ? tempmm       : 0;
            ldvm = BLKLDD(V, m);
            ldbm = BLKLDD(B, m);

            for (n = k; n < B->nt; n++) {
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                MORSE_TASK_ztpmqrt(
                    &options,
                    MorseLeft, MorseConjTrans,
                    tempmm, tempnn, tempkn, templm, ib, T->nb,
                    V(m, k), ldvm,
                    T(m, k), T->mb,
                    A(k, n), ldak,
                    B(m, n), ldbm );
            }
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
}
