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
 * @file pzgetrf_incpiv.c
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
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
//ALLOC_WS :  ib*L->nb
//WS_ADD :  ib*L->nb
#include "control/common.h"

#define A(_m_,_n_) A, _m_, _n_
#if defined(CHAMELEON_COPY_DIAG)
#define DIAG(_k_) DIAG, _k_, 0
#else
#define DIAG(_k_) A, _k_, _k_
#endif
#define L(_m_,_n_) L,  _m_,  _n_
#define IPIV(_m_,_n_) &(IPIV[(int64_t)A->mb*((int64_t)(_m_)+(int64_t)A->mt*(int64_t)(_n_))])

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling
 **/
void morse_pzgetrf_incpiv(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_desc_t *DIAG = NULL;
    MORSE_context_t *morse;
    MORSE_option_t options;
    size_t h_work_size, d_work_size;

    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;
    int minMNT = min(A->mt, A->nt);

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;
#if defined(CHAMELEON_USE_MAGMA)
    h_work_size  = sizeof(MORSE_Complex64_t)*( 2*ib + 2*L->nb )*2*A->mb;
    d_work_size  = sizeof(MORSE_Complex64_t)*(   ib           )*2*A->mb;
#else
    h_work_size  = sizeof(MORSE_Complex64_t)*( ib*L->nb );
    d_work_size  = 0;
#endif
    RUNTIME_options_ws_alloc( &options, h_work_size, d_work_size );

    /* necessary to avoid dependencies between tasks regarding the diag tile */
    DIAG = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    morse_zdesc_alloc2(*DIAG, A->mb, A->nb, min(A->m, A->n), A->nb, 0, 0, min(A->m, A->n), A->nb);

    for (k = 0; k < minMNT; k++) {
        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);
        MORSE_TASK_zgetrf_incpiv(
            &options,
            tempkm, tempkn, ib, L->nb,
            A(k, k), ldak,
            L(k, k), L->mb,
            IPIV(k, k),
            k == A->mt-1, A->nb*k);

        if ( k < (minMNT-1) ) {
#if defined(CHAMELEON_COPY_DIAG)
            MORSE_TASK_zlacpy(
                &options,
                MorseUpperLower, tempkm, tempkn, A->nb,
                A(k, k), ldak,
                DIAG(k), ldak);
#endif
        }

        for (n = k+1; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            MORSE_TASK_zgessm(
                &options,
                tempkm, tempnn, tempkm, ib, L->nb,
                IPIV(k, k),
                L(k, k), L->mb,
                DIAG(k), ldak,
                A(k, n), ldak);
        }
        for (m = k+1; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            MORSE_TASK_ztstrf(
                &options,
                tempmm, tempkn, ib, L->nb,
                A(k, k), ldak,
                A(m, k), ldam,
                L(m, k), L->mb,
                IPIV(m, k),
                m == A->mt-1, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zssssm(
                    &options,
                    A->nb, tempnn, tempmm, tempnn, A->nb, ib, L->nb,
                    A(k, n), ldak,
                    A(m, n), ldam,
                    L(m, k), L->mb,
                    A(m, k), ldam,
                    IPIV(m, k));
            }
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();

    morse_desc_mat_free(DIAG);
    free(DIAG);
}
