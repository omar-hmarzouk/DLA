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
 * @file zunglq.c
 *
 *  MORSE computational routines
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

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zunglq - Generates an M-by-N matrix Q with orthonormal rows, which is defined as the
 *  first M rows of a product of the elementary reflectors returned by MORSE_zgelqf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix Q. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix Q. N >= M.
 *
 * @param[in] K
 *          The number of rows of elementary tile reflectors whose product defines the matrix Q.
 *          M >= K >= 0.
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 * @param[in] LDQ
 *          The leading dimension of the array Q. LDQ >= max(1,M).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval MORSE_SUCCESS <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zunglq_Tile
 * @sa MORSE_zunglq_Tile_Async
 * @sa MORSE_cunglq
 * @sa MORSE_dorglq
 * @sa MORSE_sorglq
 * @sa MORSE_zgelqf
 *
 ******************************************************************************/
int MORSE_zunglq(int M, int N, int K,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_desc_t *descT,
                  MORSE_Complex64_t *Q, int LDQ)
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descA, descQ;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunglq", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zunglq", "illegal value of M");
        return -1;
    }
    if (N < M) {
        morse_error("MORSE_zunglq", "illegal value of N");
        return -2;
    }
    if (K < 0 || K > M) {
        morse_error("MORSE_zunglq", "illegal value of K");
        return -3;
    }
    if (LDA < max(1, M)) {
        morse_error("MORSE_zunglq", "illegal value of LDA");
        return -5;
    }
    if (LDQ < max(1, M)) {
        morse_error("MORSE_zunglq", "illegal value of LDQ");
        return -8;
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDQ ) */
    if (min(M, min(N, K)) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create(morse, &sequence);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K, N, sequence, &request,
                             morse_desc_mat_free(&(descA)) );
        morse_zooplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, M, N, sequence, &request,
                             morse_desc_mat_free(&(descA)); morse_desc_mat_free(&(descQ)));
/*    } else {*/
/*        morse_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K, N,*/
/*                            sequence, &request);*/
/*        morse_ziplap2tile( descQ, Q, NB, NB, LDQ, N, 0, 0, M, N,*/
/*                            sequence, &request);*/
/*    }*/

    /* Call the tile interface */
    MORSE_zunglq_Tile_Async(&descA, descT, &descQ, sequence, &request);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooptile2lap(descQ, Q, NB, NB, LDQ, N,  sequence, &request);
        RUNTIME_barrier(morse);
        morse_desc_mat_free(&descA);
        morse_desc_mat_free(&descQ);
/*    } else {*/
/*        morse_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);*/
/*        morse_ziptile2lap( descQ, Q, NB, NB, LDQ, N,  sequence, &request);*/
/*        RUNTIME_barrier(morse);*/
/*    }*/
        
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 * MORSE_zunglq_Tile - Generates an M-by-N matrix Q with orthonormal rows, which is defined as the
 * first M rows of a product of the elementary reflectors returned by MORSE_zgelqf.
 * All matrices are passed through descriptors. All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Details of the LQ factorization of the original matrix A as returned by MORSE_zgelqf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by MORSE_zgelqf.
 *
 * @param[out] Q
 *          On exit, the M-by-N matrix Q.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zunglq
 * @sa MORSE_zunglq_Tile_Async
 * @sa MORSE_cunglq_Tile
 * @sa MORSE_dorglq_Tile
 * @sa MORSE_sorglq_Tile
 * @sa MORSE_zgelqf_Tile
 *
 ******************************************************************************/
int MORSE_zunglq_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *Q)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunglq_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    MORSE_zunglq_Tile_Async(A, T, Q, sequence, &request);
    RUNTIME_barrier(morse);
    RUNTIME_desc_getoncpu(A);
        RUNTIME_desc_getoncpu(Q);
    
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  Non-blocking equivalent of MORSE_zunglq_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa MORSE_zunglq
 * @sa MORSE_zunglq_Tile
 * @sa MORSE_cunglq_Tile_Async
 * @sa MORSE_dorglq_Tile_Async
 * @sa MORSE_sorglq_Tile_Async
 * @sa MORSE_zgelqf_Tile_Async
 *
 ******************************************************************************/
int MORSE_zunglq_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *Q,
                             MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zunglq_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zunglq_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zunglq_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(Q) != MORSE_SUCCESS) {
        morse_error("MORSE_zunglq_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || Q->nb != Q->mb) {
        morse_error("MORSE_zunglq_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return - currently NOT equivalent to LAPACK's:
     * CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, Q, LDQ ) */
/*
    if (min(M, N) == 0)
        return MORSE_SUCCESS;
*/
    if (morse->householder == MORSE_FLAT_HOUSEHOLDER) {
        morse_pzunglq(A, Q, T, sequence, request);
    }
    else {
        morse_pzunglqrh(A, Q, T, MORSE_RHBLK, sequence, request);
    }

    return MORSE_SUCCESS;
}
