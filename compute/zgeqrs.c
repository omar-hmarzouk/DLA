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
 * @file zgeqrs.c
 *
 *  MORSE computational routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
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
 *  MORSE_zgeqrs - Compute a minimum-norm solution min || A*X - B || using the RQ factorization
 *  A = R*Q computed by MORSE_zgeqrf.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= M >= 0.
 *
 * @param[in] NRHS
 *          The number of columns of B. NRHS >= 0.
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by MORSE_zgeqrf.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= M.
 *
 * @param[in] descT
 *          Auxiliary factorization data, computed by MORSE_zgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_zgeqrs_Tile
 * @sa MORSE_zgeqrs_Tile_Async
 * @sa MORSE_cgeqrs
 * @sa MORSE_dgeqrs
 * @sa MORSE_sgeqrs
 * @sa MORSE_zgeqrf
 *
 *******************************************************************************/
int MORSE_zgeqrs(int M, int N, int NRHS,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_desc_t *descT,
                  MORSE_Complex64_t *B, int LDB)
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descA, descB;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgeqrs", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_zgeqrs", "illegal value of M");
        return -1;
    }
    if (N < 0 || N > M) {
        morse_error("MORSE_zgeqrs", "illegal value of N");
        return -2;
    }
    if (NRHS < 0) {
        morse_error("MORSE_zgeqrs", "illegal value of N");
        return -3;
    }
    if (LDA < max(1, M)) {
        morse_error("MORSE_zgeqrs", "illegal value of LDA");
        return -5;
    }
    if (LDB < max(1, max(1, M))) {
        morse_error("MORSE_zgeqrs", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (min(M, min(N, NRHS)) == 0) {
        return MORSE_SUCCESS;
    }

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGELS, M, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create(morse, &sequence);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N,    sequence, &request,
                             morse_desc_mat_free(&(descA)) );
        morse_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS, sequence, &request,
                             morse_desc_mat_free(&(descA)); morse_desc_mat_free(&(descB)));
/*    } else {*/
/*        morse_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, M, N,   */
/*                            sequence, &request);*/
/*        morse_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, M, NRHS,*/
/*                            sequence, &request);*/
/*    }*/

    /* Call the tile interface */
    MORSE_zgeqrs_Tile_Async(&descA, descT, &descB, sequence, &request);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooptile2lap(descA, A, NB, NB, LDA, N,     sequence, &request);
        morse_zooptile2lap(descB, B, NB, NB, LDB, NRHS,  sequence, &request);
        RUNTIME_barrier(morse);
        morse_desc_mat_free(&descA);
        morse_desc_mat_free(&descB);
/*    } else {*/
/*        morse_ziptile2lap( descA, A, NB, NB, LDA, N,     sequence, &request);*/
/*        morse_ziptile2lap( descB, B, NB, NB, LDB, NRHS,  sequence, &request);*/
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
 *  MORSE_zgeqrs_Tile - Computes a minimum-norm solution using the tile QR factorization.
 *  Tile equivalent of MORSE_zgeqrf().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Details of the QR factorization of the original matrix A as returned by MORSE_zgeqrf.
 *
 * @param[in] T
 *          Auxiliary factorization data, computed by MORSE_zgeqrf.
 *
 * @param[in,out] B
 *          On entry, the m-by-nrhs right hand side matrix B.
 *          On exit, the n-by-nrhs solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zgeqrs
 * @sa MORSE_zgeqrs_Tile_Async
 * @sa MORSE_cgeqrs_Tile
 * @sa MORSE_dgeqrs_Tile
 * @sa MORSE_sgeqrs_Tile
 * @sa MORSE_zgeqrf_Tile
 *
 ******************************************************************************/
int MORSE_zgeqrs_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgeqrs_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    MORSE_zgeqrs_Tile_Async(A, T, B, sequence, &request);
    RUNTIME_barrier(morse);
    RUNTIME_desc_getoncpu(A);
        RUNTIME_desc_getoncpu(B);
    
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_zgeqrs_Tile_Async - Computes a minimum-norm solution using the tile
 *  QR factorization.
 *  Non-blocking equivalent of MORSE_zgeqrs_Tile().
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
 * @sa MORSE_zgeqrs
 * @sa MORSE_zgeqrs_Tile
 * @sa MORSE_cgeqrs_Tile_Async
 * @sa MORSE_dgeqrs_Tile_Async
 * @sa MORSE_sgeqrs_Tile_Async
 * @sa MORSE_zgeqrf_Tile_Async
 *
 ******************************************************************************/
int MORSE_zgeqrs_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B,
                             MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_desc_t *subA;
    MORSE_desc_t *subB;
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgeqrs_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgeqrs_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgeqrs_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zgeqrs_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_zgeqrs_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(M, min(N, NRHS)) == 0) {
        return MORSE_SUCCESS;
    }
*/
    if (morse->householder == MORSE_FLAT_HOUSEHOLDER) {
        morse_pzunmqr(MorseLeft, MorseConjTrans, A, B, T, sequence, request);
    }
    else {
        morse_pzunmqrrh(MorseLeft, MorseConjTrans, A, B, T, MORSE_RHBLK, sequence, request);
    }

    subB = morse_desc_submatrix(B, 0, 0, A->n, B->n);
    subA = morse_desc_submatrix(A, 0, 0, A->n, A->n);
    morse_pztrsm(MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, 1.0, subA, subB, sequence, request);
    free(subA);
    free(subB);

    return MORSE_SUCCESS;
}
