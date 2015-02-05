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
 * @file zgesv_nopiv.c
 *
 *  MORSE computational routines
 *  Release Date: November, 15th 2009
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
 *  MORSE_zgesv_nopiv - Computes the solution to a system of linear equations A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 *  The tile LU decomposition with partial tile pivoting and row interchanges is used to factor A.
 *  The factored form of A is then used to solve the system of equations A * X = B.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[in] NRHS
 *          The number of right hand sides, i.e., the number of columns of the matrix B.
 *          NRHS >= 0.
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,N).
 *
 * @param[out] descL
 *          On exit, auxiliary factorization data, related to the tile L factor,
 *          necessary to solve the system of equations.
 *
 * @param[out] IPIV
 *          On exit, the pivot indices that define the permutations (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa MORSE_zgesv_nopiv_Tile
 * @sa MORSE_zgesv_nopiv_Tile_Async
 * @sa MORSE_cgesv_nopiv
 * @sa MORSE_dgesv_nopiv
 * @sa MORSE_sgesv_nopiv
 *
 ******************************************************************************/
int MORSE_zgesv_nopiv(int N, int NRHS,
                      MORSE_Complex64_t *A, int LDA,
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
        morse_error("MORSE_zgesv_nopiv", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (N < 0) {
        morse_error("MORSE_zgesv_nopiv", "illegal value of N");
        return -1;
    }
    if (NRHS < 0) {
        morse_error("MORSE_zgesv_nopiv", "illegal value of NRHS");
        return -2;
    }
    if (LDA < max(1, N)) {
        morse_error("MORSE_zgesv_nopiv", "illegal value of LDA");
        return -4;
    }
    if (LDB < max(1, N)) {
        morse_error("MORSE_zgesv_nopiv", "illegal value of LDB");
        return -8;
    }
    /* Quick return */
    if (min(N, NRHS) == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNB */
    status = morse_tune(MORSE_FUNC_ZGESV, N, N, NRHS);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zgesv_nopiv", "morse_tune() failed");
        return status;
    }

    /* Set NT & NTRHS */
    NB    = MORSE_NB;

    morse_sequence_create(morse, &sequence);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,    sequence, &request,
                             morse_desc_mat_free(&(descA)) );
        morse_zooplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS, sequence, &request,
                             morse_desc_mat_free(&(descA)); morse_desc_mat_free(&(descB)));
/*    } else {*/
/*        morse_ziplap2tile( descA, A, NB, NB, LDA, N,    0, 0, N, N,   */
/*                            sequence, &request);*/
/*        morse_ziplap2tile( descB, B, NB, NB, LDB, NRHS, 0, 0, N, NRHS,*/
/*                            sequence, &request);*/
/*    }*/

    /* Call the tile interface */
    MORSE_zgesv_nopiv_Tile_Async(&descA, &descB, sequence, &request);

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
 *  MORSE_zgesv_nopiv_Tile - Solves a system of linear equations using the tile LU factorization.
 *  Tile equivalent of MORSE_zgetrf_nopiv().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the N-by-N coefficient matrix A.
 *          On exit, the tile L and U factors from the factorization (not equivalent to LAPACK).
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS matrix of right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval >0 if i, U(i,i) is exactly zero. The factorization has been completed,
 *               but the factor U is exactly singular, so the solution could not be computed.
 *
 *******************************************************************************
 *
 * @sa MORSE_zgesv_nopiv
 * @sa MORSE_zgesv_nopiv_Tile_Async
 * @sa MORSE_cgesv_nopiv_Tile
 * @sa MORSE_dgesv_nopiv_Tile
 * @sa MORSE_sgesv_nopiv_Tile
 * @sa MORSE_zcgesv_Tile
 *
 ******************************************************************************/
int MORSE_zgesv_nopiv_Tile(MORSE_desc_t *A, MORSE_desc_t *B)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgesv_nopiv_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    MORSE_zgesv_nopiv_Tile_Async(A, B, sequence, &request);
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
 *  MORSE_zgesv_nopiv_Tile_Async - Solves a system of linear equations using the tile
 *  LU factorization.
 *  Non-blocking equivalent of MORSE_zgesv_nopiv_Tile().
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
 * @sa MORSE_zgesv_nopiv
 * @sa MORSE_zgesv_nopiv_Tile
 * @sa MORSE_cgesv_nopiv_Tile_Async
 * @sa MORSE_dgesv_nopiv_Tile_Async
 * @sa MORSE_sgesv_nopiv_Tile_Async
 * @sa MORSE_zcgesv_Tile_Async
 *
 ******************************************************************************/
int MORSE_zgesv_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *B,
                                 MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgesv_nopiv_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgesv_nopiv_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgesv_nopiv_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zgesv_nopiv_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_zgesv_nopiv_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb || B->nb != B->mb) {
        morse_error("MORSE_zgesv_nopiv_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (min(N, NRHS) == 0)
        return MORSE_SUCCESS;
*/

    morse_pzgetrf_nopiv(A, sequence, request);

    morse_pztrsm(MorseLeft, MorseLower, MorseNoTrans, MorseUnit, (MORSE_Complex64_t)1.0, A, B, sequence, request);

    morse_pztrsm(MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, (MORSE_Complex64_t)1.0, A, B, sequence, request);

    return MORSE_SUCCESS;
}
