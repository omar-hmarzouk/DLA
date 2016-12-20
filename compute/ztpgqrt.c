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
 * @file ztpgqrt.c
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

/**
 ******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_ztpgqrt - Generates a partial Q matrix formed with a blocked QR
 *  factorization of a "triangular-pentagonal" matrix C, which is composed of an
 *  unused triangular block and a pentagonal block V, using the compact
 *  representation for Q. See MORSE_ztpqrt() to generate V.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrices B, and V. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrices B, and A. N >= 0.
 *
 * @param[in] K
 *          The number of elementary reflectors whose product defines
 *          the matrix Q in the matrix V.
 *
 * @param[in] L
 *          The number of rows of the upper trapezoidal part of V.
 *          MIN(M,N) >= L >= 0.  See Further Details.
 *
 * @param[in] V
 *          The i-th row must contain the vector which defines the
 *          elementary reflector H(i), for i = 1,2,...,k, as returned by
 *          MORSE_ztpqrt() in the first k rows of its array argument V.
 *          V is matrx of size M-by-K. The first M-L rows
 *          are rectangular, and the last L rows are upper trapezoidal.
 *
 * @param[in] LDV
 *          The leading dimension of the array V. LDV >= max(1,K).
 *
 * @param[int] descT
 *          On exit, auxiliary factorization data, required by MORSE_zgeqrs to
 *          solve the system of equations, or by any function to apply the Q.
 *
 * @param[in,out] A
 *          A is COMPLEX*16 array, dimension (LDA,N)
 *          On entry, the K-by-N matrix A.
 *          On exit, A is overwritten by the corresponding block of
 *          Q*A.  See Further Details.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,K).
 *
 * @param[in,out] B
 *          On entry, the pentagonal M-by-N matrix B.
 *          On exit, B contains Q.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.  LDB >= max(1,M).
 *
 * @par Further Details:
 * =====================
 *
 *  The input matrix Q is a (K+M)-by-N matrix
 *
 *               Q = [ A ]
 *                   [ B ]
 *
 *  where A is an identity matrix, and B is a M-by-N matrix of 0.
 *  V a matrix of householder reflectors with a pentagonal shape consisting of a
 *  (M-L)-by-K rectangular matrix V1 on top of a L-by-N
 *  Upper trapezoidal matrix V2:
 *
 *               V = [ V1 ]  <- (M-L)-by-N rectangular
 *                   [ V2 ]  <-     L-by-N upper trapezoidal.
 *
 *  The upper trapezoidal matrix V2 consists of the first L rows of a
 *  K-by-K upper triangular matrix, where 0 <= L <= MIN(M,K).  If L=0,
 *  V is rectangular M-by-K; if M=L=K, V is upper triangular.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa MORSE_ztpgqrt_Tile
 * @sa MORSE_ztpgqrt_Tile_Async
 * @sa MORSE_ctpgqrt
 * @sa MORSE_dtpgqrt
 * @sa MORSE_stpgqrt
 * @sa MORSE_zgeqrs
 *
 ******************************************************************************/
int MORSE_ztpgqrt( int M, int N, int K, int L,
                   MORSE_Complex64_t *V, int LDV,
                   MORSE_desc_t *descT,
                   MORSE_Complex64_t *A, int LDA,
                   MORSE_Complex64_t *B, int LDB )
{
    int NB;
    int status;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    MORSE_desc_t descA, descB, descV;
    int minMK = min( M, K );

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztpgqrt", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if (M < 0) {
        morse_error("MORSE_ztpgqrt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        morse_error("MORSE_ztpgqrt", "illegal value of N");
        return -2;
    }
    if (K < 0) {
        morse_error("MORSE_ztpgqrt", "illegal value of K");
        return -3;
    }
    if ((L < 0) || ((L > minMK) && (minMK > 0))) {
        morse_error("MORSE_ztpgqrt", "illegal value of N");
        return -4;
    }
    if (LDV < max(1, M)) {
        morse_error("MORSE_ztpgqrt", "illegal value of LDV");
        return -6;
    }
    if (LDA < max(1, K)) {
        morse_error("MORSE_ztpgqrt", "illegal value of LDA");
        return -9;
    }
    if (LDB < max(1, M)) {
        morse_error("MORSE_ztpgqrt", "illegal value of LDB");
        return -11;
    }

    /* Quick return */
    if (minMK == 0)
        return MORSE_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZGELS, M, K, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_ztpgqrt", "morse_tune() failed");
        return status;
    }

    /* Set NT */
    NB = MORSE_NB;

    morse_sequence_create(morse, &sequence);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooplap2tile( descV, V, NB, NB, LDB, K, 0, 0, M, K, sequence, &request,
                            morse_desc_mat_free(&(descV)) );
        morse_zooplap2tile( descA, A, NB, NB, LDA, N, 0, 0, K, N, sequence, &request,
                            (morse_desc_mat_free(&(descV)),
                             morse_desc_mat_free(&(descA))) );
        morse_zooplap2tile( descB, B, NB, NB, LDB, N, 0, 0, M, N, sequence, &request,
                            (morse_desc_mat_free(&(descV)),
                             morse_desc_mat_free(&(descA)),
                             morse_desc_mat_free(&(descB))) );
/*    } else {*/
/*        morse_ziplap2tile( descA, A, NB, NB, LDA, N, 0, 0, M, N,*/
/*                            sequence, &request);*/
/*    }*/

    /* Call the tile interface */
    MORSE_ztpgqrt_Tile_Async(L, &descV, descT, &descA, &descB, sequence, &request);

/*    if ( MORSE_TRANSLATION == MORSE_OUTOFPLACE ) {*/
        morse_zooptile2lap(descA, A, NB, NB, LDA, N, sequence, &request);
        morse_zooptile2lap(descB, B, NB, NB, LDB, N, sequence, &request);
        morse_sequence_wait(morse, sequence);
        morse_desc_mat_free(&descV);
        morse_desc_mat_free(&descA);
        morse_desc_mat_free(&descB);
/*    } else {*/
/*        morse_ziptile2lap( descV, V, NB, NB, LDV, K,  sequence, &request);*/
/*        morse_ziptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);*/
/*        morse_ziptile2lap( descB, B, NB, NB, LDB, N,  sequence, &request);*/
/*        morse_sequence_wait(morse, sequence);*/
/*    }*/

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_ztpgqrt_Tile - Generates a partial Q matrix formed with a blocked QR
 *  factorization of a "triangular-pentagonal" matrix C, which is composed of an
 *  unused triangular block and a pentagonal block V, using the compact
 *  representation for Q. See MORSE_ztpqrt() to generate V.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array contain the min(M,N)-by-N
 *          upper trapezoidal matrix R (R is upper triangular if M >= N); the elements below the
 *          diagonal represent the unitary matrix Q as a product of elementary reflectors stored
 *          by tiles.
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by MORSE_zgeqrs to solve the system
 *          of equations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_ztpgqrt
 * @sa MORSE_ztpgqrt_Tile_Async
 * @sa MORSE_ctpgqrt_Tile
 * @sa MORSE_dtpgqrt_Tile
 * @sa MORSE_stpgqrt_Tile
 * @sa MORSE_zgeqrs_Tile
 *
 ******************************************************************************/
int MORSE_ztpgqrt_Tile( int L, MORSE_desc_t *V, MORSE_desc_t *T, MORSE_desc_t *A, MORSE_desc_t *B )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztpgqrt_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    MORSE_ztpgqrt_Tile_Async(L, V, T, A, B, sequence, &request);
    morse_sequence_wait(morse, sequence);
    RUNTIME_desc_getoncpu(A);
    RUNTIME_desc_getoncpu(B);

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/**
 *******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_ztpgqrt_Tile_Async - Generates a partial Q matrix formed with a blocked QR
 *  factorization of a "triangular-pentagonal" matrix C, which is composed of an
 *  unused triangular block and a pentagonal block V, using the compact
 *  representation for Q. See MORSE_ztpqrt() to generate V.
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
 * @sa MORSE_ztpgqrt
 * @sa MORSE_ztpgqrt_Tile
 * @sa MORSE_ctpgqrt_Tile_Async
 * @sa MORSE_dtpgqrt_Tile_Async
 * @sa MORSE_stpgqrt_Tile_Async
 * @sa MORSE_zgeqrs_Tile_Async
 *
 ******************************************************************************/
int MORSE_ztpgqrt_Tile_Async( int L, MORSE_desc_t *V, MORSE_desc_t *T, MORSE_desc_t *A, MORSE_desc_t *B,
                              MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_ztpgqrt_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_ztpgqrt_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_ztpgqrt_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(V) != MORSE_SUCCESS) {
        morse_error("MORSE_ztpgqrt_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(T) != MORSE_SUCCESS) {
        morse_error("MORSE_ztpgqrt_Tile", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_ztpgqrt_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(B) != MORSE_SUCCESS) {
        morse_error("MORSE_ztpgqrt_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        morse_error("MORSE_ztpgqrt_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (((B->m - L) % B->mb) != 0) {
        morse_error("MORSE_ztpgqrt_Tile", "Triangular part must be aligned with tiles");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* if (morse->householder == MORSE_FLAT_HOUSEHOLDER) { */
    morse_pztpgqrt(L, V, T, A, B, sequence, request);
    /* } */
    /* else { */
    /*    morse_pztpgqrtrh(A, T, MORSE_RHBLK, sequence, request); */
    /* } */

    return MORSE_SUCCESS;
}
