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
 * @file ztile.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *          If MORSE_TRANSLATION_MODE is set to MORSE_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          MORSE_TRANSLATION_MODE is set to MORSE_OUTOFPLACE,
 *          A->mat has to be allocated before.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zLapack_to_Tile_Async
 * @sa MORSE_zTile_to_Lapack
 * @sa MORSE_cLapack_to_Tile
 * @sa MORSE_dLapack_to_Tile
 * @sa MORSE_sLapack_to_Tile
 *
 ******************************************************************************/
int MORSE_zLapack_to_Tile(MORSE_Complex64_t *Af77, int LDA, MORSE_desc_t *A)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zLapack_to_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check( A ) != MORSE_SUCCESS) {
        morse_error("MORSE_zLapack_to_Tile", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    morse_sequence_create(morse, &sequence);

    morse_pzlapack_to_tile( Af77, LDA, A, sequence, &request);

    RUNTIME_barrier( morse );
    RUNTIME_desc_getoncpu( A );

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_zLapack_to_Tile_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of MORSE_zLapack_to_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *          If MORSE_TRANSLATION_MODE is set to MORSE_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          MORSE_TRANSLATION_MODE is set to MORSE_OUTOFPLACE,
 *          A->mat has to be allocated before.
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
 * @sa MORSE_zTile_to_Lapack_Async
 * @sa MORSE_zLapack_to_Tile
 * @sa MORSE_cLapack_to_Tile_Async
 * @sa MORSE_dLapack_to_Tile_Async
 * @sa MORSE_sLapack_to_Tile_Async
 *
 ******************************************************************************/
int MORSE_zLapack_to_Tile_Async(MORSE_Complex64_t *Af77, int LDA, MORSE_desc_t *A,
                                  MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zLapack_to_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check( A ) != MORSE_SUCCESS) {
        morse_error("MORSE_zLapack_to_Tile", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }

    morse_pzlapack_to_tile( Af77, LDA, A, sequence, request);

    return MORSE_SUCCESS;
}

/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If MORSE_TRANSLATION_MODE is set to MORSE_INPLACE,
 *          Af77 has to be A->mat, else if
 *          MORSE_TRANSLATION_MODE is set to MORSE_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zTile_to_Lapack_Async
 * @sa MORSE_zLapack_to_Tile
 * @sa MORSE_cTile_to_Lapack
 * @sa MORSE_dTile_to_Lapack
 * @sa MORSE_sTile_to_Lapack
 *
******************************************************************************/
int MORSE_zTile_to_Lapack(MORSE_desc_t *A, MORSE_Complex64_t *Af77, int LDA)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zTile_to_Lapack", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check( A ) != MORSE_SUCCESS) {
        morse_error("MORSE_zTile_to_Lapack", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    morse_sequence_create(morse, &sequence);

    morse_pztile_to_lapack( A, Af77, LDA, sequence, &request);
    RUNTIME_barrier( morse );
    RUNTIME_desc_getoncpu( A );
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/*******************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_zTile_to_Lapack_Async - Conversion from LAPACK layout to tile layout.
 *  Non-blocking equivalent of MORSE_zTile_to_Lapack().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations ar runtime.
 *
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If MORSE_TRANSLATION_MODE is set to MORSE_INPLACE,
 *          Af77 has to be A->mat, else if
 *          MORSE_TRANSLATION_MODE is set to MORSE_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
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
 * @sa MORSE_zLapack_to_Tile_Async
 * @sa MORSE_zTile_to_Lapack
 * @sa MORSE_cTile_to_Lapack_Async
 * @sa MORSE_dTile_to_Lapack_Async
 * @sa MORSE_sTile_to_Lapack_Async
 *
 ******************************************************************************/
int MORSE_zTile_to_Lapack_Async(MORSE_desc_t *A, MORSE_Complex64_t *Af77, int LDA,
                                MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zTile_to_Lapack", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check( A ) != MORSE_SUCCESS) {
        morse_error("MORSE_zTile_to_Lapack", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }

    morse_pztile_to_lapack( A, Af77, LDA, sequence, request );

    return MORSE_SUCCESS;
}
