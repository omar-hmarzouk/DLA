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
 * @file tile.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

/**
 *
 * @defgroup Tile
 * @brief Group routines exposed to users for matrices conversion LAPACK-Tile
 *
 */

#include "control/common.h"
#include "control/auxiliary.h"
#include "control/tile.h"

/** ***************************************************************************
 *
 * @ingroup Tile
 *
 *  MORSE_Lapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 ******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Lapack_to_Tile(void *Af77, int LDA, MORSE_desc_t *A)
{
    MORSE_context_t  *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t   request;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Lapack_to_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_Lapack_to_Tile", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    morse_sequence_create(morse, &sequence);
    switch( A->dtyp ) {
#if defined(PRECISION_s)
    case MorseRealFloat:
        morse_pslapack_to_tile(Af77, LDA, A, sequence, &request);
        break;
#endif

#if defined(PRECISION_d)
    case MorseRealDouble:
        morse_pdlapack_to_tile(Af77, LDA, A, sequence, &request);
        break;
#endif

#if defined(PRECISION_c)
    case MorseComplexFloat:
        morse_pclapack_to_tile(Af77, LDA, A, sequence, &request);
        break;
#endif

#if defined(PRECISION_z)
    case MorseComplexDouble:
        morse_pzlapack_to_tile(Af77, LDA, A, sequence, &request);
        break;
#endif

    default:
        morse_error("MORSE_Lapack_to_Tile", "Type unknown");
    }
    RUNTIME_barrier(morse);
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/** ***************************************************************************
 *
 * @ingroup Tile
 *
 *  MORSE_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 ******************************************************************************
 *
 * @param[out] A
 *          Descriptor of the MORSE matrix in tile layout.
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Tile_to_Lapack(MORSE_desc_t *A, void *Af77, int LDA)
{
    MORSE_context_t  *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t   request;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Tile_to_Lapack", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_Tile_to_Lapack", "invalid descriptor");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    morse_sequence_create(morse, &sequence);
    switch( A->dtyp ) {
#if defined(PRECISION_s)
    case MorseRealFloat:
        morse_pstile_to_lapack(A, Af77, LDA, sequence, &request);
        break;
#endif

#if defined(PRECISION_d)
    case MorseRealDouble:
        morse_pdtile_to_lapack(A, Af77, LDA, sequence, &request);
        break;
#endif

#if defined(PRECISION_c)
    case MorseComplexFloat:
        morse_pctile_to_lapack(A, Af77, LDA, sequence, &request);
        break;
#endif

#if defined(PRECISION_z)
    case MorseComplexDouble:
        morse_pztile_to_lapack(A, Af77, LDA, sequence, &request);
        break;
#endif

    default:
        morse_error("MORSE_Tile_to_Lapack", "Type unknown");
    }
    RUNTIME_barrier(morse);
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}
