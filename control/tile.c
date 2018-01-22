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
    MORSE_desc_t B;
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

    /* Create the B descriptor to handle the Lapack format matrix */
    B = morse_desc_init_user(
        A->dtyp, A->mb, A->nb, A->bsiz,
        LDA, A->n, 0, 0, A->m, A->n, 1, 1,
        morse_getaddr_cm, morse_getblkldd_cm, NULL );
    B.mat  = Af77;
    B.styp = MorseCM;

    RUNTIME_desc_create( &B );

    morse_sequence_create(morse, &sequence);

    morse_pzlacpy( MorseUpperLower, &B, A, sequence, &request );

    RUNTIME_desc_flush( &B, sequence );
    RUNTIME_desc_flush(  A, sequence );
    RUNTIME_sequence_wait( morse, sequence );

    RUNTIME_desc_destroy( &B );

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
 *          LAPACK matrix (only needed on proc 0).
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
    MORSE_desc_t      B;
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

    /* Create the B descriptor to handle the Lapack format matrix */
    B = morse_desc_init_user(
        A->dtyp, A->mb, A->nb, A->bsiz,
        LDA, A->n, 0, 0, A->m, A->n, 1, 1,
        morse_getaddr_cm, morse_getblkldd_cm, NULL );
    B.mat  = Af77;
    B.styp = MorseCM;

    RUNTIME_desc_create( &B );

    morse_sequence_create(morse, &sequence);
    morse_pzlacpy( MorseUpperLower, A, &B, sequence, &request );

    RUNTIME_desc_flush(  A, sequence );
    RUNTIME_desc_flush( &B, sequence );
    RUNTIME_sequence_wait( morse, sequence );

    RUNTIME_desc_destroy( &B );

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}
