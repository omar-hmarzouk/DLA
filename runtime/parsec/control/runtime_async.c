/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
#include <stdlib.h>
#include "runtime/parsec/include/morse_parsec.h"

/*******************************************************************************
 *  Create a sequence
 **/
int RUNTIME_sequence_create(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    dague_context_t    *dague = (dague_context_t *)morse->schedopt;
    dague_dtd_handle_t *dague_dtd_handle = dague_dtd_handle_new((dague_context_t *)morse->schedopt);

    dague_enqueue(dague, (dague_handle_t*) dague_dtd_handle);
    sequence->schedopt = dague_dtd_handle;

#if defined (OVERLAP)
    dague_context_start(dague);
#endif
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Destroy a sequence
 **/
int RUNTIME_sequence_destroy(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    dague_context_t    *dague = (dague_context_t *)morse->schedopt;
    dague_dtd_handle_t *dague_dtd_handle = (dague_dtd_handle_t *) sequence->schedopt;
    (void)morse;

    assert( dague_dtd_handle );

    dague_dtd_context_wait_on_handle(dague, dague_dtd_handle);
    dague_dtd_handle_destruct(dague_dtd_handle);
    sequence->schedopt = NULL;
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Wait for the completion of a sequence
 **/
int RUNTIME_sequence_wait(MORSE_context_t *morse, MORSE_sequence_t *sequence )
{
    dague_context_t    *dague = (dague_context_t *)morse->schedopt;
    dague_dtd_handle_t *dague_dtd_handle = (dague_dtd_handle_t *) sequence->schedopt;

    assert( dague_dtd_handle );

    dague_dtd_handle_wait(dague, dague_dtd_handle);

    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Terminate a sequence
 **/
void RUNTIME_sequence_flush(void *schedopt, MORSE_sequence_t *sequence, MORSE_request_t *request, int status)
{
    dague_context_t *dague = (dague_context_t *)schedopt;
    (void)schedopt;
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}

