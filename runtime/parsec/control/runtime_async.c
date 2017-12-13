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
#include "chameleon_parsec.h"

/*******************************************************************************
 *  Create a sequence
 **/
int RUNTIME_sequence_create(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    parsec_context_t  *parsec        = (parsec_context_t *)morse->schedopt;
    parsec_taskpool_t *parsec_dtd_tp = parsec_dtd_taskpool_new();

    parsec_enqueue(parsec, (parsec_taskpool_t *) parsec_dtd_tp);
    sequence->schedopt = parsec_dtd_tp;

    parsec_context_start(parsec);

    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Destroy a sequence
 **/
int RUNTIME_sequence_destroy(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    parsec_context_t  *parsec = (parsec_context_t *)morse->schedopt;
    parsec_taskpool_t *parsec_dtd_tp = (parsec_taskpool_t *) sequence->schedopt;
    (void)morse;

    assert( parsec_dtd_tp );

    // TODO: switch to a partial wait
    //parsec_dtd_taskpool_wait(parsec, parsec_dtd_tp);
    parsec_context_wait(parsec);

    parsec_taskpool_free(parsec_dtd_tp);

    sequence->schedopt = NULL;
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Wait for the completion of a sequence
 **/
int RUNTIME_sequence_wait(MORSE_context_t *morse, MORSE_sequence_t *sequence )
{
    parsec_context_t  *parsec = (parsec_context_t *)morse->schedopt;
    parsec_taskpool_t *parsec_dtd_tp = (parsec_taskpool_t *) sequence->schedopt;

    assert( parsec_dtd_tp );

    parsec_dtd_taskpool_wait(parsec, parsec_dtd_tp);

    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Terminate a sequence
 **/
void RUNTIME_sequence_flush(void *schedopt, MORSE_sequence_t *sequence, MORSE_request_t *request, int status)
{
    parsec_context_t *parsec = (parsec_context_t *)schedopt;
    (void)schedopt;
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}

