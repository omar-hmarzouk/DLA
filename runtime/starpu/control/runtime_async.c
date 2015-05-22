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
 * @file runtime_async.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "runtime/starpu/include/morse_starpu.h"

/*******************************************************************************
 *  Create a sequence
 **/
int RUNTIME_sequence_create( MORSE_context_t *morse, MORSE_sequence_t *sequence )
{
    (void)morse;
    (void)sequence;
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Destroy a sequence
 **/
int RUNTIME_sequence_destroy( MORSE_context_t *morse, MORSE_sequence_t *sequence )
{
    (void)morse;
    (void)sequence;
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Wait for the completion of a sequence
 **/
int RUNTIME_sequence_wait( MORSE_context_t *morse, MORSE_sequence_t *sequence )
{
    (void)morse;
    (void)sequence;
    starpu_task_wait_for_all();
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *  Terminate a sequence
 **/
void RUNTIME_sequence_flush( void *schedopt, MORSE_sequence_t *sequence, MORSE_request_t *request, int status)
{
    (void)schedopt;
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    starpu_task_wait_for_all();
    return;
}
