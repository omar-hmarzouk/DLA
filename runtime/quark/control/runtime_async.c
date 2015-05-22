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
 *  MORSE async routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version
 * @author Jakub Kurzak
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2011-10-29
 *
 **/
#include <stdlib.h>
#include "runtime/quark/include/morse_quark.h"

/** *****************************************************************************
 *  Create a sequence
 **/
int RUNTIME_sequence_create(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    sequence->schedopt =(void*) QUARK_Sequence_Create((Quark*)(morse->schedopt));

    if (sequence->schedopt == NULL) {
        morse_error("MORSE_Sequence_Create", "QUARK_Sequence_Create() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }
    sequence->status = MORSE_SUCCESS;
    return MORSE_SUCCESS;
}

/** *****************************************************************************
 *  Destroy a sequence
 **/
int RUNTIME_sequence_destroy(MORSE_context_t *morse, MORSE_sequence_t *sequence)
{
    QUARK_Sequence_Destroy((Quark*)(morse->schedopt), (Quark_Sequence *)(sequence->schedopt));
    return MORSE_SUCCESS;
}

/** *****************************************************************************
 *  Wait for the completion of a sequence
 **/
int RUNTIME_sequence_wait(MORSE_context_t *morse, MORSE_sequence_t *sequence )
{
    QUARK_Sequence_Wait( (Quark*)(morse->schedopt), (Quark_Sequence *)(sequence->schedopt));
    return MORSE_SUCCESS;
}

/** *****************************************************************************
 *  Terminate a sequence
 **/
void RUNTIME_sequence_flush(void *quark, MORSE_sequence_t *sequence, MORSE_request_t *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    QUARK_Sequence_Cancel((Quark*) quark, (Quark_Sequence *)(sequence->schedopt));
}
