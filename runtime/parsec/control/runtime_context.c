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
 *  Create new context
 **/
void RUNTIME_context_create(MORSE_context_t *morse)
{
    /* In case of PaRSEC, this is done in init */
    morse->scheduler = RUNTIME_SCHED_PARSEC;
    return;
}

/*******************************************************************************
 *  Clean the context
 **/

void RUNTIME_context_destroy(MORSE_context_t *morse)
{
    (void)morse;
    return;
}

/*******************************************************************************
 *
 */
void RUNTIME_enable(MORSE_enum lever)
{
    switch (lever)
    {
        case MORSE_PROFILING_MODE:
            // TODO: check correctly for this
            //dague_profiling_start();
            break;
        default:
            return;
    }
    return;
}

/*******************************************************************************
 *
 **/
void RUNTIME_disable(MORSE_enum lever)
{
    switch (lever)
    {
        case MORSE_PROFILING_MODE:
            // TODO: check correctly for this
            //dague_profiling_stop();
            break;
        default:
            return;
    }
    return;
}
