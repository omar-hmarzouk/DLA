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
 * @file runtime_context.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2011-10-29 
 *
 **/
#include <stdlib.h>
#include "runtime/quark/include/morse_quark.h"

/*******************************************************************************
 *  Create new context
 **/
void RUNTIME_context_create(MORSE_context_t *morse)
{
    morse->scheduler = RUNTIME_SCHED_QUARK;
    /* Will require the static initialization if we want to use it in this code */
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
            fprintf(stderr, "Profiling is not available with Quark\n");
            break;
        case MORSE_BOUND:
            fprintf(stderr, "Bound computation is not available with Quark\n");
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
            fprintf(stderr, "Profiling is not available with Quark\n");
            break;
        case MORSE_BOUND:
            fprintf(stderr, "Bound computation is not available with Quark\n");
            break;
        default:
            return;
    }
    return;
}
