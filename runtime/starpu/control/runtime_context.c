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
 * @version 0.9.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "runtime/starpu/include/morse_starpu.h"

/*******************************************************************************
 *  Create new context
 **/
void RUNTIME_context_create( MORSE_context_t *morse )
{
    starpu_conf_t *conf;

    morse->scheduler = CHAMELEON_SCHED_STARPU;
    morse->schedopt = (void*) malloc (sizeof(struct starpu_conf));
    conf = morse->schedopt;

    starpu_conf_init( conf );

    return;
}

/*******************************************************************************
 *  Clean the context
 **/

void RUNTIME_context_destroy( MORSE_context_t *morse )
{
    free(morse->schedopt);
    return;
}

/*******************************************************************************
 *
 */
void RUNTIME_enable( MORSE_enum lever )
{
    switch (lever)
    {
        case MORSE_PROFILING_MODE:
            starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
            break;
        case MORSE_BOUND:
            starpu_bound_start(0, 0);
            break;
        default:
            return;
    }
    return;
}

/*******************************************************************************
 *
 **/
void RUNTIME_disable( MORSE_enum lever )
{
    switch (lever)
    {
        case MORSE_PROFILING_MODE:
            starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
            break;
        case MORSE_BOUND:
            starpu_bound_stop();
            break;
        default:
            return;
    }
    return;
}
