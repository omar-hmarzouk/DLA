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
 * @file runtime_control.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version
 * @author Vijay Joshi
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "runtime/quark/include/morse_quark.h"

/*******************************************************************************
 * Thread rank.
 **/
int RUNTIME_rank(MORSE_context_t *morse)
{
    return QUARK_Thread_Rank((Quark*)(morse->schedopt));
}

/*******************************************************************************
 * Initialize MORSE
 **/
int RUNTIME_init_scheduler(MORSE_context_t *morse, int nworkers, int ncudas, int nthreads_per_worker)
{
    int hres = 0;
    if ( ncudas > 0 )
        morse_warning( "RUNTIME_init_scheduler(quark)", "GPUs are not supported for now");

    if ( nthreads_per_worker > 0 )
        morse_warning( "RUNTIME_init_scheduler(quark)", "Multi-threaded kernels are not supported for now");

    morse->schedopt = (void*)QUARK_New(nworkers);

    return hres;
}

/*******************************************************************************
 * Barrier MORSE.
 **/
void RUNTIME_barrier(MORSE_context_t *morse)
{
    QUARK_Barrier((Quark*)(morse->schedopt));
}

/*******************************************************************************
 * Finalize MORSE
 */
void RUNTIME_finalize_scheduler(MORSE_context_t *morse)
{
    QUARK_Delete((Quark*)(morse->schedopt));
    return;
}

/*******************************************************************************
 *  To suspend the processing of new tasks by workers
 **/
void RUNTIME_pause( MORSE_context_t *morse )
{
    (void)morse;
    return;
}

/*******************************************************************************
 *  This is the symmetrical call to RUNTIME_pause,
 *  used to resume the workers polling for new tasks.
 **/
void RUNTIME_resume( MORSE_context_t *morse )
{
    (void)morse;
    return;
}
