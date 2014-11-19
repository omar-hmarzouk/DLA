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
 * @file morse_simulate.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-02
 *
 **/
#ifndef _MORSE_SIMULATE_H_
#define _MORSE_SIMULATE_H_

/* we need this when starpu is compiled with simgrid enabled */
#if defined(CHAMELEON_SCHED_STARPU) && defined(CHAMELEON_SIMULATION)
#include <starpu_simgrid_wrap.h>
#endif

#endif
