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
 * @file codelet_dataflush.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 *
 * @date 2014-02-05
 *
 **/
#include "runtime/quark/include/morse_quark.h"

void MORSE_TASK_dataflush(MORSE_option_t *options,
                          MORSE_desc_t *A, int Am, int An)
{
    (void)options; (void)A;

    /*
     * This is useful for StarPU implementation, if it happens in Quark, it will
     * need to be done carefuly to not break both runtimes.
     */
}

void MORSE_TASK_dataflush_all()
{
}
