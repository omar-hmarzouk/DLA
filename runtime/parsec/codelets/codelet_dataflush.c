/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 **/

#include "runtime/parsec/include/morse_parsec.h"

void MORSE_TASK_dataflush(const MORSE_option_t *options,
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
