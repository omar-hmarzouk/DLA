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
 *  @file runtime_zprofiling.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 0.9.0
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @author Cedric Castagnede
 *  @date 2011-06-01
 *  @precisions normal z -> s d c
 *
 **/
#include "runtime/quark/include/morse_quark.h"

void RUNTIME_zdisplay_allprofile()
{
    morse_warning("RUNTIME_zdisplay_allprofile(quark)", "Profiling is not available with Quark");
}

void RUNTIME_zdisplay_oneprofile( MORSE_kernel_t kernel )
{
    (void)kernel;
    morse_warning("RUNTIME_zdisplay_oneprofile(quark)", "Profiling is not available with Quark\n");
}

