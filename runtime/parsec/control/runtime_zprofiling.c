/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
#include "chameleon_parsec.h"

void RUNTIME_zdisplay_allprofile()
{
    morse_warning("RUNTIME_zdisplay_allprofile(PaRSEC)", "Profiling is not available with PaRSEC");
}

void RUNTIME_zdisplay_oneprofile( MORSE_kernel_t kernel )
{
    (void)kernel;
    morse_warning("RUNTIME_zdisplay_oneprofile(PaRSEC)", "Profiling is not available with PaRSEC\n");
}

