/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
#include "runtime/PaRSEC/include/chameleon_parsec.h"

void RUNTIME_zlocality_allrestrict( uint32_t where )
{
    (void)where;
    morse_warning("RUNTIME_zlocality_allrestrict(PaRSEC)", "Kernel locality cannot be specified with PaRSEC");
}

void RUNTIME_zlocality_onerestrict( MORSE_kernel_t kernel, uint32_t where )
{
    (void)kernel;
    (void)where;
    morse_warning("RUNTIME_zlocality_onerestrict(PaRSEC)", "Kernel locality cannot be specified with PaRSEC");
}

void RUNTIME_zlocality_allrestore( )
{
    morse_warning("RUNTIME_zlocality_allrestore(PaRSEC)", "Kernel locality cannot be specified with PaRSEC");
}

void RUNTIME_zlocality_onerestore( MORSE_kernel_t kernel )
{
    (void)kernel;
    morse_warning("RUNTIME_zlocality_onerestore(PaRSEC)", "Kernel locality cannot be specified with PaRSEC");
}
