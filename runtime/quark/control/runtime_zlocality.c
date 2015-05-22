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
 *  @file runtime_zlocality.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version
 *  @author Vijay Joshi
 *  @date 2011-10-29
 *  @precisions normal z -> s d c
 *
 **/
#include "runtime/quark/include/morse_quark.h"

void RUNTIME_zlocality_allrestrict( uint32_t where )
{
    (void)where;
    morse_warning("RUNTIME_zlocality_allrestrict(quark)", "Kernel locality cannot be specified with Quark");
}

void RUNTIME_zlocality_onerestrict( MORSE_kernel_t kernel, uint32_t where )
{
    (void)kernel;
    (void)where;
    morse_warning("RUNTIME_zlocality_onerestrict(quark)", "Kernel locality cannot be specified with Quark");
}

void RUNTIME_zlocality_allrestore( )
{
    morse_warning("RUNTIME_zlocality_allrestore(quark)", "Kernel locality cannot be specified with Quark");
}

void RUNTIME_zlocality_onerestore( MORSE_kernel_t kernel )
{
    (void)kernel;
    morse_warning("RUNTIME_zlocality_onerestore(quark)", "Kernel locality cannot be specified with Quark");
}
