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
 * @file runtime_descriptor.c
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
#include <stdlib.h>
#include "runtime/quark/include/morse_quark.h"

void RUNTIME_desc_init( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

void RUNTIME_desc_create( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

void RUNTIME_desc_destroy( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

void RUNTIME_desc_submatrix( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

int RUNTIME_desc_acquire( MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

int RUNTIME_desc_release( MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

int RUNTIME_desc_getoncpu( MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

void *RUNTIME_desc_getaddr( MORSE_desc_t *desc, int m, int n )
{
    return desc->get_blkaddr( desc, m, n );
}
