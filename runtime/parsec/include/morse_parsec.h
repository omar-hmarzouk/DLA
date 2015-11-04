/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

#ifndef _MORSE_PARSEC_H_
#define _MORSE_PARSEC_H_

#include <dague.h>
#include <dague/insert_function.h>
#include <dague/data_internal.h>
#include <dague/profiling.h>

#include "control/common.h"

struct morse_parsec_desc_s;
typedef struct morse_parsec_desc_s morse_parsec_desc_t;

/*
 * Access to block pointer and leading dimension
 */
//#define RTBLKADDR( desc, type, m, n ) ( (type*)RUNTIME_desc_getaddr( desc, m, n ) )
#define RTBLKADDR( desc, type, m, n ) ( tile_manage( DAGUE_dtd_handle, (desc)->schedopt, m, n ) )

#endif /* _MORSE_PARSEC_H_ */
