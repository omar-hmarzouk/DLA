/**
 *
 * @file chameleon_parsec.h
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 **/
#ifndef _MORSE_PARSEC_H_
#define _MORSE_PARSEC_H_

#include <parsec.h>
#include <parsec/insert_function.h>

#include "control/common.h"

struct morse_parsec_desc_s;
typedef struct morse_parsec_desc_s morse_parsec_desc_t;

int morse_parsec_get_arena_index(const MORSE_desc_t *desc);

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( parsec_dtd_tile_of( (parsec_data_collection_t *) ((desc)->schedopt), \
                                                            ((parsec_data_collection_t *) (desc)->schedopt)->data_key((desc)->schedopt, m, n) ))

#define RUNTIME_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION

#endif /* _MORSE_PARSEC_H_ */
