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
 * @file auxiliary.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/
#ifndef _MORSE_AUXILIARY_H_
#define _MORSE_AUXILIARY_H_

#include "morse_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 *  Internal routines
 **/
void morse_warning      (const char *func_name, char* msg_text);
void morse_error        (const char *func_name, char* msg_text);
void morse_fatal_error  (const char *func_name, char* msg_text);
int  morse_rank         (MORSE_context_t *morse);
int  morse_tune         (MORSE_enum func, int M, int N, int NRHS);

#ifdef __cplusplus
}
#endif

#endif
