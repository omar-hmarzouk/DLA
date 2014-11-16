/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University 
 *                          of Tennessee Research Foundation. 
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 IPB. All rights reserved. 
 *
 **/

/**
 *
 * @file auxiliary.c
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
 * @date 2012-09-15
 *
 **/
#include "common.h"
#include "auxiliary.h"

#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action without severe consequences.
 *  Problems occuring while MORSE is being used correctly.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 ******************************************************************************/
void morse_warning(const char *func_name, char* msg_text)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL)
        morse_fatal_error("morse_warning", "MORSE not initialized");
    if (morse->warnings_enabled)
        fprintf(stderr, "MORSE WARNING: %s(): %s\n", func_name, msg_text);
}

/*******************************************************************************
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action with potentially severe consequences.
 *  Problems occuring due to incorrect use of MORSE.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 ******************************************************************************/
void morse_error(const char *func_name, char* msg_text)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL)
        morse_fatal_error("morse_error", "MORSE not initialized");
    if (morse->errors_enabled)
        fprintf(stderr, "MORSE ERROR: %s(): %s\n", func_name, msg_text);

}

/*******************************************************************************
 *
 *  Unexpected behavior within the library.
 *  Unrecoverable user errors.
 *  Context oblivious.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 ******************************************************************************/
void morse_fatal_error(const char *func_name, char* msg_text)
{
    fprintf(stderr, "MORSE FATAL ERROR: %s(): %s\n", func_name, msg_text);
    exit(0);
}

/*******************************************************************************
 *
 **/
int morse_element_size(int type)
{
    switch(type) {
        case MorseByte:          return          1;
        case MorseInteger:       return   sizeof(int);
        case MorseRealFloat:     return   sizeof(float);
        case MorseRealDouble:    return   sizeof(double);
        case MorseComplexFloat:  return 2*sizeof(float);
        case MorseComplexDouble: return 2*sizeof(double);
        default: morse_fatal_error("morse_element_size", "undefined type");
                 return MORSE_ERR_ILLEGAL_VALUE;

    }
}

/*******************************************************************************
 *  Returns core id
 **/
int morse_rank(MORSE_context_t *morse)
{
    return RUNTIME_rank( morse );
}

/*******************************************************************************
 *  Tune block size nb and internal block size ib
 **/
int morse_tune(MORSE_enum func, int M, int N, int NRHS)
{
    (void)func;
    (void)M;
    (void)N;
    (void)NRHS;
    morse_warning( "morse_tune", "Autotunning not available for now" );
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *
 * @ingroup Auxiliary
 *
 *  MORSE_Version - Reports MORSE version number.
 *
 *******************************************************************************
 *
 * @param[out] ver_major
 *          MORSE major version number.
 *
 * @param[out] ver_minor
 *          MORSE minor version number.
 *
 * @param[out] ver_micro
 *          MORSE micro version number.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int MORSE_Version(int *ver_major, int *ver_minor, int *ver_micro)
{
    if (! ver_major && ! ver_minor && ! ver_micro)
        return  MORSE_ERR_ILLEGAL_VALUE;

    if (ver_major)
        *ver_major = CHAMELEON_VERSION_MAJOR;

    if (ver_minor)
        *ver_minor = CHAMELEON_VERSION_MINOR;

    if (ver_micro)
        *ver_micro = CHAMELEON_VERSION_MICRO;

    return MORSE_SUCCESS;
}