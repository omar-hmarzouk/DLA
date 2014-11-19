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
 * @file morse_types.h
 *
 *  MORSE header
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 **/

#ifndef _MORSE_TYPES_H_
#define _MORSE_TYPES_H_

/** ****************************************************************************
 *  System requirements
 **/
#include <stddef.h>
#if defined( _WIN32 )
  /* This must be included before INPUT is defined below, otherwise we
     have a name clash/problem  */
  #include <windows.h>
  #include <limits.h>
#else /* _WIN32 */
  #include <inttypes.h>
#endif /* _WIN32 */


/** ****************************************************************************
 *  MORSE types
 **/
typedef int  MORSE_enum;
typedef int  MORSE_bool;
typedef long MORSE_index;
typedef long MORSE_size;


/** ****************************************************************************
 * MORSE Complex numbers
 **/
#define MORSE_HAS_COMPLEX_H 1
#if defined(_WIN32)
# include <float.h>
# if defined(__INTEL_COMPILER)
    /* Fix name conflict within the cabs prototype (_Complex) that    */
    /* conflicts with a C99 keyword.                                  */
    #define _Complex __ConflictingComplex
    #include <math.h>
    #undef _Complex
    #undef complex
    typedef float  _Complex MORSE_Complex32_t;
    typedef double _Complex MORSE_Complex64_t;
# else /* __INTEL_COMPILER */
    /* Use MS VC complex class                                        */
    #include <complex>
    typedef std::complex<float> MORSE_Complex32_t;
    typedef std::complex<double> MORSE_Complex64_t;
    /* For LAPACKE lapacke.h force usage of Windows C++ Complex types */
    #define LAPACK_COMPLEX_CUSTOM
    #define lapack_complex_float std::complex<float>
    #define lapack_complex_double std::complex<double>
    #undef MORSE_HAS_COMPLEX_H
# endif /* __INTEL_COMPILER */
# define isnan _isnan
# define isinf !_finite
#else /* _WIN32 */
    typedef float  _Complex MORSE_Complex32_t;
    typedef double _Complex MORSE_Complex64_t;
#endif /* _WIN32 */


/* Sun doesn't ship the complex.h header. Sun Studio doesn't have it and older GCC compilers don't have it either. */
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC) || defined(sun) || defined(__sun)
#undef MORSE_HAS_COMPLEX_H
#endif /* __SUNPRO_C */

#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define MORSE_DEPRECATED  __attribute__((__deprecated__))
#else
#define MORSE_DEPRECATED
#endif /* __GNUC__ */


#ifdef MORSE_HAS_COMPLEX_H
#include <complex.h>
#else /* MORSE_HAS_COMPLEX_H*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* These declarations will not clash with what C++ provides because the names in C++ are name-mangled. */
#if !defined(_WIN32)
extern double cabs(MORSE_Complex64_t z);
extern MORSE_Complex64_t conj(MORSE_Complex64_t z);
#endif /* _WIN32 */
extern float cabsf(MORSE_Complex32_t z);
extern double cimag(MORSE_Complex64_t z);
extern double creal(MORSE_Complex64_t z);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MORSE_HAS_COMPLEX_H*/

#endif /* __CHAMELEON_H__ */
