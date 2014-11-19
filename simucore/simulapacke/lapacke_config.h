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
 * @file lapacke_config.h
 *
 * MORSE lapacke header for simulation mode
 *
 * @version 1.0.0
 * @comment This is a copy of lapacke_config.h from lapack 3.5.0
 * @author Florent Pruvost
 * @date 2014-10-01
 *
 **/

#ifndef _LAPACKE_CONFIG_H_
#define _LAPACKE_CONFIG_H_

#ifdef __cplusplus
#if defined(LAPACK_COMPLEX_CPP)
#include <complex>
#endif
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>

#ifndef lapack_int
#if defined(LAPACK_ILP64)
#define lapack_int              long
#else
#define lapack_int              int
#endif
#endif

#ifndef lapack_logical
#define lapack_logical          lapack_int
#endif

#ifndef LAPACK_COMPLEX_CUSTOM

#if defined(LAPACK_COMPLEX_STRUCTURE)

typedef struct { float real, imag; } _lapack_complex_float;
typedef struct { double real, imag; } _lapack_complex_double;
#define lapack_complex_float  _lapack_complex_float
#define lapack_complex_double _lapack_complex_double
#define lapack_complex_float_real(z)  ((z).real)
#define lapack_complex_float_imag(z)  ((z).imag)
#define lapack_complex_double_real(z)  ((z).real)
#define lapack_complex_double_imag(z)  ((z).imag)

#elif defined(LAPACK_COMPLEX_C99)

#include <complex.h>
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
#define lapack_complex_float_real(z)       (creal(z))
#define lapack_complex_float_imag(z)       (cimag(z))
#define lapack_complex_double_real(z)       (creal(z))
#define lapack_complex_double_imag(z)       (cimag(z))

#elif defined(LAPACK_COMPLEX_CPP)

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#define lapack_complex_float_real(z)       ((z).real())
#define lapack_complex_float_imag(z)       ((z).imag())
#define lapack_complex_double_real(z)       ((z).real())
#define lapack_complex_double_imag(z)       ((z).imag())

#else

#include <complex.h>
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
#define lapack_complex_float_real(z)       (creal(z))
#define lapack_complex_float_imag(z)       (cimag(z))
#define lapack_complex_double_real(z)       (creal(z))
#define lapack_complex_double_imag(z)       (cimag(z))

#endif

lapack_complex_float lapack_make_complex_float( float re, float im );
lapack_complex_double lapack_make_complex_double( double re, double im );

#endif

#ifndef LAPACK_malloc
#define LAPACK_malloc( size )   malloc( size )
#endif

#ifndef LAPACK_free
#define LAPACK_free( p )        free( p )
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LAPACKE_CONFIG_H_ */
