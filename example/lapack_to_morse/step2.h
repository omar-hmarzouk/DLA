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
 * @file step2.h
 *
 *  MORSE example routines
 *  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
 *  University of Bordeaux, Bordeaux INP
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-29
 *
 **/

#ifndef STEP2_H
#define STEP2_H

/* Common include for all steps of the tutorial */
#include "lapack_to_morse.h"

/* Specific includes for step 2 */
#include <coreblas/include/lapacke.h>
#include "morse.h"

/* Integer parameters for step2 */
enum iparam_step2 {
    IPARAM_THRDNBR,        /* Number of cores                            */
    IPARAM_N,              /* Number of columns of the matrix            */
    IPARAM_NRHS,           /* Number of RHS                              */
    /* End */
    IPARAM_SIZEOF
};

/* Specific routines used in step2.c main program */

/******************************************************************************
 * Initialize integer parameters
 */
static void init_iparam(int iparam[IPARAM_SIZEOF]){
    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_NRHS          ] = 1;
 }

/******************************************************************************
 * Print how to use the program
 */
static void show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --n=X            dimension (N). (default: 500)\n"
            "  --nrhs=X         number of RHS. (default: 1)\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "\n");
}

/******************************************************************************
 * Read arguments following step2 program call
 */
static void read_args(int argc, char *argv[], int *iparam){
    int i;
    for (i = 1; i < argc && argv[i]; ++i) {
        if ( startswith( argv[i], "--help") || startswith( argv[i], "-help") ||
             startswith( argv[i], "--h") || startswith( argv[i], "-h") ) {
            show_help( argv[0] );
            exit(0);
        } else if (startswith( argv[i], "--n=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_N]) );
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NRHS]) );
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }
}

/******************************************************************************
 * Print a header message to summarize main parameters
 */
static void print_header(char *prog_name, int * iparam) {
#if defined(CHAMELEON_SIMULATION)
    double    eps = 0.;
#else
    double    eps = LAPACKE_dlamch_work( 'e' );
#endif

    printf( "#\n"
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb gpus:    %d\n"
            "# N:          %d\n"
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            0,
            iparam[IPARAM_N],
            128,
            32,
            eps );

    printf( "#      M       N  K/NRHS   seconds   Gflop/s\n");
    printf( "#%7d %7d %7d ", iparam[IPARAM_N], iparam[IPARAM_N], iparam[IPARAM_NRHS]);
    fflush( stdout );
    return;
}

/******************************************************************************
 * Macro to allocate a matrix as a 1D array
 */
#define PASTE_CODE_ALLOCATE_MATRIX(_name_, _type_, _m_, _n_)    \
    _type_ *_name_ = NULL;                                      \
    _name_ = (_type_*)malloc( (_m_) * (_n_) * sizeof(_type_) ); \
    if ( ! _name_ ) {                                           \
        fprintf(stderr, "Out of Memory for %s\n", #_name_);     \
        return -1;                                              \
    }

#endif /* STEP2_H */
