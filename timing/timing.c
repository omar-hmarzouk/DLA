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
 * @file timing.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Dulceneia Becker
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

#if defined( _WIN32 ) || defined( _WIN64 )
#define int64_t __int64
#endif

/* Define these so that the Microsoft VC compiler stops complaining
   about scanf and friends */
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#else  /* Non-Windows */
#include <unistd.h>
#include <sys/resource.h>
#endif

#include <coreblas/include/cblas.h>
#include <coreblas/include/lapacke.h>
#include <morse.h>
#include <coreblas/include/coreblas.h>
#include "flops.h"
#include "timing.h"
#include <control/auxiliary.h>

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

static int RunTest(int *iparam, _PREC *dparam, double *t_);

int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

static int
Test(int64_t n, int *iparam) {
    int      i, j, iter;
    int      thrdnbr, niter;
    int64_t  M, N, K, NRHS;
    double  *t;
#if defined(CHAMELEON_SIMULATION)
    _PREC    eps = 0.;
#else
    _PREC    eps = _LAMCH( 'e' );
#endif
    _PREC    dparam[IPARAM_DNBPARAM];
    double   fmuls, fadds, fp_per_mul, fp_per_add;
    double   sumgf, sumgf2, sumt, sd, flops, gflops;
    char    *s;
    char    *env[] = {
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "GOTO_NUM_THREADS",
        "ACML_NUM_THREADS",
        "ATLAS_NUM_THREADS",
        "BLAS_NUM_THREADS", ""
    };
    int gnuplot = 0;

/*
 * if hres = 0 then the test succeed
 * if hres = n then the test failed n times
 */
    int hres = 0;

    memset( &dparam, 0, IPARAM_DNBPARAM * sizeof(_PREC) );
    dparam[IPARAM_THRESHOLD_CHECK] = 100.0;

    thrdnbr = iparam[IPARAM_THRDNBR];
    niter   = iparam[IPARAM_NITER];

    M    = iparam[IPARAM_M];
    N    = iparam[IPARAM_N];
    K    = iparam[IPARAM_K];
    NRHS = K;
    (void)M;(void)N;(void)K;(void)NRHS;

    if ( (n < 0) || (thrdnbr < 0 ) ) {
        if (gnuplot && (MORSE_My_Mpi_Rank() == 0) ) {
            printf( "set title '%d_NUM_THREADS: ", thrdnbr );
            for (i = 0; env[i][0]; ++i) {
                s = getenv( env[i] );

                if (i) printf( " " ); /* separating space */

                for (j = 0; j < 5 && env[i][j] && env[i][j] != '_'; ++j)
                    printf( "%c", env[i][j] );

                if (s)
                    printf( "=%s", s );
                else
                    printf( "->%s", "?" );
            }
            printf( "'\n" );
            printf( "%s\n%s\n%s\n%s\n%s%s%s\n",
                    "set xlabel 'Matrix size'",
                    "set ylabel 'Gflop/s'",
                    "set key bottom",
                    gnuplot > 1 ? "set terminal png giant\nset output 'timeplot.png'" : "",
                    "plot '-' using 1:5 title '", _NAME, "' with linespoints" );
        }
        return 0;
    }

    if ( MORSE_My_Mpi_Rank() == 0)
        printf( "%7d %7d %7d ", iparam[IPARAM_M], iparam[IPARAM_N], iparam[IPARAM_K] );
    fflush( stdout );

    t = (double*)malloc(niter*sizeof(double));
    memset(t, 0, niter*sizeof(double));

    if (sizeof(_TYPE) == sizeof(_PREC)) {
        fp_per_mul = 1;
        fp_per_add = 1;
    } else {
        fp_per_mul = 6;
        fp_per_add = 2;
    }

    fadds = (double)(_FADDS);
    fmuls = (double)(_FMULS);
    flops = 1e-9 * (fmuls * fp_per_mul + fadds * fp_per_add);
    gflops = 0.0;

    if ( iparam[IPARAM_WARMUP] ) {
        RunTest( iparam, dparam, &(t[0]));
    }

    sumgf  = 0.0;
    double sumgf_upper  = 0.0;
    sumgf2 = 0.0;
    sumt   = 0.0;

    for (iter = 0; iter < niter; iter++)
    {
        if( iter == 0 ) {
          if ( iparam[IPARAM_TRACE] )
            iparam[IPARAM_TRACE] = 2;
          if ( iparam[IPARAM_DAG] )
            iparam[IPARAM_DAG] = 2;
          if ( iparam[IPARAM_PROFILE] )
            iparam[IPARAM_PROFILE] = 2;

          RunTest( iparam, dparam, &(t[iter]));

          iparam[IPARAM_TRACE] = 0;
          iparam[IPARAM_DAG] = 0;
          iparam[IPARAM_PROFILE] = 0;
        }
        else
            RunTest( iparam, dparam, &(t[iter]));

        gflops = flops / t[iter];


        double upper_gflops = 0.0;
        double tmin = 0.0;
        double integer_tmin = 0.0;
#if defined (CHAMELEON_SCHED_STARPU)
        if (iparam[IPARAM_BOUND])
        {
#if 0
            if (iparam[IPARAM_BOUNDDEPS]) {
                FILE *out = fopen("bounddeps.pl", "w");
                starpu_bound_print_lp(out);
                fclose(out);
                out = fopen("bound.dot", "w");
                starpu_bound_print_dot(out);
                fclose(out);
            } else {
                FILE *out = fopen("bound.pl", "w");
                starpu_bound_print_lp(out);
                fclose(out);
#else
            {
#endif
                starpu_bound_compute(&tmin, &integer_tmin, 0);
                upper_gflops  = (flops / (tmin / 1000.0));
                sumgf_upper += upper_gflops;
            }
        }
#endif
        sumt   += t[iter];
        sumgf  += gflops;
        sumgf2 += gflops*gflops;
    }

    gflops = sumgf / niter;
    sd = sqrt((sumgf2 - (sumgf*sumgf)/niter)/niter);

    if ( MORSE_My_Mpi_Rank() == 0) {
        printf( "%9.3f %9.2f +-%7.2f  ", sumt/niter, gflops, sd);

        if (iparam[IPARAM_BOUND] && !iparam[IPARAM_BOUNDDEPS])
            printf(" %9.2f",  sumgf_upper/niter);

        if ( iparam[IPARAM_PEAK] )
        {
            if (dparam[IPARAM_ESTIMATED_PEAK]<0.0f)
                printf("  n/a    n/a   ");
            else
                printf("  %5.2f%%  %9.2f ", 100.0f*(gflops/dparam[IPARAM_ESTIMATED_PEAK]), dparam[IPARAM_ESTIMATED_PEAK]);
        }

        if ( iparam[IPARAM_CHECK] ){
            hres = ( dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ) > dparam[IPARAM_THRESHOLD_CHECK] );

            if (hres)
                printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e FAILURE",
                    dparam[IPARAM_RES], dparam[IPARAM_ANORM], dparam[IPARAM_XNORM], dparam[IPARAM_BNORM],
                    dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ));
            else
                printf( "%8.5e %8.5e %8.5e %8.5e                       %8.5e SUCCESS",
                    dparam[IPARAM_RES], dparam[IPARAM_ANORM], dparam[IPARAM_XNORM], dparam[IPARAM_BNORM],
                    dparam[IPARAM_RES] / n / eps / (dparam[IPARAM_ANORM] * dparam[IPARAM_XNORM] + dparam[IPARAM_BNORM] ));
        }

        if ( iparam[IPARAM_INVERSE] )
            printf( " %8.5e %8.5e %8.5e     %8.5e",
                    dparam[IPARAM_RNORM], dparam[IPARAM_ANORM], dparam[IPARAM_AinvNORM],
                    dparam[IPARAM_RNORM] /((dparam[IPARAM_ANORM] * dparam[IPARAM_AinvNORM])*n*eps));

        printf("\n");

        fflush( stdout );
    }
    free(t);

    return hres;
}

static int
startswith(const char *s, const char *prefix) {
    size_t n = strlen( prefix );
    if (strncmp( s, prefix, n ))
        return 0;
    return 1;
}

static int
get_range(char *range, int *start_p, int *stop_p, int *step_p) {
    char *s, *s1, buf[21];
    int colon_count, copy_len, nbuf=20, n;
    int start=1000, stop=10000, step=1000;

    colon_count = 0;
    for (s = strchr( range, ':'); s; s = strchr( s+1, ':'))
        colon_count++;

    if (colon_count == 0) { /* No colon in range. */
        if (sscanf( range, "%d", &start ) < 1 || start < 1)
            return -1;
        step = start / 10;
        if (step < 1) step = 1;
        stop = start + 10 * step;

    } else if (colon_count == 1) { /* One colon in range.*/
        /* First, get the second number (after colon): the stop value. */
        s = strchr( range, ':' );
        if (sscanf( s+1, "%d", &stop ) < 1 || stop < 1)
            return -1;

        /* Next, get the first number (before colon): the start value. */
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start > stop || start < 1)
            return -1;

        /* Let's have 10 steps or less. */
        step = (stop - start) / 10;
        if (step < 1)
            step = 1;
    } else if (colon_count == 2) { /* Two colons in range. */
        /* First, get the first number (before the first colon): the start value. */
        s = strchr( range, ':' );
        n = s - range;
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, range, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &start ) < 1 || start < 1)
            return -1;

        /* Next, get the second number (after the first colon): the stop value. */
        s1 = strchr( s+1, ':' );
        n = s1 - (s + 1);
        copy_len = n > nbuf ? nbuf : n;
        strncpy( buf, s+1, copy_len );
        buf[copy_len] = 0;
        if (sscanf( buf, "%d", &stop ) < 1 || stop < start)
            return -1;

        /* Finally, get the third number (after the second colon): the step value. */
        if (sscanf( s1+1, "%d", &step ) < 1 || step < 1)
            return -1;
    } else

        return -1;

    *start_p = start;
    *stop_p = stop;
    *step_p = step;

    return 0;
}

static void
show_help(char *prog_name) {
    printf( "Usage:\n%s [options]\n\n", prog_name );
    printf( "Options are:\n"
            "  --help           Show this help\n"
            "\n"
            "  --threads=X      Number of CPU workers (default: _SC_NPROCESSORS_ONLN)\n"
            "  --gpus=X         Number of GPU workers (default: 0)\n"
            "\n"
            "  --[a]sync        Enable/Disable synchronous calls in wrapper function such as POTRI. (default: async)\n"
            "  --[no]check      Check result (default: nocheck)\n"
            "  --[no]inv        Check on inverse (default: noinv)\n"
            "  --[no]warmup     Perform a warmup run to pre-load libraries (default: warmup)\n"
            "  --[no]trace      Enable/Disable trace generation (default: notrace)\n"
            "  --[no]dag        Enable/Disable DAG generation (default: nodag)\n"
            "                   Generates a dot_dag_file.dot.\n"
             "  --[no]profile   Print profiling informations (default: noprofile)\n"
             "  --nocpu         All GPU kernels are exclusively executed on GPUs (default: 0)\n"
/*             "  --inplace        Enable layout conversion inplace for lapack interface timers (default: enable)\n" */
/*             "  --outplace       Enable layout conversion out of place for lapack interface timers (default: disable)\n" */
/*             "  --[no]atun       Activate autotuning (default: noatun)\n" */
            "\n"
            "  --n_range=R      Range of N values\n"
            "                   with R=Start:Stop:Step (default: 500:5000:500)\n"
            "  --m=X            dimension (M) of the matrices (default: N)\n"
            "  --k=X            dimension (K) of the matrices (default: 1)\n"
            "  --nrhs=X         Number of right-hand size (default: 1)\n"
            "  --nb=N           Nb size. (default: 128)\n"
            "  --ib=N           IB size. (default: 32)\n"
            "\n"
            "  --niter=N        Number of iterations performed for each test (default: 1)\n"
            "\n"
            "  --rhblk=N        If N > 0, enable Householder mode for QR and LQ factorization\n"
            "                   N is the size of each subdomain (default: 0)\n"
/*             "\n" */
/*             " Options specific to the conversion format timings xgetri and xgecfi:\n" */
/*             "  --ifmt           Input format. (default: 0)\n" */
/*             "  --ofmt           Output format. (default: 1)\n" */
/*             "                   The possible values are:\n" */
/*             "                     0 - MorseCM, Column major\n" */
/*             "                     1 - MorseCCRB, Column-Colum rectangular block\n" */
/*             "                     2 - MorseCRRB, Column-Row rectangular block\n" */
/*             "                     3 - MorseRCRB, Row-Colum rectangular block\n" */
/*             "                     4 - MorseRRRB, Row-Row rectangular block\n" */
/*             "                     5 - MorseRM, Row Major\n" */
/*             "  --thrdbypb       Number of threads per subproblem for inplace transformation (default: 1)\n" */
            "\n");
}


static void
print_header(char *prog_name, int * iparam) {
    const char *bound_header   = iparam[IPARAM_BOUND]   ? "   thGflop/s" : "";
    const char *check_header   = iparam[IPARAM_CHECK]   ? "     ||Ax-b||       ||A||       ||x||       ||b|| ||Ax-b||/N/eps/(||A||||x||+||b||)  RETURN" : "";
    const char *inverse_header = iparam[IPARAM_INVERSE] ? " ||I-A*Ainv||       ||A||    ||Ainv||       ||Id - A*Ainv||/((||A|| ||Ainv||).N.eps)" : "";
    const char *peak_header    = iparam[IPARAM_PEAK]    ? "  (% of peak)  peak" : "";
#if defined(CHAMELEON_SIMULATION)
    _PREC    eps = 0.;
#else
    _PREC    eps = _LAMCH( 'e' );
#endif

    printf( "#\n"
            "# CHAMELEON %d.%d.%d, %s\n"
            "# Nb threads: %d\n"
            "# Nb GPUs:    %d\n"
#if defined(CHAMELEON_USE_MPI)
            "# Nb mpi:     %d\n"
            "# PxQ:        %dx%d\n"
#endif
            "# NB:         %d\n"
            "# IB:         %d\n"
            "# eps:        %e\n"
            "#\n",
            CHAMELEON_VERSION_MAJOR,
            CHAMELEON_VERSION_MINOR,
            CHAMELEON_VERSION_MICRO,
            prog_name,
            iparam[IPARAM_THRDNBR],
            iparam[IPARAM_NCUDAS],
#if defined(CHAMELEON_USE_MPI)
            iparam[IPARAM_NMPI],
            iparam[IPARAM_P], iparam[IPARAM_Q],
#endif
            iparam[IPARAM_NB],
            iparam[IPARAM_IB],
            eps );

    printf( "#     M       N  K/NRHS   seconds   Gflop/s Deviation%s%s%s\n",
            bound_header, peak_header, iparam[IPARAM_INVERSE] ? inverse_header : check_header);
    return;
}

static void
get_thread_count(int *thrdnbr) {
#if defined WIN32 || defined WIN64
    sscanf( getenv( "NUMBER_OF_PROCESSORS" ), "%d", thrdnbr );
#else
    *thrdnbr = sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

int
main(int argc, char *argv[]) {
    int i, m, mx, nx;
    int nbnode = 1;
    int start =  500;
    int stop  = 5000;
    int step  =  500;
    int iparam[IPARAM_SIZEOF];
    int success = 0;

    memset(iparam, 0, IPARAM_SIZEOF*sizeof(int));

    iparam[IPARAM_THRDNBR       ] = -1;
    iparam[IPARAM_THRDNBR_SUBGRP] = 1;
    iparam[IPARAM_SCHEDULER     ] = 0;
    iparam[IPARAM_M             ] = -1;
    iparam[IPARAM_N             ] = 500;
    iparam[IPARAM_K             ] = 1;
    iparam[IPARAM_LDA           ] = -1;
    iparam[IPARAM_LDB           ] = -1;
    iparam[IPARAM_LDC           ] = -1;
    iparam[IPARAM_MB            ] = 128;
    iparam[IPARAM_NB            ] = 128;
    iparam[IPARAM_IB            ] = 32;
    iparam[IPARAM_NITER         ] = 1;
    iparam[IPARAM_WARMUP        ] = 1;
    iparam[IPARAM_CHECK         ] = 0;
    iparam[IPARAM_VERBOSE       ] = 0;
    iparam[IPARAM_AUTOTUNING    ] = 0;
    iparam[IPARAM_INPUTFMT      ] = 0;
    iparam[IPARAM_OUTPUTFMT     ] = 0;
    iparam[IPARAM_TRACE         ] = 0;
    iparam[IPARAM_DAG           ] = 0;
    iparam[IPARAM_ASYNC         ] = 1;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_RHBLK         ] = 0;
    iparam[IPARAM_MX            ] = -1;
    iparam[IPARAM_NX            ] = -1;
    iparam[IPARAM_RHBLK         ] = 0;
    iparam[IPARAM_INPLACE       ] = MORSE_OUTOFPLACE;

    iparam[IPARAM_INVERSE       ] = 0;
    iparam[IPARAM_NCUDAS        ] = 0;
    iparam[IPARAM_NMPI          ] = 1;
    iparam[IPARAM_P             ] = 1;
    iparam[IPARAM_Q             ] = 1;
    iparam[IPARAM_PROFILE       ] = 0;
    iparam[IPARAM_PRINT_ERRORS  ] = 0;
    iparam[IPARAM_PEAK          ] = 0;
    iparam[IPARAM_PARALLEL_TASKS] = 0;
    iparam[IPARAM_NO_CPU        ] = 0;
    iparam[IPARAM_BOUND         ] = 0;
    iparam[IPARAM_BOUNDDEPS     ] = 0;
    iparam[IPARAM_BOUNDDEPSPRIO ] = 0;

    for (i = 1; i < argc && argv[i]; ++i) {
        if ( startswith( argv[i], "--help") || startswith( argv[i], "-help") ||
             startswith( argv[i], "--h") || startswith( argv[i], "-h") ) {
            show_help( argv[0] );
            return EXIT_SUCCESS;
        } else if (startswith( argv[i], "--threads=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR]) );
        } else if (startswith( argv[i], "--gpus=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NCUDAS]) );
        } else if (startswith( argv[i], "--check" )) {
            iparam[IPARAM_CHECK] = 1;
        } else if (startswith( argv[i], "--nocheck" )) {
            iparam[IPARAM_CHECK] = 0;
        } else if (startswith( argv[i], "--inv" )) {
            iparam[IPARAM_INVERSE] = 1;
        } else if (startswith( argv[i], "--noinv" )) {
            iparam[IPARAM_INVERSE] = 0;
        } else if (startswith( argv[i], "--warmup" )) {
            iparam[IPARAM_WARMUP] = 1;
        } else if (startswith( argv[i], "--nowarmup" )) {
            iparam[IPARAM_WARMUP] = 0;
/*         } else if (startswith( argv[i], "--atun" )) { */
/*             iparam[IPARAM_AUTOTUNING] = 1; */
/*         } else if (startswith( argv[i], "--noatun" )) { */
/*             iparam[IPARAM_AUTOTUNING] = 0; */
        } else if (startswith( argv[i], "--trace" )) {
            iparam[IPARAM_TRACE] = 1;
        } else if (startswith( argv[i], "--notrace" )) {
            iparam[IPARAM_TRACE] = 0;
        } else if (startswith( argv[i], "--dag" )) {
            iparam[IPARAM_DAG] = 1;
        } else if (startswith( argv[i], "--nodag" )) {
            iparam[IPARAM_DAG] = 0;
        } else if (startswith( argv[i], "--sync" )) {
            iparam[IPARAM_ASYNC] = 0;
        } else if (startswith( argv[i], "--async" )) {
            iparam[IPARAM_ASYNC] = 1;
        } else if (startswith( argv[i], "--n_range=" )) {
            get_range( strchr( argv[i], '=' ) + 1, &start, &stop, &step );
        } else if (startswith( argv[i], "--m=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_M]) );
        } else if (startswith( argv[i], "--nb=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NB]) );
            iparam[IPARAM_MB] = iparam[IPARAM_NB];
        } else if (startswith( argv[i], "--nrhs=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_K]) );
        } else if (startswith( argv[i], "--k=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_K]) );
        } else if (startswith( argv[i], "--ib=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_IB]) );
        } else if (startswith( argv[i], "--niter=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &iparam[IPARAM_NITER] );
        } else if (startswith( argv[i], "--mx=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_MX]) );
        } else if (startswith( argv[i], "--nx=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_NX]) );
        } else if (startswith( argv[i], "--rhblk=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_RHBLK]) );
/*         } else if (startswith( argv[i], "--inplace" )) { */
/*             iparam[IPARAM_INPLACE] = MORSE_INPLACE; */
/*         } else if (startswith( argv[i], "--outplace" )) { */
/*             iparam[IPARAM_INPLACE] = MORSE_OUTOFPLACE; */
/*         } else if (startswith( argv[i], "--ifmt=" )) { */
/*             sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_INPUTFMT]) ); */
/*         } else if (startswith( argv[i], "--ofmt=" )) { */
/*             sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_OUTPUTFMT]) ); */
/*         } else if (startswith( argv[i], "--thrdbypb=" )) { */
/*             sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_THRDNBR_SUBGRP]) ); */
        } else if (startswith( argv[i], "--profile" )) {
            iparam[IPARAM_PROFILE] = 1;
        } else if (startswith( argv[i], "--peak" )) {
            iparam[IPARAM_PEAK] = 1;
        } else if (startswith( argv[i], "--noprofile" )) {
            iparam[IPARAM_PROFILE] = 0;
        } else if (startswith( argv[i], "--printerrors" )) {
            iparam[IPARAM_PRINT_ERRORS] = 1;
/*         } else if (startswith( argv[i], "--parallel=" )) { */
/*             sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_PARALLEL_TASKS]) ); */
/*         } else if (startswith( argv[i], "--noparallel" )) { */
/*             iparam[IPARAM_PARALLEL_TASKS] = 0; */
        } else if (startswith( argv[i], "--nocpu" )) {
            iparam[IPARAM_NO_CPU] = 1;
        } else if (startswith( argv[i], "--bounddepsprio" )) {
            iparam[IPARAM_BOUND] = 1;
            iparam[IPARAM_BOUNDDEPS] = 1;
            iparam[IPARAM_BOUNDDEPSPRIO] = 1;
        } else if (startswith( argv[i], "--bounddeps" )) {
            iparam[IPARAM_BOUND] = 1;
            iparam[IPARAM_BOUNDDEPS] = 1;
        } else if (startswith( argv[i], "--bound" )) {
            iparam[IPARAM_BOUND] = 1;
        } else if (startswith( argv[i], "--p=" )) {
            sscanf( strchr( argv[i], '=' ) + 1, "%d", &(iparam[IPARAM_P]) );
        } else {
            fprintf( stderr, "Unknown option: %s\n", argv[i] );
        }
    }

    if ( iparam[IPARAM_THRDNBR] == -1 ) {
      get_thread_count( &(iparam[IPARAM_THRDNBR]) );
      iparam[IPARAM_THRDNBR] -= iparam[IPARAM_NCUDAS];
    }

    m  = iparam[IPARAM_M];
    mx = iparam[IPARAM_MX];
    nx = iparam[IPARAM_NX];

    /* Initialize MORSE */
    MORSE_Init( iparam[IPARAM_THRDNBR],
                iparam[IPARAM_NCUDAS] );

    MORSE_Disable(MORSE_AUTOTUNING);
    MORSE_Set(MORSE_TILE_SIZE,        iparam[IPARAM_NB] );
    MORSE_Set(MORSE_INNER_BLOCK_SIZE, iparam[IPARAM_IB] );

    /* Householder mode */
    if (iparam[IPARAM_RHBLK] < 1) {
        MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_FLAT_HOUSEHOLDER);
    } else {
        MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_TREE_HOUSEHOLDER);
        MORSE_Set(MORSE_HOUSEHOLDER_SIZE, iparam[IPARAM_RHBLK]);
    }

    if (iparam[IPARAM_PROFILE] == 1)
        MORSE_Enable(MORSE_PROFILING_MODE);

    if (iparam[IPARAM_PRINT_ERRORS] == 1)
        MORSE_Enable(MORSE_ERRORS);

#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_size( MPI_COMM_WORLD, &nbnode );
    iparam[IPARAM_NMPI] = nbnode;
    /* Check P */
    if ( (iparam[IPARAM_P] > 1) &&
         (nbnode % iparam[IPARAM_P] != 0) ) {
      fprintf(stderr, "ERROR: %d doesn't divide the number of node %d\n",
              iparam[IPARAM_P], nbnode );
      return EXIT_FAILURE;
    }
#endif
    iparam[IPARAM_Q] = nbnode / iparam[IPARAM_P];

    /* Layout conversion */
    MORSE_Set(MORSE_TRANSLATION_MODE, iparam[IPARAM_INPLACE]);

    if ( MORSE_My_Mpi_Rank() == 0 )
        print_header( argv[0], iparam);

    if (step < 1) step = 1;

    Test( -1, iparam ); /* print header */
    for (i = start; i <= stop; i += step)
    {
        if ( nx > 0 ) {
            iparam[IPARAM_M] = i;
            iparam[IPARAM_N] = max(1, i/nx);
        } else if ( mx > 0 ) {
            iparam[IPARAM_M] = max(1, i/mx);
            iparam[IPARAM_N] = i;
        } else {
            if ( m == -1 )
                iparam[IPARAM_M] = i;
            iparam[IPARAM_N] = i;
        }
        success += Test( iparam[IPARAM_N], iparam );
    }

    MORSE_Finalize();
    return success;
}

