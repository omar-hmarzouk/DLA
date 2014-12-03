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
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  MORSE_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "MORSE_zlange"
#define _FMULS FMULS_LANGE(M, N)
#define _FADDS FADDS_LANGE(M, N)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
	int hres = 0;
	int n;
	double normmorse, normlapack, result;
	double eps;
	int   norm[4]   = { MorseMaxNorm, MorseOneNorm, MorseInfNorm, MorseFrobeniusNorm };
	char *normstr[4]  = { "Max", "One", "Inf", "Fro" };
    PASTE_CODE_IPARAM_LOCALS( iparam );

    eps = LAPACKE_dlamch_work('e');

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, MORSE_Complex64_t, M, N );

    MORSE_zplrnt( M, N, A, LDA, 3436 );

    /* MORSE ZLANGE */
    START_TIMING();
    normmorse = MORSE_zlange(MorseInfNorm, M, N, A, LDA);
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        double *work = (double*) malloc(max(M,N)*sizeof(double));
    	normlapack = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseInfNorm), M, N, A, LDA, work);
    	result = fabs(normmorse - normlapack);
        switch(norm[2]) {
        case MorseMaxNorm:
            /* result should be perfectly equal */
            break;
        case MorseInfNorm:
            /* Sum order on the line can differ */
            result = result / (double)N;
            break;
        case MorseOneNorm:
            /* Sum order on the column can differ */
            result = result / (double)M;
            break;
        case MorseFrobeniusNorm:
            /* Sum oreder on every element can differ */
            result = result / ((double)M * (double)N);
            break;
        }
        if ( MORSE_My_Mpi_Rank() == 0 ) {
        	dparam[IPARAM_ANORM] = normlapack;
            dparam[IPARAM_BNORM] = 0.;
            dparam[IPARAM_XNORM] = 1.;
            dparam[IPARAM_RES] = result;
        }
        free( work );
    }

    free( A );

    return hres;
}
