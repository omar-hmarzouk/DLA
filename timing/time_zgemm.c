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

#define _NAME  "MORSE_zgemm"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEMM(M, N, K)
#define _FADDS FADDS_GEMM(M, N, K)

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, morse_time_t *t_) 
{
    MORSE_Complex64_t alpha, beta;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    LDB = max(K, iparam[IPARAM_LDB]);
    LDC = max(M, iparam[IPARAM_LDC]);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A,      1, MORSE_Complex64_t, LDA, K   );
    PASTE_CODE_ALLOCATE_MATRIX( B,      1, MORSE_Complex64_t, LDB, N   );
    PASTE_CODE_ALLOCATE_MATRIX( C,      1, MORSE_Complex64_t, LDC, N   );
    PASTE_CODE_ALLOCATE_MATRIX( C2, check, MORSE_Complex64_t, LDC, N   );

    MORSE_zplrnt( M, K, A, LDA,  453 );
    MORSE_zplrnt( K, N, B, LDB, 5673 );
    MORSE_zplrnt( M, N, C, LDC,  740 );

    LAPACKE_zlarnv_work(1, ISEED, 1, &alpha);
    LAPACKE_zlarnv_work(1, ISEED, 1, &beta );

    if (check)
    {
        memcpy(C2, C, LDC*N*sizeof(MORSE_Complex64_t));
    }

    START_TIMING();
    MORSE_zgemm( MorseNoTrans, MorseNoTrans, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
//        dparam[IPARAM_RES] = z_check_gemm( MorseNoTrans, MorseNoTrans, M, N, K,
//                                           &alpha, A, LDA, B, LDB, &beta, C, C2, LDC,
//                                           &(dparam[IPARAM_ANORM]),
//                                           &(dparam[IPARAM_BNORM]),
//                                           &(dparam[IPARAM_XNORM]));

        MORSE_Complex64_t beta_const = -1.0;
        double *work = (double *)malloc(max(K,max(M, N))* sizeof(double));

        dparam[IPARAM_ANORM] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseInfNorm), M, N, C2,   LDC, work);
        dparam[IPARAM_BNORM] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseInfNorm), M, N, C, LDC, work);

        cblas_zgemm(CblasColMajor, (CBLAS_TRANSPOSE)MorseNoTrans, (CBLAS_TRANSPOSE)MorseNoTrans, M, N, K,
                    CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), C2, LDC);

        dparam[IPARAM_XNORM] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseInfNorm), M, N, C2, LDC, work);

        cblas_zaxpy(LDC * N, CBLAS_SADDR(beta_const), C, 1, C2, 1);

        dparam[IPARAM_RES] = LAPACKE_zlange_work(LAPACK_COL_MAJOR, morse_lapack_const(MorseInfNorm), M, N, C2, LDC, work);

        free(work);
        free(C2);
    }

    free( A );
    free( B );
    free( C );

    return 0;
}
