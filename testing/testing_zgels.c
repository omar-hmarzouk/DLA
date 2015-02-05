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
 * @file testing_zgels.c
 *
 *  MORSE testing routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Bilel Hadri
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <morse.h>
#include <coreblas/include/cblas.h>
#include <coreblas/include/lapacke.h>
#include <coreblas/include/coreblas.h>
#include "testing_zauxiliary.h"

#undef REAL
#define COMPLEX

enum blas_order_type {
            blas_rowmajor = 101,
            blas_colmajor = 102 };

enum blas_uplo_type  {
            blas_upper = 121,
            blas_lower = 122 };

enum blas_cmach_type {
            blas_base      = 151,
            blas_t         = 152,
            blas_rnd       = 153,
            blas_ieee      = 154,
            blas_emin      = 155,
            blas_emax      = 156,
            blas_eps       = 157,
            blas_prec      = 158,
            blas_underflow = 159,
            blas_overflow  = 160,
            blas_sfmin     = 161};

enum blas_norm_type {
            blas_one_norm       = 171,
            blas_real_one_norm  = 172,
            blas_two_norm       = 173,
            blas_frobenius_norm = 174,
            blas_inf_norm       = 175,
            blas_real_inf_norm  = 176,
            blas_max_norm       = 177,
            blas_real_max_norm  = 178 };

static void
BLAS_error(char *rname, int err, int val, int x) {
  fprintf( stderr, "%s %d %d %d\n", rname, err, val, x );
  abort();
}

static
void
BLAS_zsy_norm(enum blas_order_type order, enum blas_norm_type norm,
  enum blas_uplo_type uplo, int n, const MORSE_Complex64_t *a, int lda, double *res) {
  int i, j; double anorm, v;
  char rname[] = "BLAS_zsy_norm";

  if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

  if (norm == blas_inf_norm) {
    anorm = 0.0;
    if (blas_upper == uplo) {
      for (i = 0; i < n; ++i) {
        v = 0.0;
        for (j = 0; j < i; ++j) {
          v += cabs( a[j + i * lda] );
        }
        for (j = i; j < n; ++j) {
          v += cabs( a[i + j * lda] );
        }
        if (v > anorm)
          anorm = v;
      }
    } else {
      BLAS_error( rname, -3, norm, 0 );
      return;
    }
  } else {
    BLAS_error( rname, -2, norm, 0 );
    return;
  }

  if (res) *res = anorm;
}

static
void
BLAS_zge_norm(enum blas_order_type order, enum blas_norm_type norm,
  int m, int n, const MORSE_Complex64_t *a, int lda, double *res) {
  int i, j; float anorm, v;
  char rname[] = "BLAS_zge_norm";

  if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

  if (norm == blas_frobenius_norm) {
    anorm = 0.0f;
    for (j = n; j; --j) {
      for (i = m; i; --i) {
        v = a[0];
        anorm += v * v;
        a++;
      }
      a += lda - m;
    }
    anorm = sqrt( anorm );
  } else if (norm == blas_inf_norm) {
    anorm = 0.0f;
    for (i = 0; i < m; ++i) {
      v = 0.0f;
      for (j = 0; j < n; ++j) {
        v += cabs( a[i + j * lda] );
      }
      if (v > anorm)
        anorm = v;
    }
  } else {
    BLAS_error( rname, -2, norm, 0 );
    return;
  }

  if (res) *res = anorm;
}

static
double
BLAS_dpow_di(double x, int n) {
  double rv = 1.0;

  if (n < 0) {
    n = -n;
    x = 1.0 / x;
  }

  for (; n; n >>= 1, x *= x) {
    if (n & 1)
      rv *= x;
  }

  return rv;
}

static
double
BLAS_dfpinfo(enum blas_cmach_type cmach) {
  double eps = 1.0, r = 1.0, o = 1.0, b = 2.0;
  int t = 53, l = 1024, m = -1021;
  char rname[] = "BLAS_dfpinfo";

  if ((sizeof eps) == sizeof(float)) {
    t = 24;
    l = 128;
    m = -125;
  } else {
    t = 53;
    l = 1024;
    m = -1021;
  }

  /* for (i = 0; i < t; ++i) eps *= half; */
  eps = BLAS_dpow_di( b, -t );
  /* for (i = 0; i >= m; --i) r *= half; */
  r = BLAS_dpow_di( b, m-1 );

  o -= eps;
  /* for (i = 0; i < l; ++i) o *= b; */
  o = (o * BLAS_dpow_di( b, l-1 )) * b;

  switch (cmach) {
    case blas_eps: return eps;
    case blas_sfmin: return r;
    default:
      BLAS_error( rname, -1, cmach, 0 );
      break;
  }
  return 0.0;
}

static int check_orthogonality(int, int, int, MORSE_Complex64_t*, double);
static int check_factorization(int, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, MORSE_Complex64_t*, double);
static int check_solution(int, int, int, MORSE_Complex64_t*, int, MORSE_Complex64_t*, MORSE_Complex64_t*, int, double);

int testing_zgels(int argc, char **argv)
{
    int hres = 0;
    int mode = 0;

    if ( argc < 1 ){
        goto usage;
    } else {
        mode = atoi(argv[0]);
    }

    /* Check for number of arguments*/
    if ( ((mode == 0) && (argc != 6)) ||
         ((mode != 0) && (argc != 7)) ){
      usage:
        USAGE("GELS", "MODE M N LDA NRHS LDB [RH]",
              "   - MODE : 0: flat, 1: tree (RH needed)\n"
              "   - M    : number of rows of the matrix A\n"
              "   - N    : number of columns of the matrix A\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - NRHS : number of RHS\n"
              "   - LDB  : leading dimension of the matrix B\n"
              "   - RH   : Size of each subdomains\n");
        return -1;
    }

    int M     = atoi(argv[1]);
    int N     = atoi(argv[2]);
    int LDA   = atoi(argv[3]);
    int NRHS  = atoi(argv[4]);
    int LDB   = atoi(argv[5]);
    int rh;

    int K = min(M, N);
    double eps;
    int info_ortho, info_solution, info_factorization;
    int i,j;
    int LDAxN = LDA*N;
    int LDBxNRHS = LDB*NRHS;

    MORSE_Complex64_t *A1 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *A2 = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B1 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *B2 = (MORSE_Complex64_t *)malloc(LDB*NRHS*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Q  = (MORSE_Complex64_t *)malloc(LDA*N*sizeof(MORSE_Complex64_t));
    MORSE_desc_t *T;

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)||(!Q)){
        printf("Out of Memory \n ");
        return -2;
    }

    if ( mode ) {
        rh = atoi(argv[6]);

        MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_TREE_HOUSEHOLDER);
        MORSE_Set(MORSE_HOUSEHOLDER_SIZE, rh);
    }

    MORSE_Alloc_Workspace_zgels(M, N, &T);
    memset(T->mat, 0, (T->llm*T->lln)*sizeof(MORSE_Complex64_t));
    eps = BLAS_dfpinfo( blas_eps );

    /*----------------------------------------------------------
    *  TESTING ZGELS
    */

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i] ;

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i] ;

    memset((void*)Q, 0, LDA*N*sizeof(MORSE_Complex64_t));
    for (i = 0; i < K; i++)
        Q[LDA*i+i] = 1.0;

    /* MORSE ZGELS */
    MORSE_zgels(MorseNoTrans, M, N, NRHS, A2, LDA, T, B2, LDB);

    /* MORSE ZGELS */
    if (M >= N)
       /* Building the economy-size Q */
       MORSE_zungqr(M, N, K, A2, LDA, T, Q, LDA);
    else
       /* Building the economy-size Q */
       MORSE_zunglq(M, N, K, A2, LDA, T, Q, LDA);

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZGELS ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    /* Check the orthogonality, factorization and the solution */
    info_ortho = check_orthogonality(M, N, LDA, Q, eps);
    info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
    info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution == 0)&(info_factorization == 0)&(info_ortho == 0)) {
        printf("***************************************************\n");
        printf(" ---- TESTING ZGELS ...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else {
        printf("************************************************\n");
        printf(" - TESTING ZGELS ... FAILED !\n");    hres++;
        printf("************************************************\n");
    }

    /*-------------------------------------------------------------
    *  TESTING ZGEQRF + ZGEQRS or ZGELQF + ZGELQS
    */

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
             B2[LDB*j+i] = B1[LDB*j+i];

    memset((void*)Q, 0, LDA*N*sizeof(MORSE_Complex64_t));
    for (i = 0; i < K; i++)
        Q[LDA*i+i] = 1.0;

    if (M >= N) {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGEQRF + ZGEQRS ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n", eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        /* Morse routines */
        MORSE_zgeqrf(M, N, A2, LDA, T);
        MORSE_zungqr(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zgeqrs(M, N, NRHS, A2, LDA, T, B2, LDB);

        /* Check the orthogonality, factorization and the solution */
        info_ortho = check_orthogonality(M, N, LDA, Q, eps);
        info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
        info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

        if ((info_solution == 0)&(info_factorization == 0)&(info_ortho == 0)) {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGEQRF + ZGEQRS ............ PASSED !\n");
            printf("***************************************************\n");
        }
        else{
            printf("***************************************************\n");
            printf(" - TESTING ZGEQRF + ZGEQRS ... FAILED !\n");
            printf("***************************************************\n");
        }
    }
    else  {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGELQF + ZGELQS ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n", eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        /* Morse routines */
        MORSE_zgelqf(M, N, A2, LDA, T);
        MORSE_zunglq(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zgelqs(M, N, NRHS, A2, LDA, T, B2, LDB);

       /* Check the orthogonality, factorization and the solution */
       info_ortho = check_orthogonality(M, N, LDA, Q, eps);
       info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
       info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

       if ( (info_solution == 0) & (info_factorization == 0) & (info_ortho == 0) ) {
          printf("***************************************************\n");
          printf(" ---- TESTING ZGELQF + ZGELQS ............ PASSED !\n");
          printf("***************************************************\n");
       }
       else {
          printf("***************************************************\n");
          printf(" - TESTING ZGELQF + ZGELQS ... FAILED !\n");
          printf("***************************************************\n");
        }
    }

    /*----------------------------------------------------------
    *  TESTING ZGEQRF + ZORMQR + ZTRSM
    */

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i];

    /* Initialize B1 and B2 */
    memset(B2, 0, LDB*NRHS*sizeof(MORSE_Complex64_t));
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
            B2[LDB*j+i] = B1[LDB*j+i];

    /* MORSE ZGEQRF+ ZUNMQR + ZTRSM */
    memset((void*)Q, 0, LDA*N*sizeof(MORSE_Complex64_t));
    for (i = 0; i < K; i++)
        Q[LDA*i+i] = 1.0;

    if (M >= N) {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGEQRF + ZUNMQR + ZTRSM  ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n",eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        MORSE_zgeqrf(M, N, A2, LDA, T);
        MORSE_zungqr(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zunmqr(MorseLeft, MorseConjTrans, M, NRHS, N, A2, LDA, T, B2, LDB);
        MORSE_ztrsm(MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, N, NRHS, 1.0, A2, LDA, B2, LDB);
    }
    else {
        printf("\n");
        printf("------ TESTS FOR CHAMELEON ZGELQF + ZUNMLQ + ZTRSM  ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", M, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n",eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");

        MORSE_zgelqf(M, N, A2, LDA, T);
        MORSE_ztrsm(MorseLeft, MorseLower, MorseNoTrans, MorseNonUnit, M, NRHS, 1.0, A2, LDA, B2, LDB);
        MORSE_zunglq(M, N, K, A2, LDA, T, Q, LDA);
        MORSE_zunmlq(MorseLeft, MorseConjTrans, N, NRHS, M, A2, LDA, T, B2, LDB);
    }

    /* Check the orthogonality, factorization and the solution */
    info_ortho = check_orthogonality(M, N, LDA, Q, eps);
    info_factorization = check_factorization(M, N, A1, A2, LDA, Q, eps);
    info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ( (info_solution == 0) & (info_factorization == 0) & (info_ortho == 0) ) {
        if (M >= N) {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGEQRF + ZUNMQR + ZTRSM .... PASSED !\n");
            printf("***************************************************\n");
        }
        else {
            printf("***************************************************\n");
            printf(" ---- TESTING ZGELQF + ZTRSM + ZUNMLQ .... PASSED !\n");
            printf("***************************************************\n");
        }
    }
    else {
        if (M >= N) {
            printf("***************************************************\n");
            printf(" - TESTING ZGEQRF + ZUNMQR + ZTRSM ... FAILED !\n");
            printf("***************************************************\n");
        }
        else {
            printf("***************************************************\n");
            printf(" - TESTING ZGELQF + ZTRSM + ZUNMLQ ... FAILED !\n");
            printf("***************************************************\n");
        }
    }

    free(A1); free(A2); free(B1); free(B2); free(Q);
    MORSE_Dealloc_Workspace( &T );

    return hres;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */

static int check_orthogonality(int M, int N, int LDQ, MORSE_Complex64_t *Q, double eps)
{
    double alpha, beta;
    double normQ;
    int info_ortho;
    int i;
    int minMN = min(M, N);

    double *work = (double *)malloc(minMN*sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    /* Build the idendity matrix USE DLASET?*/
    MORSE_Complex64_t *Id = (MORSE_Complex64_t *) malloc(minMN*minMN*sizeof(MORSE_Complex64_t));
    memset((void*)Id, 0, minMN*minMN*sizeof(MORSE_Complex64_t));
    for (i = 0; i < minMN; i++)
        Id[i*minMN+i] = (MORSE_Complex64_t)1.0;

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans, M, N, alpha, Q, LDQ, beta, Id, M);

    BLAS_zsy_norm( blas_colmajor, blas_inf_norm, blas_upper, minMN, Id, minMN, &normQ );

    printf("============\n");
    printf("Checking the orthogonality of Q \n");
    printf("||Id-Q'*Q||_oo / (N*eps) = %e \n", normQ/(minMN*eps));

    if ( isnan(normQ / (minMN * eps)) || isinf(normQ / (minMN * eps)) || (normQ / (minMN * eps) > 60.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho=1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho=0;
    }

    free(work); free(Id);

    return info_ortho;
}

/*------------------------------------------------------------
 *  Check the factorization QR
 */

static int check_factorization(int M, int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, MORSE_Complex64_t *Q, double eps )
{
    double Anorm, Rnorm;
    MORSE_Complex64_t alpha, beta;
    int info_factorization;
    int i,j;

    MORSE_Complex64_t *Ql       = (MORSE_Complex64_t *)malloc(M*N*sizeof(MORSE_Complex64_t));
    MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(M*N*sizeof(MORSE_Complex64_t));
    double *work              = (double *)malloc(max(M,N)*sizeof(double));

    alpha=1.0;
    beta=0.0;

    if (M >= N) {
        /* Extract the R */
        MORSE_Complex64_t *R = (MORSE_Complex64_t *)malloc(N*N*sizeof(MORSE_Complex64_t));
        memset((void*)R, 0, N*N*sizeof(MORSE_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'u', M, N, A2, LDA, R, N);

        /* Perform Ql=Q*R */
        memset((void*)Ql, 0, M*N*sizeof(MORSE_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, CBLAS_SADDR(alpha), Q, LDA, R, N, CBLAS_SADDR(beta), Ql, M);
        free(R);
    }
    else {
        /* Extract the L */
        MORSE_Complex64_t *L = (MORSE_Complex64_t *)malloc(M*M*sizeof(MORSE_Complex64_t));
        memset((void*)L, 0, M*M*sizeof(MORSE_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,'l', M, N, A2, LDA, L, M);

    /* Perform Ql=LQ */
        memset((void*)Ql, 0, M*N*sizeof(MORSE_Complex64_t));
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M, CBLAS_SADDR(alpha), L, M, Q, LDA, CBLAS_SADDR(beta), Ql, M);
        free(L);
    }

    /* Compute the Residual */
    for (i = 0; i < M; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*M+i] = A1[j*LDA+i]-Ql[j*M+i];

    BLAS_zge_norm( blas_colmajor, blas_inf_norm, M, N, Residual, M, &Rnorm );
    BLAS_zge_norm( blas_colmajor, blas_inf_norm, M, N, A2, LDA, &Anorm );

    if (M >= N) {
        printf("============\n");
        printf("Checking the QR Factorization \n");
        printf("-- ||A-QR||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }
    else {
        printf("============\n");
        printf("Checking the LQ Factorization \n");
        printf("-- ||A-LQ||_oo/(||A||_oo.N.eps) = %e \n",Rnorm/(Anorm*N*eps));
    }

    if (isnan(Rnorm / (Anorm * N *eps)) || isinf(Rnorm / (Anorm * N *eps)) || (Rnorm / (Anorm * N * eps) > 60.0) ) {
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else {
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(work); free(Ql); free(Residual);

    return info_factorization;
}

/*--------------------------------------------------------------
 * Check the solution
 */

static int check_solution(int M, int N, int NRHS, MORSE_Complex64_t *A1, int LDA, MORSE_Complex64_t *B1, MORSE_Complex64_t *B2, int LDB, double eps)
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    MORSE_Complex64_t alpha, beta;
    double result;
    double *work = (double *)malloc(max(M, N)* sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    BLAS_zge_norm( blas_colmajor, blas_inf_norm, M, N,    A1, LDA, &Anorm );
    BLAS_zge_norm( blas_colmajor, blas_inf_norm, M, NRHS, B2, LDB, &Xnorm );
    BLAS_zge_norm( blas_colmajor, blas_inf_norm, N, NRHS, B1, LDB, &Bnorm );

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, CBLAS_SADDR(alpha), A1, LDA, B2, LDB, CBLAS_SADDR(beta), B1, LDB);

    if (M >= N) {
       MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(M*NRHS*sizeof(MORSE_Complex64_t));
       memset((void*)Residual, 0, M*NRHS*sizeof(MORSE_Complex64_t));
       cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A1, LDA, B1, LDB, CBLAS_SADDR(beta), Residual, M);
       BLAS_zge_norm( blas_colmajor, blas_inf_norm, M, NRHS, Residual, M, &Rnorm );
       free(Residual);
    }
    else {
       MORSE_Complex64_t *Residual = (MORSE_Complex64_t *)malloc(N*NRHS*sizeof(MORSE_Complex64_t));
       memset((void*)Residual, 0, N*NRHS*sizeof(MORSE_Complex64_t));
       cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A1, LDA, B1, LDB, CBLAS_SADDR(beta), Residual, N);
       BLAS_zge_norm( blas_colmajor, blas_inf_norm, N, NRHS, Residual, N, &Rnorm );
       free(Residual);
    }

    if (getenv("MORSE_TESTING_VERBOSE"))
       printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*N*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }
    free(work);
    return info_solution;
}
