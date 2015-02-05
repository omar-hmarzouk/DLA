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
 * @file simucblas.c
 *
 * MORSE fake cblas functions for simulation mode
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-01
 *
 **/

#include "simucore/simucblas/cblas.h"
#include <assert.h>

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
float  cblas_sdsdot(const int N, const float alpha, const float *X,
                    const int incX, const float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
double cblas_dsdot(const int N, const float *X, const int incX, const float *Y,
                   const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
float  cblas_sdot(const int N, const float  *X, const int incX,
                  const float  *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

/*
 * Functions having prefixes Z and C only
 */
void   cblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void   cblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void   cblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void   cblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/*
 * Functions having prefixes S D SC DZ
 */
float  cblas_snrm2(const int N, const float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
float  cblas_sasum(const int N, const float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

double cblas_dnrm2(const int N, const double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
double cblas_dasum(const int N, const double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

float  cblas_scnrm2(const int N, const void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
float  cblas_scasum(const int N, const void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

double cblas_dznrm2(const int N, const void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
double cblas_dzasum(const int N, const void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX cblas_isamax(const int N, const float  *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
CBLAS_INDEX cblas_idamax(const int N, const double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
CBLAS_INDEX cblas_icamax(const int N, const void   *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
CBLAS_INDEX cblas_izamax(const int N, const void   *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void cblas_sswap(const int N, float *X, const int incX, 
                 float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_scopy(const int N, const float *X, const int incX, 
                 float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_dswap(const int N, double *X, const int incX, 
                 double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dcopy(const int N, const double *X, const int incX, 
                 double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_cswap(const int N, void *X, const int incX, 
                 void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ccopy(const int N, const void *X, const int incX, 
                 void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_caxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_zswap(const int N, void *X, const int incX, 
                 void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zcopy(const int N, const void *X, const int incX, 
                 void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/* 
 * Routines with S and D prefix only
 */
void cblas_srotg(float *a, float *b, float *c, float *s){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_srot(const int N, float *X, const int incX,
                float *Y, const int incY, const float c, const float s){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_srotm(const int N, float *X, const int incX,
                 float *Y, const int incY, const float *P){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_drotg(double *a, double *b, double *c, double *s){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_drot(const int N, double *X, const int incX,
                double *Y, const int incY, const double c, const double  s){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_drotm(const int N, double *X, const int incX,
                 double *Y, const int incY, const double *P){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/* 
 * Routines with S D C Z CS and ZD prefixes
 */
void cblas_sscal(const int N, const float alpha, float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dscal(const int N, const double alpha, double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cscal(const int N, const void *alpha, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zscal(const int N, const void *alpha, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_csscal(const int N, const float alpha, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zdscal(const int N, const double alpha, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_sgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const float alpha,
                 const float *A, const int lda, const float *X,
                 const int incX, const float beta, float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_strmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *A, const int lda, 
                 float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_stbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda, 
                 float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_stpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_strsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *A, const int lda, float *X,
                 const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_stbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const float *A, const int lda,
                 float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_stpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const float *Ap, float *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const double alpha,
                 const double *A, const int lda, const double *X,
                 const int incX, const double beta, double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda, 
                 double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda, 
                 double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda, double *X,
                 const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const double *A, const int lda,
                 double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const double *Ap, double *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_cgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, 
                 void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda, 
                 void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_zgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const int KL, const int KU, const void *alpha,
                 const void *A, const int lda, const void *X,
                 const int incX, const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, 
                 void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda, 
                 void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *A, const int lda, void *X,
                 const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const int K, const void *A, const int lda,
                 void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const int N, const void *Ap, void *X, const int incX){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/* 
 * Routines with S and D prefixes only
 */
void cblas_ssymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ssbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const float alpha, const float *A,
                 const int lda, const float *X, const int incX,
                 const float beta, float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_sspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *Ap,
                 const float *X, const int incX,
                 const float beta, float *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_sger(const enum CBLAS_ORDER order, const int M, const int N,
                const float alpha, const float *X, const int incX,
                const float *Y, const int incY, float *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ssyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_sspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const float *X,
                const int incX, float *Ap){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ssyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *X,
                 const int incX, const float *Y, const int incY, float *A,
                 const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_sspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const float alpha, const float *X,
                 const int incX, const float *Y, const int incY, float *A){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_dsymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dsbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *Ap,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double *X, const int incX,
                const double *Y, const int incY, double *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dsyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dsyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *X,
                 const int incX, const double *Y, const int incY, double *A,
                 const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const double alpha, const double *X,
                 const int incX, const double *Y, const int incY, double *A){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/* 
 * Routines with C and Z prefixes only
 */
void cblas_chemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_chbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_chpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cgeru(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cgerc(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const void *X, const int incX,
                void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_chpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const float alpha, const void *X,
                const int incX, void *A){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_chpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *Ap){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_zhemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zhbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const int K, const void *alpha, const void *A,
                 const int lda, const void *X, const int incX,
                 const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zhpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const int N, const void *alpha, const void *Ap,
                 const void *X, const int incX,
                 const void *beta, void *Y, const int incY){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zgeru(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zgerc(const enum CBLAS_ORDER order, const int M, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const void *X, const int incX,
                void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zhpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const void *X,
                const int incX, void *A){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *A, const int lda){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zhpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const int N,
                 const void *alpha, const void *X, const int incX,
                 const void *Y, const int incY, void *Ap){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *B, const int ldb, const float beta,
                 float *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const float *A, const int lda,
                 const float beta, float *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const float alpha, const float *A, const int lda,
                  const float *B, const int ldb, const float beta,
                  float *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 float *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *B, const int ldb, const double beta,
                 double *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const double alpha, const double *A, const int lda,
                  const double *B, const int ldb, const double beta,
                  double *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const void *alpha, const void *A, const int lda,
                 const void *beta, void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const void *beta,
                  void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 void *B, const int ldb){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}


/* 
 * Routines with prefixes C and Z only
 */
void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const float alpha, const void *A, const int lda,
                 const float beta, void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const float beta,
                  void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *B, const int ldb, const void *beta,
                 void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const void *A, const int lda,
                 const double beta, void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
                  const void *alpha, const void *A, const int lda,
                  const void *B, const int ldb, const double beta,
                  void *C, const int ldc){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}

void cblas_xerbla(int p, const char *rout, const char *form, ...){assert(0 && "Error: you should not enter to this function, simucblas is a fake library for use with a simulator like Simgrid. ");}
