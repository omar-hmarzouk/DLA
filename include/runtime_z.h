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
 * @file runtime_z.h
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifndef _RUNTIME_ZBLAS_H_
#define _RUNTIME_ZBLAS_H_

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif


/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by MORSE) - alphabetical order
 **/
void MORSE_TASK_zdasum(MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, int M, int N,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zgeadd(MORSE_option_t *options,
                      int m, int n, MORSE_Complex64_t alpha,
                      MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zbrdalg(MORSE_option_t *options,
                        MORSE_enum uplo,
                        int N, int NB,
                        MORSE_desc_t *A,
                        MORSE_desc_t *C, int Cm, int Cn,
                        MORSE_desc_t *S, int Sm, int Sn,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void MORSE_TASK_zgelqt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zgemm(MORSE_option_t *options,
                      MORSE_enum transA, MORSE_enum transB,
                      int m, int n, int k, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zgemm2( MORSE_option_t *options,
                        MORSE_enum transA, MORSE_enum transB,
                        int m, int n, int k, int nb,
                        MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_desc_t *B, int Bm, int Bn, int ldb,
                        MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zgemm_f2(MORSE_option_t *options,
                         MORSE_enum transA, MORSE_enum transB,
                         int m, int n, int k, int nb,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                         MORSE_desc_t *B, int Bm, int Bn, int ldb,
                         MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc,
                         MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                         MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2);
void MORSE_TASK_zgemm_p2(MORSE_option_t *options,
                         MORSE_enum transA, MORSE_enum transB,
                         int m, int n, int k, int nb,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_Complex64_t **B, int ldb,
                         MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zgemm_p2f1(MORSE_option_t *options,
                           MORSE_enum transA, MORSE_enum transB,
                           int m, int n, int k, int nb,
                           MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                           const MORSE_Complex64_t **B, int ldb,
                           MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc,
                           MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1);
void MORSE_TASK_zgemm_p3(MORSE_option_t *options,
                         MORSE_enum transA, MORSE_enum transB,
                         int m, int n, int k, int nb,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                         MORSE_desc_t *B, int Bm, int Bn, int ldb,
                         MORSE_Complex64_t beta, MORSE_Complex64_t **C, int ldc);
void MORSE_TASK_zgeqrt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zgessm(MORSE_option_t *options,
                       int m, int n, int k, int ib, int nb,
                       int *IPIV,
                       MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       MORSE_desc_t *D, int Dm, int Dn, int ldd,
                       MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_zgessq( MORSE_option_t *options,
                        int m, int n,
                        MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn);
void MORSE_TASK_zgetrf(MORSE_option_t *options,
                       int m, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       int *IPIV,

                       MORSE_bool check_info, int iinfo);
void MORSE_TASK_zgetrf_incpiv(MORSE_option_t *options,
                              int m, int n, int ib, int nb,
                              MORSE_desc_t *A, int Am, int An, int lda,
                              MORSE_desc_t *L, int Lm, int Ln, int ldl,
                              int *IPIV,
                              MORSE_bool check_info, int iinfo);
void MORSE_TASK_zgetrf_reclap(MORSE_option_t *options,
                              int m, int n, int nb,
                              MORSE_desc_t *A, int Am, int An, int lda,
                              int *IPIV,

                              MORSE_bool check_info, int iinfo,
                              int nbthread);
void MORSE_TASK_zgetrf_rectil(MORSE_option_t *options,
                              MORSE_desc_t A, MORSE_desc_t *Amn, int Amnm, int Amnn, int size,
                              int *IPIV,

                              MORSE_bool check_info, int iinfo,
                              int nbthread);
void MORSE_TASK_zgetrip(MORSE_option_t *options,
                        int m, int n, MORSE_desc_t *A, int Am, int An, int szeA);
void MORSE_TASK_zgetrip_f1(MORSE_option_t *options,
                           int m, int n, MORSE_desc_t *A, int Am, int An, int szeA,
                           MORSE_desc_t *fake, int fakem, int faken, int szeF, int paramF);
void MORSE_TASK_zgetrip_f2(MORSE_option_t *options,
                           int m, int n, MORSE_desc_t *A, int Am, int An, int szeA,
                           MORSE_desc_t *fake1, int fake1m, int fake1n, int szeF1, int paramF1,
                           MORSE_desc_t *fake2, int fake2m, int fake2n, int szeF2, int paramF2);
void MORSE_TASK_zhemm(MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zhegst(MORSE_option_t *options,
                       int itype, MORSE_enum uplo, int N,
                       MORSE_desc_t *A, int Am, int An, int LDA,
                       MORSE_desc_t *B, int Bm, int Bn, int LDB,

                       int iinfo);
void MORSE_TASK_zherk(MORSE_option_t *options,
                      MORSE_enum uplo, MORSE_enum trans,
                      int n, int k, int nb,
                      double alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      double beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zher2k(MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int LDB,
                       double beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zherfb(MORSE_option_t *options,
                       MORSE_enum uplo,
                       int n, int k, int ib, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zlacpy(MORSE_option_t *options,
                       MORSE_enum uplo, int m, int n, int mb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zlange(MORSE_option_t *options,
                       MORSE_enum norm, int M, int N, int NB,
                       MORSE_desc_t *A, int Am, int An, int LDA,
                       MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlange_max(MORSE_option_t *options,
                           MORSE_desc_t *A, int Am, int An,
                           MORSE_desc_t *B, int Bm, int Bn);
#ifdef COMPLEX
void MORSE_TASK_zhessq( MORSE_option_t *options,
                        MORSE_enum uplo, int n,
                        MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void MORSE_TASK_zlanhe(MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                       MORSE_desc_t *A, int Am, int An, int LDA,
                       MORSE_desc_t *B, int Bm, int Bn);
#endif
void MORSE_TASK_zlansy(MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                       MORSE_desc_t *A, int Am, int An, int LDA,
                       MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlantr(MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag,
                       int M, int N, int NB,
                       MORSE_desc_t *A, int Am, int An, int LDA,
                       MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_zlaset(MORSE_option_t *options,
                       MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                       MORSE_Complex64_t beta, MORSE_desc_t *tileA, int tileAm, int tileAn, int ldtilea);
void MORSE_TASK_zlaset2(MORSE_option_t *options,
                        MORSE_enum uplo, int n1, int n2, MORSE_Complex64_t alpha,
                        MORSE_desc_t *tileA, int tileAm, int tileAn, int ldtilea);
void MORSE_TASK_zlaswp(MORSE_option_t *options,
                       int n, MORSE_desc_t *A, int Am, int An, int lda,
                       int i1,  int i2, int *ipiv, int inc);
void MORSE_TASK_zlaswp_f2(MORSE_option_t *options,
                          int n, MORSE_desc_t *A, int Am, int An, int lda,
                          int i1,  int i2, int *ipiv, int inc,
                          MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                          MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2);
void MORSE_TASK_zlaswp_ontile(MORSE_option_t *options,
                              MORSE_desc_t descA, MORSE_desc_t *A, int Am, int An,
                              int i1,  int i2, int *ipiv, int inc, MORSE_Complex64_t *fakepanel);
void MORSE_TASK_zlaswp_ontile_f2(MORSE_option_t *options,
                                 MORSE_desc_t descA, MORSE_desc_t *A, int Am, int An,
                                 int i1,  int i2, int *ipiv, int inc,
                                 MORSE_desc_t *fake1, int fake1m, int fake1n, int szefake1, int flag1,
                                 MORSE_desc_t *fake2, int fake2m, int fake2n, int szefake2, int flag2);
void MORSE_TASK_zlaswpc_ontile(MORSE_option_t *options,
                               MORSE_desc_t descA, MORSE_desc_t *A, int Am, int An,
                               int i1,  int i2, int *ipiv, int inc, MORSE_Complex64_t *fakepanel);
void MORSE_TASK_zlatro(MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans, int m, int n, int mb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_zlauum(MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_zplghe(MORSE_option_t *options,
                       double bump, int m, int n, MORSE_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_TASK_zplgsy(MORSE_option_t *options,
                       MORSE_Complex64_t bump, int m, int n, MORSE_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_TASK_zplrnt(MORSE_option_t *options,
                       int m, int n, MORSE_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed );
void MORSE_TASK_zpotrf(MORSE_option_t *options,
                       MORSE_enum uplo, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,

                       int iinfo);
void MORSE_TASK_zshift( MORSE_option_t *options,
                        int s, int m, int n, int L,
                        MORSE_Complex64_t *A);
void MORSE_TASK_zshiftw(MORSE_option_t *options,
                        int s, int cl, int m, int n, int L,
                        MORSE_desc_t *A, int Am, int An, MORSE_Complex64_t *W);
void MORSE_TASK_zssssm(MORSE_option_t *options,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *L1, int L1m, int L1n, int ldl1,
                       MORSE_desc_t *L2, int L2m, int L2n, int ldl2,
                       const int *IPIV);
void MORSE_TASK_zsymm(MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zsyrk(MORSE_option_t *options,
                      MORSE_enum uplo, MORSE_enum trans,
                      int n, int k, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zsyr2k(MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum trans,
                       int n, int k, int nb,
                       MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int LDB,
                       MORSE_Complex64_t beta, MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zsyssq( MORSE_option_t *options,
                        MORSE_enum uplo, int n,
                        MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void MORSE_TASK_zswpab(MORSE_option_t *options,
                       int i, int n1, int n2,
                       MORSE_desc_t *A, int Am, int An, int szeA);
void MORSE_TASK_zswptr_ontile(MORSE_option_t *options,
                              MORSE_desc_t descA, MORSE_desc_t *Aij, int Aijm, int Aijn,
                              int i1,  int i2, int *ipiv, int inc,
                              MORSE_desc_t *Akk, int Akkm, int Akkn, int ldak);
void MORSE_TASK_ztrdalg(MORSE_option_t *options,
                        MORSE_enum uplo,
                        int N, int NB,
                        MORSE_desc_t *A,
                        MORSE_desc_t *C, int Cm, int Cn,
                        MORSE_desc_t *S, int Sm, int Sn,
                        int i, int j, int m, int grsiz, int BAND,
                        int *PCOL, int *ACOL, int *MCOL);
void MORSE_TASK_ztrasm(MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn);
void MORSE_TASK_ztrmm(MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_ztrmm_p2(MORSE_option_t *options,
                         MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                         int m, int n, int nb,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                         MORSE_Complex64_t **B, int ldb);
void MORSE_TASK_ztrsm(MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, int n, int nb,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A, int Am, int An, int lda,
                      MORSE_desc_t *B, int Bm, int Bn, int ldb);
void MORSE_TASK_ztrssq( MORSE_option_t *options,
                        MORSE_enum uplo, MORSE_enum diag,
                        int m, int n,
                        MORSE_desc_t *A, int Am, int An, int lda,
                        MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn );
void MORSE_TASK_ztrtri(MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum diag, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,

                       int iinfo);
void MORSE_TASK_ztslqt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmlq(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmlq_hetra1(MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmlq_corner(MORSE_option_t *options,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              MORSE_desc_t *A3, int A3m, int A3n, int lda3,
                              MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmqr(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmqr_hetra1(MORSE_option_t *options,
                              MORSE_enum side, MORSE_enum trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsmqr_corner(MORSE_option_t *options,
                              int m1, int n1, int m2, int n2, int m3, int n3, int k, int ib, int nb,
                              MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                              MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                              MORSE_desc_t *A3, int A3m, int A3n, int lda3,
                              MORSE_desc_t *V, int Vm, int Vn, int ldv,
                              MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztsqrt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_ztstrf(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *U, int Um, int Un, int ldu,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *L, int Lm, int Ln, int ldl,
                       int *IPIV,

                       MORSE_bool check_info, int iinfo);
void MORSE_TASK_zttmqr(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zttqrt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zttmlq(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m1, int n1, int m2, int n2, int k, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zttlqt(MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt);
void MORSE_TASK_zpamm(MORSE_option_t *options,
                       int op, MORSE_enum side, MORSE_enum storev,
                       int m, int n, int k, int l,
                       MORSE_desc_t *A1, int A1m, int A1n, int lda1,
                       MORSE_desc_t *A2, int A2m, int A2n, int lda2,
                       MORSE_desc_t *V, int Vm, int Vn, int ldv,
                       MORSE_desc_t *W, int Wm, int Wn, int ldw);
void MORSE_TASK_zplssq( MORSE_option_t *options,
                        MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                        MORSE_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn );
void MORSE_TASK_zplssq2( MORSE_option_t *options,
                         MORSE_desc_t *RESULT, int RESULTm, int RESULTn );
void MORSE_TASK_zunmlq(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int ib,  int nb, int k,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       MORSE_desc_t *C, int Cm, int Cn, int ldc);
void MORSE_TASK_zunmqr(MORSE_option_t *options,
                       MORSE_enum side, MORSE_enum trans,
                       int m, int n, int k, int ib, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *T, int Tm, int Tn, int ldt,
                       MORSE_desc_t *C, int Cm, int Cn, int ldc);



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
