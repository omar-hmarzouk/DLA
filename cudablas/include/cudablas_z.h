/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file cudablas_z.h
 *
 *  MORSE cudablas headers
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 **/
#ifndef _MORSE_CUDA_ZBLAS_H_
#define _MORSE_CUDA_ZBLAS_H_

#include "morse_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif

/** ****************************************************************************
 *  Declarations of cuda kernels - alphabetical order
 **/
#if defined(CHAMELEON_USE_MAGMA)
int CUDA_zgelqt(
        magma_int_t m, magma_int_t n, magma_int_t nb,
        magmaDoubleComplex *da, magma_int_t ldda,
        magmaDoubleComplex *v,  magma_int_t ldv,
        magmaDoubleComplex *dt, magma_int_t lddt,
        magmaDoubleComplex *t,  magma_int_t ldt,
        magmaDoubleComplex *dd,
        magmaDoubleComplex *d,  magma_int_t ldd,
        magmaDoubleComplex *tau,
        magmaDoubleComplex *hwork,
        magmaDoubleComplex *dwork,
        CUstream stream);
int CUDA_zgemerge(
        magma_side_t side, magma_diag_t diag,
        magma_int_t M, magma_int_t N,
        magmaDoubleComplex *A, magma_int_t LDA,
        magmaDoubleComplex *B, magma_int_t LDB,
        CUstream stream);
int CUDA_zgemm_V2(
        MORSE_enum transa, MORSE_enum transb,
        int m, int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zgemm(
        MORSE_enum transa, MORSE_enum transb,
        int m, int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zgeqrt(
        magma_int_t m, magma_int_t n, magma_int_t nb,
        magmaDoubleComplex *da, magma_int_t ldda,
        magmaDoubleComplex *v,  magma_int_t ldv,
        magmaDoubleComplex *dt, magma_int_t lddt,
        magmaDoubleComplex *t,  magma_int_t ldt,
        magmaDoubleComplex *dd,
        magmaDoubleComplex *d,  magma_int_t ldd,
        magmaDoubleComplex *tau,
        magmaDoubleComplex *hwork,
        magmaDoubleComplex *dwork,
        CUstream stream);
int CUDA_zgessm(
        char storev, magma_int_t m, magma_int_t n, magma_int_t k, magma_int_t ib,
        magma_int_t *ipiv,
        cuDoubleComplex *dL1, magma_int_t lddl1,
        cuDoubleComplex *dL,  magma_int_t lddl,
        cuDoubleComplex *dA,  magma_int_t ldda,
        magma_int_t *info);
int CUDA_zgetrf_incpiv(
        char storev, magma_int_t m, magma_int_t n, magma_int_t ib,
        cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda,
        cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
        magma_int_t *ipiv,
        cuDoubleComplex *dwork, magma_int_t lddwork,
        magma_int_t *info);
int CUDA_zgetrf_nopiv(
        magma_int_t m, magma_int_t n,
        cuDoubleComplex *dA, magma_int_t ldda,
        magma_int_t *info);
int CUDA_zhemm_V2(
        MORSE_enum side, MORSE_enum uplo,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zhemm(
        MORSE_enum side, MORSE_enum uplo,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zher2k_V2(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        double *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zher2k(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        double *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zherk_V2(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        double *alpha,
        const cuDoubleComplex *A, int lda,
        double *beta,
        cuDoubleComplex *B, int ldb,
        CUstream stream);
int CUDA_zherk(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        double *alpha,
        const cuDoubleComplex *A, int lda,
        double *beta,
        cuDoubleComplex *B, int ldb,
        CUstream stream);
int CUDA_zlauum(
        char uplo,  magma_int_t n,
        cuDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
int CUDA_zparfb(
        magma_side_t side, magma_trans_t trans,
        magma_direct_t direct, magma_storev_t storev,
        magma_int_t M1, magma_int_t N1,
        magma_int_t M2, magma_int_t N2,
        magma_int_t K, magma_int_t L,
        magmaDoubleComplex *A1, magma_int_t LDA1,
        magmaDoubleComplex *A2, magma_int_t LDA2,
        const magmaDoubleComplex *V, magma_int_t LDV,
        const magmaDoubleComplex *T, magma_int_t LDT,
        magmaDoubleComplex *WORK, magma_int_t LDWORK,
        magmaDoubleComplex *WORKC, magma_int_t LDWORKC,
        CUstream stream);
int CUDA_zpotrf(
        magma_uplo_t uplo,  magma_int_t n,
        magmaDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
int CUDA_zssssm(
        magma_storev_t storev, magma_int_t m1, magma_int_t n1,
        magma_int_t m2, magma_int_t n2, magma_int_t k, magma_int_t ib,
        magmaDoubleComplex *dA1, magma_int_t ldda1,
        magmaDoubleComplex *dA2, magma_int_t ldda2,
        magmaDoubleComplex *dL1, magma_int_t lddl1,
        magmaDoubleComplex *dL2, magma_int_t lddl2,
        magma_int_t *IPIV, magma_int_t *info);
int CUDA_zsymm_V2(
        MORSE_enum side, MORSE_enum uplo,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zsymm(
        MORSE_enum side, MORSE_enum uplo,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zsyr2k_V2(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zsyr2k(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zsyrk_V2(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_zsyrk(
        MORSE_enum uplo, MORSE_enum trans,
        int n, int k,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *beta,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_ztrmm_V2(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        const cuDoubleComplex *B, int ldb,
        cuDoubleComplex *C, int ldc,
        CUstream stream);
int CUDA_ztrmm(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        CUstream stream);
int CUDA_ztrsm_V2(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        const cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        CUstream stream);
int CUDA_ztrsm(
        MORSE_enum side, MORSE_enum uplo,
        MORSE_enum transa, MORSE_enum diag,
        int m, int n,
        cuDoubleComplex *alpha,
        const cuDoubleComplex *A, int lda,
        cuDoubleComplex *B, int ldb,
        CUstream stream);
int CUDA_ztrtri(
        magma_uplo_t uplo,  magma_diag_t diag, magma_int_t n,
        magmaDoubleComplex *dA, magma_int_t ldda, magma_int_t *info);
int CUDA_ztslqt(
        magma_int_t m, magma_int_t n, magma_int_t nb,
        magmaDoubleComplex *da1, magma_int_t ldda1,
        magmaDoubleComplex *da2, magma_int_t ldda2,
        magmaDoubleComplex *a2,  magma_int_t lda2,
        magmaDoubleComplex *dt,  magma_int_t lddt,
        magmaDoubleComplex *t,  magma_int_t ldt,
        magmaDoubleComplex *dd,
        magmaDoubleComplex *d,  magma_int_t ldd,
        magmaDoubleComplex *tau,
        magmaDoubleComplex *hwork,
        magmaDoubleComplex *dwork,
        CUstream stream);
int CUDA_ztsmlq(
        magma_side_t side, magma_trans_t trans,
        magma_int_t M1, magma_int_t N1,
        magma_int_t M2, magma_int_t N2,
        magma_int_t K, magma_int_t IB,
        magmaDoubleComplex *A1, magma_int_t LDA1,
        magmaDoubleComplex *A2, magma_int_t LDA2,
        const magmaDoubleComplex *V, magma_int_t LDV,
        const magmaDoubleComplex *T, magma_int_t LDT,
              magmaDoubleComplex *WORK,  magma_int_t LDWORK,
              magmaDoubleComplex *WORKC, magma_int_t LDWORKC,
        CUstream stream);
int CUDA_ztsmqr(
        magma_side_t side, magma_trans_t trans,
        magma_int_t M1, magma_int_t N1,
        magma_int_t M2, magma_int_t N2,
        magma_int_t K, magma_int_t IB,
        magmaDoubleComplex *A1, magma_int_t LDA1,
        magmaDoubleComplex *A2, magma_int_t LDA2,
        const magmaDoubleComplex *V, magma_int_t LDV,
        const magmaDoubleComplex *T, magma_int_t LDT,
            magmaDoubleComplex *WORK,  magma_int_t LDWORK,
            magmaDoubleComplex *WORKC, magma_int_t LDWORKC,
        CUstream stream);
int CUDA_ztsqrt(
        magma_int_t m, magma_int_t n, magma_int_t nb,
        magmaDoubleComplex *da1, magma_int_t ldda1,
        magmaDoubleComplex *da2, magma_int_t ldda2,
        magmaDoubleComplex *a2,  magma_int_t lda2,
        magmaDoubleComplex *dt,  magma_int_t lddt,
        magmaDoubleComplex *t,  magma_int_t ldt,
        magmaDoubleComplex *dd,
        magmaDoubleComplex *d,  magma_int_t ldd,
        magmaDoubleComplex *tau,
        magmaDoubleComplex *hwork,
        magmaDoubleComplex *dwork,
        CUstream stream);
int CUDA_ztstrf(
        char storev, magma_int_t m, magma_int_t n, magma_int_t ib, magma_int_t nb,
        cuDoubleComplex *hU, magma_int_t ldhu, cuDoubleComplex *dU, magma_int_t lddu,
        cuDoubleComplex *hA, magma_int_t ldha, cuDoubleComplex *dA, magma_int_t ldda,
        cuDoubleComplex *hL, magma_int_t ldhl, cuDoubleComplex *dL, magma_int_t lddl,
        magma_int_t *ipiv,
        cuDoubleComplex *hwork, magma_int_t ldhwork,
        cuDoubleComplex *dwork, magma_int_t lddwork,
        magma_int_t *info);
int CUDA_zunmlqt(
        magma_side_t side, magma_trans_t trans,
        magma_int_t M, magma_int_t N, magma_int_t K, magma_int_t IB,
        const magmaDoubleComplex *A,    magma_int_t LDA,
        const magmaDoubleComplex *T,    magma_int_t LDT,
        magmaDoubleComplex *C,    magma_int_t LDC,
        magmaDoubleComplex *WORK, magma_int_t LDWORK );
int CUDA_zunmqrt(
        magma_side_t side, magma_trans_t trans,
        magma_int_t M, magma_int_t N, magma_int_t K, magma_int_t IB,
        const magmaDoubleComplex *A,    magma_int_t LDA,
        const magmaDoubleComplex *T,    magma_int_t LDT,
        magmaDoubleComplex *C,    magma_int_t LDC,
        magmaDoubleComplex *WORK, magma_int_t LDWORK );
#endif

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
