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
 * @file cuda_zgeqrt.c
 *
 *  MORSE cudablas kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @author Florent Pruvost
 * @date 2015-09-16
 * @precisions normal z -> c d s
 *
 **/
#include "cudablas/include/cudablas.h"

#if defined(CHAMELEON_USE_MAGMA)
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
        CUstream stream)
{
#define da_ref(a_1,a_2) ( da+(a_2)*(ldda) + (a_1))
#define v_ref(a_1,a_2)  ( v+(a_2)*(ldv) + (a_1))
#define dt_ref(a_1,a_2) ( dt+(a_2)*(lddt) + (a_1))
#define t_ref(a_1,a_2)  ( t+(a_2)*(ldt) + (a_1))

    int i, k, ib, lddwork, old_i, old_ib, rows, cols;
    double _Complex one=1.;

//    int lwkopt = n * nb;
//    hwork[0] = *((magmaDoubleComplex*) &lwkopt);

    if (m < 0) {
        return -1;
    } else if (n < 0) {
        return -2;
    } else if (ldda < max(1,m)) {
        return -4;
    }

    k = min(m,n);
    if (k == 0) {
        hwork[0] = *((magmaDoubleComplex*) &one);
        return MAGMA_SUCCESS;
    }

    lddwork= k;

    /* lower parts of little T must be zero: memset to 0 for simplicity */
    memset(t_ref(0,0), 0, nb*nb*sizeof(magmaDoubleComplex));
    cudaMemsetAsync(dt_ref(0,0), 0, nb*n*sizeof(magmaDoubleComplex), stream);

    /* copy first panel of A on the host */
//  cublasGetMatrix(m, min(nb,n), sizeof(magmaDoubleComplex),
//                    da_ref(0, 0), ldda,
//                    v, ldv);
    cudaMemcpy( v, da_ref(0,0),
                m*min(nb,n)*sizeof(magmaDoubleComplex),
                cudaMemcpyDeviceToHost);

    /* Use blocked code initially */
    for (i = 0; i < k; i += nb) {

        ib = min(k-i, nb);
        if (i+nb>=n) ib = min(n-i, nb);
        rows = m-i;

        if (i>0){

            /* copy panel of A from device to host */
//          cublasGetMatrix(m, ib, sizeof(magmaDoubleComplex),
//                            da_ref(0, i), ldda,
//                            v, ldv);
            /* copy panel of A from device to host */
            cudaMemcpy( v, da_ref(0,i),
                        m*ib*sizeof(magmaDoubleComplex),
                        cudaMemcpyDeviceToHost);

            /* Apply H' to A(i:m,i+2*ib:n) from the left */
            cols = n-old_i-2*old_ib;
            if (cols > 0){
                magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                  m-old_i, cols, old_ib,
                                  da_ref(old_i, old_i), ldda, dt_ref(0,old_i), lddt,
                                  da_ref(old_i, old_i+2*old_ib), ldda,
                                  dwork, cols);
            }

            /* copy the upper diag tile into d_A */
            CUDA_zgemerge(MagmaLeft, MagmaUnit, old_ib, old_ib,
                          dd, ldd, da_ref(old_i, old_i), ldda, stream);

        }

        /* Form the triangular factor of the block reflector on the host
         H = H(i) H(i+1) . . . H(i+ib-1) */
        CORE_zgeqrt(rows, ib, ib,
                    (double _Complex*) v_ref(i,0), ldv,
                    (double _Complex*) t_ref(0,0), ldt,
                    (double _Complex*) tau+i,
                    (double _Complex*) hwork);

        if ( i + ib < n ){
            /* put 0s in the upper triangular part of a panel (and 1s on the
             diagonal); copy the upper triangular in d */
            CORE_zgesplit(MorseLeft, MorseUnit, min(rows,ib), ib,
                          (double _Complex*) v_ref(i,0), ldv,
                          (double _Complex*) d, ldd);

            /* copy from host to device a tile diag */
            cublasSetMatrix( min(rows,ib), ib, sizeof(magmaDoubleComplex),
                             d, ldd, dd, ldd );
        }

        /* Send the triangular factor T to the GPU */
        cublasSetMatrix( ib, ib, sizeof(magmaDoubleComplex),
                         t_ref(0,0), ldt, dt_ref(0,i), lddt );

        /* A panel (with zeros in upper tri of its diag) is ready to be used
         in input of zlarfb_gpu: we send the panel to the gpu */
        cublasSetMatrix( rows, ib, sizeof(magmaDoubleComplex),
                         v_ref(i,0), ldv, da_ref(i,i), ldda );

        if (i + ib < n) {

            if (i+2*ib < n){
                cols = ib;
            }
            else{
                cols = n-i-ib;
            }
            /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
            magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                              rows, cols, ib, da_ref(i,i), ldda, dt_ref(0,i),
                              lddt, da_ref(i,i+ib), ldda, dwork, cols);
            old_i = i;
            old_ib = ib;
            if (i+nb>=k){
                /* Apply H' to A(i:m,i+2*ib:n) from the left */
                cols = n-old_i-2*old_ib;
                if (cols > 0){
                    magma_zlarfb_gpu( MagmaLeft, MagmaConjTrans, MagmaForward, MagmaColumnwise,
                                      rows, cols, old_ib,
                                      da_ref(old_i, old_i), ldda, dt_ref(0,old_i), lddt,
                                      da_ref(old_i, old_i+2*old_ib), ldda,
                                      dwork, cols);
                }
                /* copy the upper diag tile into d_A */
                CUDA_zgemerge(MagmaLeft, MagmaUnit, old_ib, old_ib,
                              dd, ldd, da_ref(old_i, old_i), ldda, stream);
            }
        }

    }

#undef da_ref
#undef v_ref
#undef dt_ref
#undef t_ref

    return MORSE_SUCCESS;
}

#if defined(CHAMELEON_USE_CUBLAS_V2)
int CUDA_zgemerge(
        magma_side_t side, magma_diag_t diag,
        magma_int_t M, magma_int_t N,
        magmaDoubleComplex *A, magma_int_t LDA,
        magmaDoubleComplex *B, magma_int_t LDB,
        CUstream stream)
{
    int i, j;
    magmaDoubleComplex *cola, *colb;
    cublasHandle_t handle;
    cublasStatus_t stat;

    stat = cublasCreate(&handle);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        printf ("CUBLAS initialization failed\n");
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    stat = cublasSetStream(handle, stream);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        printf ("cublasSetStream failed\n");
        assert( stat == CUBLAS_STATUS_SUCCESS );
    }

    if (M < 0) {
        return -1;
    }
    if (N < 0) {
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        return -7;
    }

    if (side == MagmaLeft){
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(handle, i+1, cola, 1, colb, 1);
            cudaMemcpyAsync(colb , cola,
                            (i+1)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }else{
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(handle, M-i, cola + i, 1, colb + i, 1);
            cudaMemcpyAsync(colb+i , cola+i,
                            (M-i)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }

    cublasDestroy(handle);

    return MORSE_SUCCESS;
}
#else /* CHAMELEON_USE_CUBLAS_V2 */
int CUDA_zgemerge(
        magma_side_t side, magma_diag_t diag,
        magma_int_t M, magma_int_t N,
        magmaDoubleComplex *A, magma_int_t LDA,
        magmaDoubleComplex *B, magma_int_t LDB,
        CUstream stream)
{
    int i, j;
    magmaDoubleComplex *cola, *colb;

    if (M < 0) {
        return -1;
    }
    if (N < 0) {
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        return -7;
    }

    if (side == MagmaLeft){
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(i+1, cola, 1, colb, 1);
            cudaMemcpyAsync(colb , cola,
                            (i+1)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }else{
        for(i=0; i<N; i++){
            cola = A + i*LDA;
            colb = B + i*LDB;
//            cublasZcopy(M-i, cola + i, 1, colb + i, 1);
            cudaMemcpyAsync(colb+i , cola+i,
                            (M-i)*sizeof(cuDoubleComplex),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }

    return MORSE_SUCCESS;
}
#endif /* CHAMELEON_USE_CUBLAS_V2 */
#endif
