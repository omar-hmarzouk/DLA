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
 * @file cuda_zgelqt.c
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
        CUstream stream)
{
#define da_ref(a_1,a_2) ( da+(a_2)*(ldda) + (a_1))
#define v_ref(a_1,a_2)  ( v+(a_2)*(ldv) + (a_1))
#define dt_ref(a_1,a_2) ( dt+(a_2)*(lddt) + (a_1))
#define t_ref(a_1,a_2)  ( t+(a_2)*(ldt) + (a_1))

    int i, k, ib, lddwork, old_i, old_ib, rows, cols;
    double _Complex one=1.;

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

    lddwork= m;

    /* lower parts of little T must be zero: memset to 0 for simplicity */
    memset(t_ref(0,0), 0, nb*n*sizeof(magmaDoubleComplex));
    cudaMemset(dt_ref(0,0), 0, nb*n*sizeof(magmaDoubleComplex));

    /* copy first panel of A on the host */
    cublasGetMatrix(min(m, nb), n, sizeof(magmaDoubleComplex),
              da_ref(0, 0), ldda,
              v, ldv);

    /* Use blocked code initially */
    for (i = 0; i < k; i += nb) {

    ib = min(k-i, nb);
    if (i+nb >= m) ib = min(m-i, nb);
    cols = n-i;

    if (i > 0){

      /* copy panel of A from device to host */
      cublasGetMatrix(ib, n, sizeof(magmaDoubleComplex),
                      da_ref(i, 0), ldda,
                      v, ldv);

      /* Apply H' to A(i+2*ib:m, i:n) from the right */
      rows = m-old_i-2*old_ib;
      if (rows > 0){
          magma_zlarfb_gpu( MagmaRight, MagmaNoTrans, MagmaForward, MagmaRowwise,
                            rows, n-old_i, old_ib,
                            da_ref(old_i, old_i), ldda, dt_ref(0,old_i), lddt,
                            da_ref(old_i+2*old_ib, old_i), ldda,
                            dwork, lddwork);
      }

      /* copy the lower diag tile into d_A */
      CUDA_zgemerge(MagmaRight, MagmaUnit, old_ib, old_ib,
                    dd, ldd, da_ref(old_i, old_i), ldda, stream);

    }

    /* Form the triangular factor of the block reflector on the host
    H = H'(i+ib-1) . . . H(i+1) H(i) */
    CORE_zgelqt(ib, cols, ib,
              (double _Complex*) v_ref(0,i), ldv,
              (double _Complex*) t_ref(0,0), ldt,
              (double _Complex*) tau+i,
              (double _Complex*) hwork);

    if ( i + ib < m ){
      /* put 0s in the lower triangular part of a panel (and 1s on the
       diagonal); copy the lower triangular in d */
      CORE_zgesplit(MorseRight, MorseUnit, ib, min(cols,ib),
                    (double _Complex*) v_ref(0,i), ldv,
                    (double _Complex*) d, ldd);

      /* copy from host to device a tile diag */
      cublasSetMatrix( ib, min(cols,ib), sizeof(magmaDoubleComplex),
                       d, ldd, dd, ldd );
    }

    /* Send the triangular factor T to the GPU */
    cublasSetMatrix( ib, ib, sizeof(magmaDoubleComplex),
                   t_ref(0,0), ldt, dt_ref(0,i), lddt );

    /* A panel (with zeros in lower tri of its diag) is ready to be used
    in input of zlarfb_gpu: we send the panel to the gpu */
    cublasSetMatrix( ib, cols, sizeof(magmaDoubleComplex),
                   v_ref(0,i), ldv, da_ref(i,i), ldda );

    if (i + ib < m) {

      if (i+2*ib < m){
          rows = ib;
      }
      else{
          rows = m-i-ib;
      }
      /* Apply H' to A(i+ib:i+2*ib, i:n) from the right */
      magma_zlarfb_gpu( MagmaRight, MagmaNoTrans, MagmaForward, MagmaRowwise,
                        rows, cols, ib, da_ref(i,i), ldda, dt_ref(0,i),
                        lddt, da_ref(i+ib,i), ldda, dwork, lddwork);
      old_i = i;
      old_ib = ib;
      if (i+nb >= k){
          /* Apply H' to A(i+2*ib:m, i:n) from the right */
          rows = m-old_i-2*old_ib;
          if (rows > 0){
              magma_zlarfb_gpu( MagmaRight, MagmaNoTrans, MagmaForward, MagmaRowwise,
                                rows, cols, old_ib,
                                da_ref(old_i, old_i), ldda, dt_ref(0,old_i), lddt,
                                da_ref(old_i+2*old_ib, old_i), ldda,
                                dwork, lddwork);
          }
          /* copy the upper diag tile into d_A */
          CUDA_zgemerge(MagmaRight, MagmaUnit, old_ib, old_ib,
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
#endif
