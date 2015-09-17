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
 * @file cuda_zunmlqt.c
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
int CUDA_zunmlqt(
        magma_side_t side, magma_trans_t trans,
        magma_int_t M, magma_int_t N, magma_int_t K, magma_int_t IB,
        const magmaDoubleComplex *A,    magma_int_t LDA,
        const magmaDoubleComplex *T,    magma_int_t LDT,
        magmaDoubleComplex *C,    magma_int_t LDC,
        magmaDoubleComplex *WORK, magma_int_t LDWORK )
{
    int i, kb;
    int i1, i3;
    int nq, nw;
    int ic = 0;
    int jc = 0;
    int ni = N;
    int mi = M;

    /* Check input arguments */
    if ((side != MagmaLeft) && (side != MagmaRight)) {
        return -1;
    }
    /*
     * NQ is the order of Q and NW is the minimum dimension of WORK
     */
    if (side == MagmaLeft) {
        nq = M;
        nw = N;
    }
    else {
        nq = N;
        nw = M;
    }

    if ((trans != MagmaNoTrans) && (trans != MagmaConjTrans)) {
        return -2;
    }
    if (M < 0) {
        return -3;
    }
    if (N < 0) {
        return -4;
    }
    if ((K < 0) || (K > nq)) {
        return -5;
    }
    if ((IB < 0) || ( (IB == 0) && ((M > 0) && (N > 0)) )) {
        return -6;
    }
    if ((LDA < max(1,K)) && (K > 0)) {
        return -8;
    }
    if ((LDC < max(1,M)) && (M > 0)) {
        return -12;
    }
    if ((LDWORK < max(1,nw)) && (nw > 0)) {
        return -14;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (K == 0))
        return MAGMA_SUCCESS;

    if (((side == MagmaLeft) && (trans == MagmaNoTrans))
        || ((side == MagmaRight) && (trans != MagmaNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ( ( K-1 ) / IB )*IB;
        i3 = -IB;
    }

    if( trans == MorseNoTrans) {
        trans = MorseConjTrans;
    }
    else {
        trans = MorseNoTrans;
    }

    for(i = i1; (i >- 1) && (i < K); i+=i3 ) {
        kb = min(IB, K-i);

        if (side == MagmaLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N - i;
            jc = i;
        }

        magma_zlarfb_gpu( side, trans, MagmaForward, MagmaRowwise,
                          mi, ni, kb,
                          A + LDA * i  + i,  LDA,
                          T + LDT * i,       LDT,
                          C + LDC * jc + ic, LDC,
                          WORK, LDWORK);
    }
    return MORSE_SUCCESS;
}
#endif
