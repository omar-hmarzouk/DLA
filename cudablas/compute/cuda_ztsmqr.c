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
 * @file cuda_ztsmqr.c
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
        CUstream stream)
{
    int i, i1, i3;
    int NQ, NW;
    int kb;
    int ic = 0;
    int jc = 0;
    int mi = M1;
    int ni = N1;

    /* Check input arguments */
    if ((side != MagmaLeft) && (side != MagmaRight)) {
        return -1;
    }

    /* NQ is the order of Q */
    if (side == MagmaLeft) {
        NQ = M2;
        NW = IB;
    }
    else {
        NQ = N2;
        NW = M1;
    }

    if ((trans != MagmaNoTrans) && (trans != MagmaConjTrans)) {
        return -2;
    }
    if (M1 < 0) {
        return -3;
    }
    if (N1 < 0) {
        return -4;
    }
    if ( (M2 < 0) ||
         ( (M2 != M1) && (side == MagmaRight) ) ){
        return -5;
    }
    if ( (N2 < 0) ||
         ( (N2 != N1) && (side == MagmaLeft) ) ){
        return -6;
    }
    if ((K < 0) ||
        ( (side == MagmaLeft)  && (K > M1) ) ||
        ( (side == MagmaRight) && (K > N1) ) ) {
        return -7;
    }
    if (IB < 0) {
        return -8;
    }
    if (LDA1 < max(1,M1)){
        return -10;
    }
    if (LDA2 < max(1,M2)){
        return -12;
    }
    if (LDV < max(1,NQ)){
        return -14;
    }
    if (LDT < max(1,IB)){
        return -16;
    }
    if (LDWORK < max(1,NW)){
        return -18;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0) || (IB == 0))
        return MAGMA_SUCCESS;

    if (((side == MagmaLeft)  && (trans != MagmaNoTrans))
        || ((side == MagmaRight) && (trans == MagmaNoTrans))) {
        i1 = 0;
        i3 = IB;
    }
    else {
        i1 = ((K-1) / IB)*IB;
        i3 = -IB;
    }

    for(i = i1; (i > -1) && (i < K); i += i3) {
        kb = min(IB, K-i);

        if (side == MagmaLeft) {
            /*
             * H or H' is applied to C(i:m,1:n)
             */
            mi = M1 - i;
            ic = i;
        }
        else {
            /*
             * H or H' is applied to C(1:m,i:n)
             */
            ni = N1 - i;
            jc = i;
        }
        /*
         * Apply H or H' (NOTE: CORE_zparfb used to be CORE_ztsrfb)
         */
        CUDA_zparfb(
                side, trans, MagmaForward, MagmaColumnwise,
                mi, ni, M2, N2, kb, 0,
                A1 + LDA1*jc+ic, LDA1,
                A2, LDA2,
                V + LDV*i, LDV,
                T + LDT*i, LDT,
                WORK, LDWORK, WORKC, LDWORKC, stream );
    }
    return MORSE_SUCCESS;
}
#endif
