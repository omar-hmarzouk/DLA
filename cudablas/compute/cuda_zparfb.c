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
 * @file cuda_zparfb.c
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

int
CUDA_zparfb(MORSE_enum side, MORSE_enum trans,
            MORSE_enum direct, MORSE_enum storev,
            int M1, int N1, int M2, int N2, int K, int L,
                  cuDoubleComplex *A1, int LDA1,
                  cuDoubleComplex *A2, int LDA2,
            const cuDoubleComplex *V, int LDV,
            const cuDoubleComplex *T, int LDT,
                  cuDoubleComplex *WORK, int LDWORK,
                  cuDoubleComplex *WORKC, int LDWORKC,
            CUBLAS_STREAM_PARAM )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex zzero = make_cuDoubleComplex(0.0, 0.0);
    cuDoubleComplex zone  = make_cuDoubleComplex(1.0, 0.0);
    cuDoubleComplex mzone = make_cuDoubleComplex(-1.0, 0.0);
#else
    double zzero = 0.0;
    double zone  = 1.0;
    double mzone = -1.0;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    int j;
    MORSE_enum transW;
    MORSE_enum transA2;

    CUBLAS_GET_STREAM;

    /* Check input arguments */
    if ((side != MorseLeft) && (side != MorseRight)) {
        return -1;
    }
    if ((trans != MorseNoTrans) && (trans != MorseConjTrans)) {
        return -2;
    }
    if ((direct != MorseForward) && (direct != MorseBackward)) {
        return -3;
    }
    if ((storev != MorseColumnwise) && (storev != MorseRowwise)) {
        return -4;
    }
    if (M1 < 0) {
        return -5;
    }
    if (N1 < 0) {
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == MorseRight) && (M1 != M2) ) ) {
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == MorseLeft) && (N1 != N2) ) ) {
        return -8;
    }
    if (K < 0) {
        return -9;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
        return MORSE_SUCCESS;

    if (direct == MorseForward) {

        if (side == MorseLeft) {

            /*
             * Column or Rowwise / Forward / Left
             * ----------------------------------
             *
             * Form  H * A  or  H' * A  where  A = ( A1 )
             *                                     ( A2 )
             */

            /*
             * W = A1 + V' * A2:
             *      W = A1
             *      W = W + V' * A2
             *
             */
            cudaMemcpy2DAsync( WORK, LDWORK * sizeof(cuDoubleComplex),
                               A1,   LDA1   * sizeof(cuDoubleComplex),
                               K * sizeof(cuDoubleComplex), N1,
                               cudaMemcpyDeviceToDevice, stream );

            transW  = storev == MorseColumnwise ? MorseConjTrans : MorseNoTrans;
            transA2 = storev == MorseColumnwise ? MorseNoTrans : MorseConjTrans;

            cublasZgemm(CUBLAS_HANDLE
                        morse_lapack_const(transW), morse_lapack_const(MorseNoTrans),
                        K, N1, M2,
                        CUBLAS_SADDR(zone),
                        V     /* K*M2  */, LDV,
                        A2    /* M2*N1 */, LDA2,
                        CUBLAS_SADDR(zone),
                        WORK  /* K*N1  */, LDWORK);

            if (WORKC == NULL) {
                /* W = op(T) * W */
                cublasZtrmm( CUBLAS_HANDLE
                             morse_lapack_const(MorseLeft), morse_lapack_const(MorseUpper),
                             morse_lapack_const(trans), morse_lapack_const(MorseNonUnit),
                             K, N2,
                             CUBLAS_SADDR(zone),
                             T,    LDT,
                             WORK, LDWORK);


                /* A1 = A1 - W = A1 - op(T) * W */
                for(j = 0; j < N1; j++) {
                    cublasZaxpy(CUBLAS_HANDLE
                                K, CUBLAS_SADDR(mzone),
                                (WORK + LDWORK*j), 1,
                                (A1 + LDA1*j),     1);
                }

                /* A2 = A2 - op(V) * W  */
                cublasZgemm(CUBLAS_HANDLE
                            morse_lapack_const(transA2), morse_lapack_const(MorseNoTrans),
                            M2, N2, K,
                            CUBLAS_SADDR(mzone), V    /* M2*K  */, LDV,
                                                 WORK /* K*N2  */, LDWORK,
                            CUBLAS_SADDR(zone),  A2   /* m2*N2 */, LDA2);

            } else {
                /* Wc = V * op(T) */
                cublasZgemm( CUBLAS_HANDLE
                             morse_lapack_const(transA2), morse_lapack_const(trans),
                             M2, K, K,
                             CUBLAS_SADDR(zone),  V, LDV,
                                                  T, LDT,
                             CUBLAS_SADDR(zzero), WORKC, LDWORKC );

                /* A1 = A1 - opt(T) * W */
                cublasZgemm( CUBLAS_HANDLE
                             morse_lapack_const(trans), morse_lapack_const(MorseNoTrans),
                             K, N1, K,
                             CUBLAS_SADDR(mzone), T,    LDT,
                                                  WORK, LDWORK,
                             CUBLAS_SADDR(zone),  A1,   LDA1 );

                /* A2 = A2 - Wc * W */
                cublasZgemm( CUBLAS_HANDLE
                             morse_lapack_const(MorseNoTrans), morse_lapack_const(MorseNoTrans),
                             M2, N2, K,
                             CUBLAS_SADDR(mzone), WORKC, LDWORKC,
                                                  WORK,  LDWORK,
                             CUBLAS_SADDR(zone),  A2,    LDA2 );
            }
        }
        else {
            /*
             * Column or Rowwise / Forward / Right
             * -----------------------------------
             *
             * Form  H * A  or  H' * A  where A  = ( A1 A2 )
             *
             */

            /*
             * W = A1 + A2 * V':
             *      W = A1
             *      W = W + A2 * V'
             *
             */
            cudaMemcpy2DAsync( WORK, LDWORK * sizeof(cuDoubleComplex),
                               A1,   LDA1   * sizeof(cuDoubleComplex),
                               M1 * sizeof(cuDoubleComplex), K,
                               cudaMemcpyDeviceToDevice, stream );

            transW  = storev == MorseColumnwise ? MorseNoTrans : MorseConjTrans;
            transA2 = storev == MorseColumnwise ? MorseConjTrans : MorseNoTrans;

            cublasZgemm(CUBLAS_HANDLE
                        morse_lapack_const(MorseNoTrans), morse_lapack_const(transW),
                        M1, K, N2,
                        CUBLAS_SADDR(zone), A2   /* M1*N2 */, LDA2,
                                            V    /* N2*K  */, LDV,
                        CUBLAS_SADDR(zone), WORK /* M1*K  */, LDWORK);

            if (WORKC == NULL) {
                /* W = W * op(T) */
                cublasZtrmm( CUBLAS_HANDLE
                             morse_lapack_const(MorseRight), morse_lapack_const(MorseUpper),
                             morse_lapack_const(trans), morse_lapack_const(MorseNonUnit),
                             M2, K,
                             CUBLAS_SADDR(zone),
                             T,    LDT,
                             WORK, LDWORK);


                /* A1 = A1 - W = A1 - W * op(T) */
                for(j = 0; j < K; j++) {
                    cublasZaxpy(CUBLAS_HANDLE
                                M1, CUBLAS_SADDR(mzone),
                                (WORK + LDWORK*j), 1,
                                (A1 + LDA1*j), 1);
                }

                /* A2 = A2 - W * op(V)  */
                cublasZgemm(CUBLAS_HANDLE
                            morse_lapack_const(MorseNoTrans), morse_lapack_const(transA2),
                            M2, N2, K,
                            CUBLAS_SADDR(mzone), WORK /* M2*K  */, LDWORK,
                                                 V    /* K*N2  */, LDV,
                            CUBLAS_SADDR(zone),  A2   /* M2*N2 */, LDA2);

            } else {
                /* A1 = A1 - W * opt(T) */
                cublasZgemm( CUBLAS_HANDLE
                             morse_lapack_const(MorseNoTrans), morse_lapack_const(trans),
                             M1, K, K,
                             CUBLAS_SADDR(mzone), WORK, LDWORK,
                                                  T,    LDT,
                             CUBLAS_SADDR(zone),  A1,   LDA1 );

                /* Wc = op(T) * V */
                cublasZgemm( CUBLAS_HANDLE
                             morse_lapack_const(trans), morse_lapack_const(transA2),
                             K, N2, K,
                             CUBLAS_SADDR(zone),  T,     LDT,
                                                  V,     LDV,
                             CUBLAS_SADDR(zzero), WORKC, LDWORKC );

                /* A2 = A2 - W * Wc */
                cublasZgemm( CUBLAS_HANDLE
                             morse_lapack_const(MorseNoTrans), morse_lapack_const(MorseNoTrans),
                             M2, N2, K,
                             CUBLAS_SADDR(mzone), WORK,  LDWORK,
                                                  WORKC, LDWORKC,
                             CUBLAS_SADDR(zone),  A2,    LDA2 );
            }
        }
    }
    else {
        fprintf(stderr, "Not implemented (Backward / Left or Right)");
        return MORSE_ERR_NOT_SUPPORTED;
    }

    return MORSE_SUCCESS;
}
