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

#if defined(CHAMELEON_USE_MAGMA)
#if defined(CHAMELEON_USE_CUBLAS_V2)
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
        CUstream stream)

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
    magma_trans_t transW;
    magma_trans_t transA2;
    cublasHandle_t handle;
    cublasStatus_t stat;
    cublasOperation_t cublasTrans;
    cublasOperation_t cublasTransW;
    cublasOperation_t cublasTransA2;

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

    if (trans == MagmaNoTrans){
        cublasTrans = CUBLAS_OP_N;
    }else if(trans == MagmaTrans){
        cublasTrans = CUBLAS_OP_T;
    }else if(trans == MagmaConjTrans){
        cublasTrans = CUBLAS_OP_C;
    }else{
        fprintf(stderr, "Error in CUDA_zparfb: bad trans parameter %d\n", trans);
    }

    /* Check input arguments */
    if ((side != MagmaLeft) && (side != MagmaRight)) {
        return -1;
    }
    if ((trans != MagmaNoTrans) && (trans != MagmaConjTrans)) {
        return -2;
    }
    if ((direct != MagmaForward) && (direct != MagmaBackward)) {
        return -3;
    }
    if ((storev != MagmaColumnwise) && (storev != MagmaRowwise)) {
        return -4;
    }
    if (M1 < 0) {
        return -5;
    }
    if (N1 < 0) {
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == MagmaRight) && (M1 != M2) ) ) {
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == MagmaLeft) && (N1 != N2) ) ) {
        return -8;
    }
    if (K < 0) {
        return -9;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
        return MAGMA_SUCCESS;

    if (direct == MagmaForward) {

        if (side == MagmaLeft) {

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

            transW  = storev == MorseColumnwise ? MagmaConjTrans : MagmaNoTrans;
            transA2 = storev == MorseColumnwise ? MagmaNoTrans : MagmaConjTrans;

            if (transW == MagmaNoTrans){
                cublasTransW = CUBLAS_OP_N;
            }else if(transW == MagmaTrans){
                cublasTransW = CUBLAS_OP_T;
            }else if(transW == MagmaConjTrans){
                cublasTransW = CUBLAS_OP_C;
            }else{
                fprintf(stderr, "Error in CUDA_zparfb: bad transW parameter %d\n", transW);
            }
            if (transA2 == MagmaNoTrans){
                cublasTransA2 = CUBLAS_OP_N;
            }else if(transA2 == MagmaTrans){
                cublasTransA2 = CUBLAS_OP_T;
            }else if(transA2 == MagmaConjTrans){
                cublasTransA2 = CUBLAS_OP_C;
            }else{
                fprintf(stderr, "Error in CUDA_zparfb: bad transA2 parameter %d\n", transA2);
            }

            cublasZgemm(handle, cublasTransW, CUBLAS_OP_N,
                        K, N1, M2,
                        (const cuDoubleComplex *) &zone,
                        (const cuDoubleComplex*)V     /* K*M2  */, LDV,
                        (const cuDoubleComplex*)A2    /* M2*N1 */, LDA2,
                        (const cuDoubleComplex *) &zone,
                        (cuDoubleComplex*)WORK  /* K*N1  */, LDWORK);

            if (WORKC == NULL) {
                /* W = op(T) * W */
                cublasZtrmm( handle,
                    CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER,
                    cublasTrans, CUBLAS_DIAG_NON_UNIT,
                    K, N2,
                    (const cuDoubleComplex *) &zone,
                    (const cuDoubleComplex*)T, LDT,
                    (cuDoubleComplex*)WORK, LDWORK,
                    (cuDoubleComplex*)WORK, LDWORK);


                /* A1 = A1 - W = A1 - op(T) * W */
                for(j = 0; j < N1; j++) {
                    cublasZaxpy(handle, K, (const cuDoubleComplex *) &mzone,
                                (const cuDoubleComplex*)(WORK + LDWORK*j), 1,
                                (cuDoubleComplex*)(A1 + LDA1*j), 1);
                }

                /* A2 = A2 - op(V) * W  */
                cublasZgemm(handle, cublasTransA2, CUBLAS_OP_N,
                            M2, N2, K,
                            (const cuDoubleComplex *) &mzone,
                            (const cuDoubleComplex*)V     /* M2*K  */, LDV,
                            (const cuDoubleComplex*)WORK  /* K*N2  */, LDWORK,
                            (const cuDoubleComplex *) &zone,
                            (cuDoubleComplex*)A2    /* m2*N2 */, LDA2);

            } else {
                /* Wc = V * op(T) */
                cublasZgemm( handle, cublasTransA2, cublasTrans,
                             M2, K, K,
                             (const cuDoubleComplex *) &zone,  V,     LDV,
                                    T,     LDT,
                             (const cuDoubleComplex *) &zzero, WORKC, LDWORKC );

                /* A1 = A1 - opt(T) * W */
                cublasZgemm( handle, cublasTrans, CUBLAS_OP_N,
                             K, N1, K,
                             (const cuDoubleComplex *) &mzone,
                             (const cuDoubleComplex *)T,    LDT,
                             (const cuDoubleComplex *)WORK, LDWORK,
                             (const cuDoubleComplex *) &zone,
                             (cuDoubleComplex*)A1,   LDA1 );

                /* A2 = A2 - Wc * W */
                cublasZgemm( handle, CUBLAS_OP_N, CUBLAS_OP_N,
                             M2, N2, K,
                             (const cuDoubleComplex *) &mzone,
                             (const cuDoubleComplex *)WORKC, LDWORKC,
                             (const cuDoubleComplex *)WORK,  LDWORK,
                             (const cuDoubleComplex *) &zone,
                             (cuDoubleComplex *)A2,    LDA2 );
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

            transW  = storev == MorseColumnwise ? MagmaNoTrans : MagmaConjTrans;
            transA2 = storev == MorseColumnwise ? MagmaConjTrans : MagmaNoTrans;

            if (transW == MagmaNoTrans){
                cublasTransW = CUBLAS_OP_N;
            }else if(transW == MagmaTrans){
                cublasTransW = CUBLAS_OP_T;
            }else if(transW == MagmaConjTrans){
                cublasTransW = CUBLAS_OP_C;
            }else{
                fprintf(stderr, "Error in CUDA_zparfb: bad transW parameter %d\n", transW);
            }
            if (transA2 == MagmaNoTrans){
                cublasTransA2 = CUBLAS_OP_N;
            }else if(transA2 == MagmaTrans){
                cublasTransA2 = CUBLAS_OP_T;
            }else if(transA2 == MagmaConjTrans){
                cublasTransA2 = CUBLAS_OP_C;
            }else{
                fprintf(stderr, "Error in CUDA_zparfb: bad transA2 parameter %d\n", transA2);
            }

            cublasZgemm(handle, CUBLAS_OP_N, cublasTransW,
                        M1, K, N2,
                        (const cuDoubleComplex *) &zone,
                        (const cuDoubleComplex*)A2    /* M1*N2 */, LDA2,
                        (const cuDoubleComplex*)V     /* N2*K  */, LDV,
                        (const cuDoubleComplex *) &zone,
                        (cuDoubleComplex*)WORK  /* M1*K  */, LDWORK);

            if (WORKC == NULL) {
                /* W = W * op(T) */
                cublasZtrmm( handle,
                    CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER,
                    cublasTrans, CUBLAS_DIAG_NON_UNIT,
                    M2, K,
                    (const cuDoubleComplex *) &zone,
                    (const cuDoubleComplex*)T, LDT,
                    (cuDoubleComplex*)WORK, LDWORK,
                    (cuDoubleComplex*)WORK, LDWORK);


                /* A1 = A1 - W = A1 - W * op(T) */
                for(j = 0; j < K; j++) {
                    cublasZaxpy(handle, M1, (const cuDoubleComplex *) &mzone,
                        (const cuDoubleComplex*)(WORK + LDWORK*j), 1,
                        (cuDoubleComplex*)(A1 + LDA1*j), 1);
                }

                /* A2 = A2 - W * op(V)  */
                cublasZgemm(handle, CUBLAS_OP_N, cublasTransA2,
                            M2, N2, K,
                            (const cuDoubleComplex *) &mzone,
                            (const cuDoubleComplex*)WORK  /* M2*K  */, LDWORK,
                            (const cuDoubleComplex*)V     /* K*N2  */, LDV,
                            (const cuDoubleComplex *) &zone,
                            (cuDoubleComplex*)A2    /* M2*N2 */, LDA2);

            } else {
                /* A1 = A1 - W * opt(T) */
                cublasZgemm( handle, CUBLAS_OP_N, cublasTrans,
                    M1, K, K,
                    (const cuDoubleComplex *) &mzone,
                    (const cuDoubleComplex *)WORK, LDWORK,
                    (const cuDoubleComplex *)T,    LDT,
                    (const cuDoubleComplex *) &zone,
                    (cuDoubleComplex *)A1,   LDA1 );

                /* Wc = op(T) * V */
                cublasZgemm( handle, cublasTrans, cublasTransA2,
                             K, N2, K,
                             (const cuDoubleComplex *) &zone,
                             (const cuDoubleComplex *)T,     LDT,
                             (const cuDoubleComplex *)V,     LDV,
                             (const cuDoubleComplex *) &zzero,
                             (cuDoubleComplex *)WORKC, LDWORKC );

                /* A2 = A2 - W * Wc */
                cublasZgemm( handle, CUBLAS_OP_N, CUBLAS_OP_N,
                             M2, N2, K,
                             (const cuDoubleComplex *) &mzone,
                             (const cuDoubleComplex *)WORK,  LDWORK,
                             (const cuDoubleComplex *)WORKC, LDWORKC,
                             (const cuDoubleComplex *) &zone,
                             (cuDoubleComplex *)A2,    LDA2 );
            }
        }
    }
    else {
        fprintf(stderr, "Not implemented (Backward / Left or Right)");
        return MAGMA_ERR_NOT_SUPPORTED;
    }

    cublasDestroy(handle);

    return MORSE_SUCCESS;
}
#else /* CHAMELEON_USE_CUBLAS_V2 */
magma_int_t
CUDA_zparfb(magma_side_t side, magma_trans_t trans,
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
                       CUstream stream)
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
    magma_trans_t transW;
    magma_trans_t transA2;

    /* Check input arguments */
    if ((side != MagmaLeft) && (side != MagmaRight)) {
        return -1;
    }
    if ((trans != MagmaNoTrans) && (trans != MagmaConjTrans)) {
        return -2;
    }
    if ((direct != MagmaForward) && (direct != MagmaBackward)) {
        return -3;
    }
    if ((storev != MagmaColumnwise) && (storev != MagmaRowwise)) {
        return -4;
    }
    if (M1 < 0) {
        return -5;
    }
    if (N1 < 0) {
        return -6;
    }
    if ((M2 < 0) ||
        ( (side == MagmaRight) && (M1 != M2) ) ) {
        return -7;
    }
    if ((N2 < 0) ||
        ( (side == MagmaLeft) && (N1 != N2) ) ) {
        return -8;
    }
    if (K < 0) {
        return -9;
    }

    /* Quick return */
    if ((M1 == 0) || (N1 == 0) || (M2 == 0) || (N2 == 0) || (K == 0))
        return MAGMA_SUCCESS;

    if (direct == MagmaForward) {

        if (side == MagmaLeft) {

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

            transW  = storev == MorseColumnwise ? MagmaConjTrans : MagmaNoTrans;
            transA2 = storev == MorseColumnwise ? MagmaNoTrans : MagmaConjTrans;

            cublasZgemm(morse_lapack_const(transW), morse_lapack_const(MagmaNoTrans),
                        K, N1, M2,
                        zone,
                        (cuDoubleComplex*)V     /* K*M2  */, LDV,
                        (cuDoubleComplex*)A2    /* M2*N1 */, LDA2,
                        zone,
                        (cuDoubleComplex*)WORK  /* K*N1  */, LDWORK);

            if (WORKC == NULL) {
                /* W = op(T) * W */
                cublasZtrmm( morse_lapack_const(MagmaLeft), morse_lapack_const(MagmaUpper),
                             morse_lapack_const(trans), morse_lapack_const(MagmaNonUnit),
                             K, N2,
                             zone,
                             (cuDoubleComplex*)T, LDT,
                             (cuDoubleComplex*)WORK, LDWORK);


                /* A1 = A1 - W = A1 - op(T) * W */
                for(j = 0; j < N1; j++) {
                    cublasZaxpy(K, mzone,
                                (cuDoubleComplex*)(WORK + LDWORK*j), 1,
                                (cuDoubleComplex*)(A1 + LDA1*j), 1);
                }

                /* A2 = A2 - op(V) * W  */
                cublasZgemm(morse_lapack_const(transA2), morse_lapack_const(MagmaNoTrans),
                            M2, N2, K,
                            mzone,
                            (cuDoubleComplex*)V     /* M2*K  */, LDV,
                            (cuDoubleComplex*)WORK  /* K*N2  */, LDWORK,
                            zone,
                            (cuDoubleComplex*)A2    /* m2*N2 */, LDA2);

            } else {
                /* Wc = V * op(T) */
                cublasZgemm( morse_lapack_const(transA2), morse_lapack_const(trans),
                             M2, K, K,
                             zone,  V,     LDV,
                                    T,     LDT,
                             zzero, WORKC, LDWORKC );

                /* A1 = A1 - opt(T) * W */
                cublasZgemm( morse_lapack_const(trans), morse_lapack_const(MagmaNoTrans),
                             K, N1, K,
                             mzone, T,    LDT,
                                    WORK, LDWORK,
                             zone,  A1,   LDA1 );

                /* A2 = A2 - Wc * W */
                cublasZgemm( morse_lapack_const(MagmaNoTrans), morse_lapack_const(MagmaNoTrans),
                             M2, N2, K,
                             mzone, WORKC, LDWORKC,
                                    WORK,  LDWORK,
                             zone,  A2,    LDA2 );
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

            transW  = storev == MorseColumnwise ? MagmaNoTrans : MagmaConjTrans;
            transA2 = storev == MorseColumnwise ? MagmaConjTrans : MagmaNoTrans;

            cublasZgemm(morse_lapack_const(MagmaNoTrans), morse_lapack_const(transW),
                        M1, K, N2,
                        zone,
                        (cuDoubleComplex*)A2    /* M1*N2 */, LDA2,
                        (cuDoubleComplex*)V     /* N2*K  */, LDV,
                        zone,
                        (cuDoubleComplex*)WORK  /* M1*K  */, LDWORK);

            if (WORKC == NULL) {
                /* W = W * op(T) */
                cublasZtrmm( morse_lapack_const(MagmaRight), morse_lapack_const(MagmaUpper),
                             morse_lapack_const(trans), morse_lapack_const(MagmaNonUnit),
                             M2, K,
                             zone,
                             (cuDoubleComplex*)T, LDT,
                             (cuDoubleComplex*)WORK, LDWORK);


                /* A1 = A1 - W = A1 - W * op(T) */
                for(j = 0; j < K; j++) {
                    cublasZaxpy(M1, mzone,
                                (cuDoubleComplex*)(WORK + LDWORK*j), 1,
                                (cuDoubleComplex*)(A1 + LDA1*j), 1);
                }

                /* A2 = A2 - W * op(V)  */
                cublasZgemm(morse_lapack_const(MagmaNoTrans), morse_lapack_const(transA2),
                            M2, N2, K,
                            mzone,
                            (cuDoubleComplex*)WORK  /* M2*K  */, LDWORK,
                            (cuDoubleComplex*)V     /* K*N2  */, LDV,
                            zone,
                            (cuDoubleComplex*)A2    /* M2*N2 */, LDA2);

            } else {
                /* A1 = A1 - W * opt(T) */
                cublasZgemm( morse_lapack_const(MagmaNoTrans), morse_lapack_const(trans),
                             M1, K, K,
                             mzone, WORK, LDWORK,
                                    T,    LDT,
                             zone,  A1,   LDA1 );

                /* Wc = op(T) * V */
                cublasZgemm( morse_lapack_const(trans), morse_lapack_const(transA2),
                             K, N2, K,
                             zone,  T,     LDT,
                                    V,     LDV,
                             zzero, WORKC, LDWORKC );

                /* A2 = A2 - W * Wc */
                cublasZgemm( morse_lapack_const(MagmaNoTrans), morse_lapack_const(MagmaNoTrans),
                             M2, N2, K,
                             mzone, WORK,  LDWORK,
                                    WORKC, LDWORKC,
                             zone,  A2,    LDA2 );
            }
        }
    }
    else {
        fprintf(stderr, "Not implemented (Backward / Left or Right)");
        return MAGMA_ERR_NOT_SUPPORTED;
    }

    return MORSE_SUCCESS;
}
#endif /* CHAMELEON_USE_CUBLAS_V2 */
#endif
