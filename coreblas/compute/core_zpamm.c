/**
 *
 * @file core_zpamm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2011-06-14
 * @precisions normal z -> c d s
 *
 **/
#include "coreblas/cblas.h"
#include "coreblas/lapacke.h"
#include "coreblas.h"

static inline int CORE_zpamm_a2(MORSE_enum side, MORSE_enum trans, MORSE_enum uplo,
                                int M, int N, int K, int L,
                                int vi2, int vi3,
                                MORSE_Complex64_t *A2, int LDA2,
                                const MORSE_Complex64_t *V, int LDV,
                                MORSE_Complex64_t *W, int LDW);
static inline int CORE_zpamm_w(MORSE_enum side, MORSE_enum trans, MORSE_enum uplo,
                               int M, int N, int K, int L,
                               int vi2, int vi3,
                               const MORSE_Complex64_t *A1, int LDA1,
                                     MORSE_Complex64_t *A2, int LDA2,
                               const MORSE_Complex64_t *V, int LDV,
                                     MORSE_Complex64_t *W, int LDW);

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  ZPAMM  performs one of the matrix-matrix operations
 *
 *                    LEFT                      RIGHT
 *     OP MorseW  :  W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *     OP MorseA2 :  A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 *  where  op( V ) is one of
 *
 *     op( V ) = V   or   op( V ) = V**T   or   op( V ) = V**H,
 *
 *  A1, A2 and W are general matrices, and V is:
 *
 *        l = k: rectangle + triangle
 *        l < k: rectangle + trapezoid
 *        l = 0: rectangle
 *
 *  Size of V, both rowwise and columnwise, is:
 *
 *         ----------------------
 *          side   trans    size
 *         ----------------------
 *          left     N     M x K
 *                   T     K x M
 *          right    N     K x N
 *                   T     N x K
 *         ----------------------
 *
 *  LEFT (columnwise and rowwise):
 *
 *              |    K    |                 |         M         |
 *           _  __________   _              _______________        _
 *              |    |    |                 |             | \
 *     V:       |    |    |            V':  |_____________|___\    K
 *              |    |    | M-L             |                  |
 *           M  |    |    |                 |__________________|   _
 *              |____|    |  _
 *              \    |    |                 |    M - L    | L  |
 *                \  |    |  L
 *           _      \|____|  _
 *
 *  RIGHT (columnwise and rowwise):
 *
 *          |         K         |                   |    N    |
 *          _______________        _             _  __________   _
 *          |             | \                       |    |    |
 *     V':  |_____________|___\    N        V:      |    |    |
 *          |                  |                    |    |    | K-L
 *          |__________________|   _             K  |    |    |
 *                                                  |____|    |  _
 *          |    K - L    | L  |                    \    |    |
 *                                                    \  |    |  L
 *                                               _      \|____|  _
 *
 *  Arguments
 *  ==========
 *
 * @param[in] op
 *
 *         OP specifies which operation to perform:
 *
 *         @arg MorseW  : W  = A1 + op(V) * A2  or  W  = A1 + A2 * op(V)
 *         @arg MorseA2 : A2 = A2 - op(V) * W   or  A2 = A2 - W * op(V)
 *
 * @param[in] side
 *
 *         SIDE specifies whether  op( V ) multiplies A2
 *         or W from the left or right as follows:
 *
 *         @arg MorseLeft  : multiply op( V ) from the left
 *                            OP MorseW  :  W  = A1 + op(V) * A2
 *                            OP MorseA2 :  A2 = A2 - op(V) * W
 *
 *         @arg MorseRight : multiply op( V ) from the right
 *                            OP MorseW  :  W  = A1 + A2 * op(V)
 *                            OP MorseA2 :  A2 = A2 - W * op(V)
 *
 * @param[in] storev
 *
 *         Indicates how the vectors which define the elementary
 *         reflectors are stored in V:
 *
 *         @arg MorseColumnwise
 *         @arg MorseRowwise
 *
 * @param[in] M
 *         The number of rows of the A1, A2 and W
 *         If SIDE is MorseLeft, the number of rows of op( V )
 *
 * @param[in] N
 *         The number of columns of the A1, A2 and W
 *         If SIDE is MorseRight, the number of columns of op( V )
 *
 * @param[in] K
 *         If SIDE is MorseLeft, the number of columns of op( V )
 *         If SIDE is MorseRight, the number of rows of op( V )
 *
 * @param[in] L
 *         The size of the triangular part of V
 *
 * @param[in] A1
 *         On entry, the M-by-N tile A1.
 *
 * @param[in] LDA1
 *         The leading dimension of the array A1. LDA1 >= max(1,M).
 *
 * @param[in,out] A2
 *         On entry, the M-by-N tile A2.
 *         On exit, if OP is MorseA2 A2 is overwritten
 *
 * @param[in] LDA2
 *         The leading dimension of the tile A2. LDA2 >= max(1,M).
 *
 * @param[in] V
 *         The matrix V as described above.
 *         If SIDE is MorseLeft : op( V ) is M-by-K
 *         If SIDE is MorseRight: op( V ) is K-by-N
 *
 * @param[in] LDV
 *         The leading dimension of the array V.
 *
 * @param[in,out] W
 *         On entry, the M-by-N matrix W.
 *         On exit, W is overwritten either if OP is MorseA2 or MorseW.
 *         If OP is MorseA2, W is an input and is used as a workspace.
 *
 * @param[in] LDW
 *         The leading dimension of array WORK.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int
CORE_zpamm(int op, MORSE_enum side, MORSE_enum storev,
           int M, int N, int K, int L,
           const MORSE_Complex64_t *A1, int LDA1,
                 MORSE_Complex64_t *A2, int LDA2,
           const MORSE_Complex64_t *V, int LDV,
                 MORSE_Complex64_t *W, int LDW)
{


    int vi2, vi3, uplo, trans, info;

    /* Check input arguments */
    if ((op != MorseW) && (op != MorseA2)) {
        coreblas_error(1, "Illegal value of op");
        return -1;
    }
    if ((side != MorseLeft) && (side != MorseRight)) {
        coreblas_error(2, "Illegal value of side");
        return -2;
    }
    if ((storev != MorseColumnwise) && (storev != MorseRowwise)) {
        coreblas_error(3, "Illegal value of storev");
        return -3;
    }
    if (M < 0) {
        coreblas_error(4, "Illegal value of M");
        return -4;
    }
    if (N < 0) {
        coreblas_error(5, "Illegal value of N");
        return -5;
    }
    if (K < 0) {
        coreblas_error(6, "Illegal value of K");
        return -6;
    }
    if (L < 0) {
        coreblas_error(7, "Illegal value of L");
        return -7;
    }
    if (LDA1 < 0) {
        coreblas_error(9, "Illegal value of LDA1");
        return -9;
    }
    if (LDA2 < 0) {
        coreblas_error(11, "Illegal value of LDA2");
        return -11;
    }
    if (LDV < 0) {
        coreblas_error(13, "Illegal value of LDV");
        return -13;
    }
    if (LDW < 0) {
        coreblas_error(15, "Illegal value of LDW");
        return -15;
    }

    /* Quick return */
    if ((M == 0) || (N == 0) || (K == 0))
        return MORSE_SUCCESS;

    /*
     * TRANS is set as:
     *
     *        -------------------------------------
     *         side   direct     MorseW  MorseA2
     *        -------------------------------------
     *         left   colwise       T        N
     *                rowwise       N        T
     *         right  colwise       N        T
     *                rowwise       T        N
     *        -------------------------------------
     */

    /* Columnwise*/
    if (storev == MorseColumnwise) {
        uplo = CblasUpper;
        if (side == MorseLeft) {
            trans = op == MorseA2 ? MorseNoTrans : MorseConjTrans;
            vi2 = trans == MorseNoTrans ? M - L : K - L;
        }
        else {
            trans = op == MorseW ? MorseNoTrans : MorseConjTrans;
            vi2 = trans == MorseNoTrans ? K - L : N - L;
        }
        vi3 = LDV * L;
    }

    /* Rowwise */
    else {
        uplo = CblasLower;
        if (side == MorseLeft) {
            trans = op == MorseW ? MorseNoTrans : MorseConjTrans;
            vi2 = trans == MorseNoTrans ? K - L : M - L;
        }
        else {
            trans = op == MorseA2 ? MorseNoTrans : MorseConjTrans;
            vi2 = trans == MorseNoTrans ? N - L : K - L;
        }
        vi2 *= LDV;
        vi3  = L;
    }

    /**/
    if (op==MorseW) {
        info = CORE_zpamm_w(
                side, trans, uplo, M, N, K, L, vi2, vi3,
                A1, LDA1, A2, LDA2, V, LDV, W, LDW);
        if (info != 0)
            return info;
    } else if (op==MorseA2) {
        info = CORE_zpamm_a2(
                side, trans, uplo, M, N, K, L, vi2, vi3,
                A2, LDA2, V, LDV, W, LDW);
        if (info != 0)
            return info;
    }

    return MORSE_SUCCESS;
}

/***************************************************************************/
static inline int
CORE_zpamm_w(MORSE_enum side, MORSE_enum trans, MORSE_enum uplo,
             int M, int N, int K, int L,
             int vi2, int vi3,
             const MORSE_Complex64_t *A1, int LDA1,
                   MORSE_Complex64_t *A2, int LDA2,
             const MORSE_Complex64_t *V, int LDV,
                   MORSE_Complex64_t *W, int LDW)
{

   /*
    * W = A1 + op(V) * A2  or  W = A1 + A2 * op(V)
    */

    int j;
    static MORSE_Complex64_t zone  =  1.0;
    static MORSE_Complex64_t zzero =  0.0;

    if (side == MorseLeft) {

        if (((trans == MorseConjTrans) && (uplo == CblasUpper)) ||
            ((trans == MorseNoTrans) && (uplo == CblasLower))) {

            /*
             * W = A1 + V' * A2
             */

            /* W = A2_2 */
            LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,
                morse_lapack_const(MorseUpperLower),
                L, N,
                &A2[K-L], LDA2, W, LDW);

            /* W = V_2' * W + V_1' * A2_1 (ge+tr, top L rows of V') */
            if (L > 0) {
                /* W = V_2' * W */
                cblas_ztrmm(
                    CblasColMajor, CblasLeft, (CBLAS_UPLO)uplo,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, L, N,
                    CBLAS_SADDR(zone), &V[vi2], LDV,
                    W, LDW);

                /* W = W + V_1' * A2_1 */
                if (K > L) {
                    cblas_zgemm(
                        CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                        L, N, K-L,
                        CBLAS_SADDR(zone), V, LDV,
                        A2, LDA2,
                        CBLAS_SADDR(zone), W, LDW);
                }
            }

            /* W_2 = V_3' * A2: (ge, bottom M-L rows of V') */
            if (M > L) {
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    (M-L), N, K,
                    CBLAS_SADDR(zone), &V[vi3], LDV,
                    A2, LDA2,
                    CBLAS_SADDR(zzero), &W[L], LDW);
            }

            /* W = A1 + W */
            for(j = 0; j < N; j++) {
                cblas_zaxpy(
                        M, CBLAS_SADDR(zone),
                        &A1[LDA1*j], 1,
                        &W[LDW*j], 1);
            }
        }
        else {
            printf("Left Upper/NoTrans & Lower/ConjTrans not implemented yet\n");
            return MORSE_ERR_NOT_SUPPORTED;

        }
    }
    else { //side right

        if (((trans == MorseConjTrans) && (uplo == CblasUpper)) ||
            ((trans == MorseNoTrans) && (uplo == CblasLower))) {
            printf("Right Upper/ConjTrans & Lower/NoTrans not implemented yet\n");
            return MORSE_ERR_NOT_SUPPORTED;

        }
        else {

            /*
             * W = A1 + A2 * V
             */

            if (L > 0) {

                /* W = A2_2 */
                LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,
                    morse_lapack_const(MorseUpperLower),
                    M, L,
                    &A2[LDA2*(K-L)], LDA2, W, LDW);

                /* W = W * V_2 --> W = A2_2 * V_2 */
                cblas_ztrmm(
                    CblasColMajor, CblasRight, (CBLAS_UPLO)uplo,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, M, L,
                    CBLAS_SADDR(zone), &V[vi2], LDV,
                    W, LDW);

                /* W = W + A2_1 * V_1 */
                if (K > L) {
                    cblas_zgemm(
                        CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                        M, L, K-L,
                        CBLAS_SADDR(zone), A2, LDA2,
                        V, LDV,
                        CBLAS_SADDR(zone), W, LDW);
                }

            }

            /* W = W + A2 * V_3 */
            if (N > L) {
                cblas_zgemm(
                    CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                    M, N-L, K,
                    CBLAS_SADDR(zone), A2, LDA2,
                    &V[vi3], LDV,
                    CBLAS_SADDR(zzero), &W[LDW*L], LDW);
            }

            /* W = A1 + W */
            for (j = 0; j < N; j++) {
                cblas_zaxpy(
                        M, CBLAS_SADDR(zone),
                        &A1[LDA1*j], 1,
                        &W[LDW*j], 1);
            }
        }
    }

    return MORSE_SUCCESS;
}

/***************************************************************************/
static inline int
CORE_zpamm_a2(MORSE_enum side, MORSE_enum trans, MORSE_enum uplo,
              int M, int N, int K, int L,
              int vi2, int vi3,
                    MORSE_Complex64_t *A2, int LDA2,
              const MORSE_Complex64_t *V, int LDV,
                    MORSE_Complex64_t *W, int LDW)
{

   /*
    *  A2 = A2 + op(V) * W  or  A2 = A2 + W * op(V)
    */

    int j;
    static MORSE_Complex64_t zone  =  1.0;
    static MORSE_Complex64_t mzone  =  -1.0;

    if (side == MorseLeft) {

        if (((trans == MorseConjTrans) && (uplo == CblasUpper)) ||
            ((trans == MorseNoTrans) && (uplo == CblasLower))) {

            printf("Left Upper/ConjTrans & Lower/NoTrans not implemented yet\n");
            return MORSE_ERR_NOT_SUPPORTED;

        }
        else {  //trans

            /*
             * A2 = A2 - V * W
             */

            /* A2_1 = A2_1 - V_1  * W_1 */
            if (M > L) {
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    M-L, N, L,
                    CBLAS_SADDR(mzone), V, LDV,
                    W, LDW,
                    CBLAS_SADDR(zone), A2, LDA2);
            }

            /* W_1 = V_2 * W_1 */
            cblas_ztrmm(
                CblasColMajor, CblasLeft, (CBLAS_UPLO)uplo,
                (CBLAS_TRANSPOSE)trans, CblasNonUnit, L, N,
                CBLAS_SADDR(zone), &V[vi2], LDV,
                W, LDW);

            /* A2_2 = A2_2 - W_1 */
            for(j = 0; j < N; j++) {
                cblas_zaxpy(
                    L, CBLAS_SADDR(mzone),
                    &W[LDW*j], 1,
                    &A2[LDA2*j+(M-L)], 1);
            }

            /* A2 = A2 - V_3  * W_2 */
            if (K > L) {
                cblas_zgemm(
                    CblasColMajor, (CBLAS_TRANSPOSE)trans, CblasNoTrans,
                    M, N, (K-L),
                    CBLAS_SADDR(mzone), &V[vi3], LDV,
                    &W[L], LDW,
                    CBLAS_SADDR(zone), A2, LDA2);
            }

        }
    }
    else { //side right

        if (((trans == MorseConjTrans) && (uplo == CblasUpper)) ||
            ((trans == MorseNoTrans) && (uplo == CblasLower))) {

            /*
             * A2 = A2 - W * V'
             */

            /* A2 = A2 - W_2 * V_3' */
            if (K > L) {
                cblas_zgemm(
                    CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                    M, N, K-L,
                    CBLAS_SADDR(mzone), &W[LDW*L], LDW,
                    &V[vi3], LDV,
                    CBLAS_SADDR(zone), A2, LDA2);
            }

            /* A2_1 = A2_1 - W_1 * V_1' */
            if (N > L) {
                cblas_zgemm(
                    CblasColMajor, CblasNoTrans, (CBLAS_TRANSPOSE)trans,
                    M, N-L, L,
                    CBLAS_SADDR(mzone), W, LDW,
                    V, LDV,
                    CBLAS_SADDR(zone), A2, LDA2);
            }

            /* A2_2 =  A2_2 -  W_1 * V_2' */
            if (L > 0) {
                cblas_ztrmm(
                    CblasColMajor, CblasRight, (CBLAS_UPLO)uplo,
                    (CBLAS_TRANSPOSE)trans, CblasNonUnit, M, L,
                    CBLAS_SADDR(mzone), &V[vi2], LDV,
                    W, LDW);

                for (j = 0; j < L; j++) {
                    cblas_zaxpy(
                        M, CBLAS_SADDR(zone),
                        &W[LDW*j], 1,
                        &A2[LDA2*(N-L+j)], 1);
                }
            }

        }
        else {
            printf("Right Upper/NoTrans & Lower/ConjTrans not implemented yet\n");
            return MORSE_ERR_NOT_SUPPORTED;
        }
    }

    return MORSE_SUCCESS;
}



/***************************************************************************//**
 *
 **/
