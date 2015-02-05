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
 * @file codelet_zgeadd.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/quark/include/morse_quark.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgeadd adds to matrices together.
 *
 *       B <- alpha * A  + B
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/

void MORSE_TASK_zgeadd(MORSE_option_t *options,
                       int m, int n, MORSE_Complex64_t alpha,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GEADD;
    QUARK_Insert_Task(opt->quark, CORE_zgeadd_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(MORSE_Complex64_t),         &alpha, VALUE,
        sizeof(MORSE_Complex64_t)*lda*n,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(MORSE_Complex64_t)*ldb*n,    RTBLKADDR(B, MORSE_Complex64_t, Bm, Bn),             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}


void CORE_zgeadd_quark(Quark *quark)
{
    int M;
    int N;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_Complex64_t *B;
    int LDB;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    CORE_zgeadd(M, N, alpha, A, LDA, B, LDB);
    return;
}

