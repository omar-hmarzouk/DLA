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
 * @file core_zgesplit.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
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
#include "coreblas/include/lapacke.h"
#include "coreblas/include/coreblas.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgesplit splits a matrix A into two parts (Upper/Lower),
 *  A keeps its lower/upper part unchanged and the other part is filled with zeros .
 *  ones can be optionally set on the diagonal.
 *  The part of A which is erased is copied in B.
 *
 *
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in,out] A
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

int CORE_zgesplit(MORSE_enum side, MORSE_enum diag,
                  int M, int N,
                  MORSE_Complex64_t *A, int LDA,
                  MORSE_Complex64_t *B, int LDB)
{

	MORSE_enum uplo;
	int i, j;

    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        coreblas_error(7, "Illegal value of LDB");
        return -7;
    }

    if (side == MorseLeft){
    	uplo = MorseUpper;
    } else{
    	uplo = MorseLower;
    }

    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR,
                        morse_lapack_const(uplo),
                        M, N, A, LDA, B, LDB);

    LAPACKE_zlaset_work(LAPACK_COL_MAJOR,
                        morse_lapack_const(uplo),
                        M, N, 0., 1., A, LDA);

    return MORSE_SUCCESS;
}



