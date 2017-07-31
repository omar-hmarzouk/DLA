/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 * @file codelet_zlaset.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 **/
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zlaset - Sets the elements of the matrix A on the diagonal
 *  to beta and on the off-diagonals to alpha
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = MorseUpper: Upper part of A is set;
 *          = MorseLower: Lower part of A is set;
 *          = MorseUpperLower: ALL elements of A are set.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 **/
static int
CORE_zlaset_parsec(dague_execution_unit_t *context, dague_execution_context_t * this_task)
{
    MORSE_enum *uplo;
    int *M;
    int *N;
    MORSE_Complex64_t *alpha;
    MORSE_Complex64_t *beta;
    MORSE_Complex64_t *A;
    int *LDA;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &uplo,
        UNPACK_VALUE, &M,
        UNPACK_VALUE, &N,
        UNPACK_VALUE, &alpha,
        UNPACK_VALUE, &beta,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &LDA );

    CORE_zlaset(*uplo, *M, *N, *alpha, *beta, A, *LDA);

    return 0;
}

void MORSE_TASK_zlaset(const MORSE_option_t *options,
                       MORSE_enum uplo, int M, int N,
                       MORSE_Complex64_t alpha, MORSE_Complex64_t beta,
                       const MORSE_desc_t *A, int Am, int An, int LDA)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zlaset_parsec, "laset",
        sizeof(MORSE_enum),              &uplo,        VALUE,
        sizeof(int),                     &M,           VALUE,
        sizeof(int),                     &N,           VALUE,
        sizeof(MORSE_Complex64_t),       &alpha,       VALUE,
        sizeof(MORSE_Complex64_t),       &beta,        VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     OUTPUT | REGION_FULL,
        sizeof(int),                     &LDA,         VALUE,
        0);
}
