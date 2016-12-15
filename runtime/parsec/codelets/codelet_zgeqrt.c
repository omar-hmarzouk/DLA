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
 * @file codelet_zgeqrt.c
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
#include "runtime/parsec/include/morse_parsec.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgeqrt computes a QR factorization of a complex M-by-N tile A:
 *  A = Q * R.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(1) H(2) . . . H(k), where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
 *  and tau in TAU(i).
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the tile A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the tile A.  N >= 0.
 *
 * @param[in] IB
 *         The inner-blocking size.  IB >= 0.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, the elements on and above the diagonal of the array
 *         contain the min(M,N)-by-N upper trapezoidal tile R (R is
 *         upper triangular if M >= N); the elements below the diagonal,
 *         with the array TAU, represent the unitary tile Q as a
 *         product of elementary reflectors (see Further Details).
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 * @param[out] T
 *         The IB-by-N triangular factor T of the block reflector.
 *         T is upper triangular by block (economic storage);
 *         The rest of the array is not referenced.
 *
 * @param[in] LDT
 *         The leading dimension of the array T. LDT >= IB.
 *
 * @param[out] TAU
 *         The scalar factors of the elementary reflectors (see Further
 *         Details).
 *
 * @param[out] WORK
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
static int
CORE_zgeqrt_parsec (dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *m;
    int *n;
    int *ib;
    MORSE_Complex64_t *A;
    int *lda;
    MORSE_Complex64_t *T;
    int *ldt;
    MORSE_Complex64_t *TAU;
    MORSE_Complex64_t *WORK;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &m,
        UNPACK_VALUE, &n,
        UNPACK_VALUE, &ib,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &lda,
        UNPACK_DATA,  &T,
        UNPACK_VALUE, &ldt,
        UNPACK_SCRATCH, &TAU,
        UNPACK_SCRATCH, &WORK );

    CORE_zgeqrt(*m, *n, *ib, A, *lda, T, *ldt, TAU, WORK);

    return 0;
}

void MORSE_TASK_zgeqrt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zgeqrt_parsec, "geqrt",
        sizeof(int),           &m,                             VALUE,
        sizeof(int),           &n,                             VALUE,
        sizeof(int),           &ib,                            VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT | REGION_FULL,
        sizeof(int),           &lda,                           VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ),     OUTPUT | REGION_FULL,
        sizeof(int),           &ldt,                           VALUE,
        sizeof(MORSE_Complex64_t)*nb,       NULL,              SCRATCH,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,              SCRATCH,
        0);
}
