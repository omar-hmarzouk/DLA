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
 * @file codelet_zgelqt.c
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
#include "coreblas/coreblas_z.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zgelqt - computes a LQ factorization of a complex M-by-N tile A: A = L * Q.
 *
 *  The tile Q is represented as a product of elementary reflectors
 *
 *    Q = H(k)' . . . H(2)' H(1)', where k = min(M,N).
 *
 *  Each H(i) has the form
 *
 *    H(i) = I - tau * v * v'
 *
 *  where tau is a complex scalar, and v is a complex vector with
 *  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
 *  A(i,i+1:n), and tau in TAU(i).
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
 *         On exit, the elements on and below the diagonal of the array
 *         contain the M-by-min(M,N) lower trapezoidal tile L (L is
 *         lower triangular if M <= N); the elements above the diagonal,
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
CORE_zgelqt_parsec(parsec_execution_stream_t *context, parsec_task_t *this_task)
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

    parsec_dtd_unpack_args(
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

    CORE_zgelqt(*m, *n, *ib, A, *lda, T, *ldt, TAU, WORK);

    return 0;
}

void MORSE_TASK_zgelqt(const MORSE_option_t *options,
                       int m, int n, int ib, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *T, int Tm, int Tn, int ldt)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zgelqt_parsec, options->priority, "gelqt",
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(int),                        &ib,    VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INOUT,
        sizeof(int),           &lda,                VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ),     OUTPUT,
        sizeof(int),           &ldt,                VALUE,
        sizeof(MORSE_Complex64_t)*nb,       NULL,   SCRATCH,
        sizeof(MORSE_Complex64_t)*ib*nb,    NULL,   SCRATCH,
        0);
}
