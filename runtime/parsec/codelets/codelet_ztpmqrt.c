/**
 *
 * @copyright (c) 2009-2016 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                          Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_ztpqrt.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @date 2016-12-15
 * @precisions normal z -> s d c
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

static int
CORE_ztpmqrt_parsec(dague_execution_unit_t    *context,
                    dague_execution_context_t *this_task)
{
    MORSE_enum *side;
    MORSE_enum *trans;
    int *M;
    int *N;
    int *K;
    int *L;
    int *ib;
    const MORSE_Complex64_t *V;
    int *ldv;
    const MORSE_Complex64_t *T;
    int *ldt;
    MORSE_Complex64_t *A;
    int *lda;
    MORSE_Complex64_t *B;
    int *ldb;
    MORSE_Complex64_t *WORK;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE,   &side,
        UNPACK_VALUE,   &trans,
        UNPACK_VALUE,   &M,
        UNPACK_VALUE,   &N,
        UNPACK_VALUE,   &K,
        UNPACK_VALUE,   &L,
        UNPACK_VALUE,   &ib,
        UNPACK_DATA,    &V,
        UNPACK_VALUE,   &ldv,
        UNPACK_DATA,    &T,
        UNPACK_VALUE,   &ldt,
        UNPACK_DATA,    &A,
        UNPACK_VALUE,   &lda,
        UNPACK_DATA,    &B,
        UNPACK_VALUE,   &ldb,
        UNPACK_SCRATCH, &WORK );

    CORE_ztpmqrt( *side, *trans, *M, *N, *K, *L, *ib,
                  V, *ldv, T, *ldt, A, *lda, B, *ldb, WORK );

    return 0;
}

void MORSE_TASK_ztpmqrt( const MORSE_option_t *options,
                         MORSE_enum side, MORSE_enum trans,
                         int M, int N, int K, int L, int ib, int nb,
                         const MORSE_desc_t *V, int Vm, int Vn, int ldv,
                         const MORSE_desc_t *T, int Tm, int Tn, int ldt,
                         const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb )
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_ztpmqrt_parsec, "tpmqrt",
        sizeof(MORSE_enum), &side,  VALUE,
        sizeof(MORSE_enum), &trans, VALUE,
        sizeof(int),        &M,     VALUE,
        sizeof(int),        &N,     VALUE,
        sizeof(int),        &K,     VALUE,
        sizeof(int),        &L,     VALUE,
        sizeof(int),        &ib,    VALUE,
        PASSED_BY_REF,       RTBLKADDR( V, MORSE_Complex64_t, Vm, Vn ), INPUT | REGION_FULL,
        sizeof(int),        &ldv,   VALUE,
        PASSED_BY_REF,       RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ), INPUT | REGION_FULL,
        sizeof(int),        &ldt,   VALUE,
        PASSED_BY_REF,       RTBLKADDR( A, MORSE_Complex64_t, Am, An ), INOUT | REGION_FULL,
        sizeof(int),        &lda,   VALUE,
        PASSED_BY_REF,       RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), INOUT | REGION_FULL,
        sizeof(int),        &ldb,   VALUE,
        sizeof(MORSE_Complex64_t)*ib*nb, NULL, SCRATCH,
        0);
}
