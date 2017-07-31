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
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"

static void
CORE_ztpqrt_quark( Quark *quark )
{
    int M;
    int N;
    int L;
    int ib;
    MORSE_Complex64_t *A;
    int lda;
    MORSE_Complex64_t *B;
    int ldb;
    MORSE_Complex64_t *T;
    int ldt;
    MORSE_Complex64_t *WORK;

    quark_unpack_args_11( quark, M, N, L, ib,
                          A, lda, B, ldb, T, ldt, WORK );

    CORE_ztpqrt( M, N, L, ib,
                 A, lda, B, ldb, T, ldt, WORK );
}

void MORSE_TASK_ztpqrt( const MORSE_option_t *options,
                         int M, int N, int L, int ib, int nb,
                         const MORSE_desc_t *A, int Am, int An, int lda,
                         const MORSE_desc_t *B, int Bm, int Bn, int ldb,
                         const MORSE_desc_t *T, int Tm, int Tn, int ldt )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_TSQRT;

    QUARK_Insert_Task(
        opt->quark, CORE_ztpqrt_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                         &M,   VALUE,
        sizeof(int),                         &N,   VALUE,
        sizeof(int),                         &L,   VALUE,
        sizeof(int),                         &ib,  VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,      RTBLKADDR( A, MORSE_Complex64_t, Am, An ), INOUT | QUARK_REGION_U | QUARK_REGION_D,
        sizeof(int),                         &lda, VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,      RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), INOUT,
        sizeof(int),                         &ldb, VALUE,
        sizeof(MORSE_Complex64_t)*nb*ib,      RTBLKADDR( T, MORSE_Complex64_t, Tm, Tn ), OUTPUT,
        sizeof(int),                         &ldt, VALUE,
        sizeof(MORSE_Complex64_t)*(ib+1)*nb,  NULL, SCRATCH,
        0);
}
