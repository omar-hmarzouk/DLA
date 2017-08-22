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
 * @file codelet_zaxpy.c
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

static int
CORE_zaxpy_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *M;
    MORSE_Complex64_t *alpha;
    MORSE_Complex64_t *A;
    int *incA;
    MORSE_Complex64_t *B;
    int *incB;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &M,
        UNPACK_VALUE, &alpha,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &incA,
        UNPACK_DATA,  &B,
        UNPACK_VALUE, &incB );

    CORE_zaxpy(*M, *alpha, A, *incA, B, *incB);

    return 0;
}

void MORSE_TASK_zaxpy(const MORSE_option_t *options,
                      int M, MORSE_Complex64_t *alpha,
                      const MORSE_desc_t *A, int Am, int An, int incA,
                      const MORSE_desc_t *B, int Bm, int Bn, int incB)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zaxpy_parsec, "axpy",
        sizeof(int),               &M,     VALUE,
        sizeof(MORSE_Complex64_t), &alpha, VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, MORSE_Complex64_t, Am, An ), INPUT | REGION_FULL,
        sizeof(int),               &incA, VALUE,
        PASSED_BY_REF,  RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), INOUT | REGION_FULL,
        sizeof(int),               &incB, VALUE,
        0);
}
