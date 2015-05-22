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
 * @file codelet_ztrasm.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Inria Bordeaux - Sud-Ouest, LaBRI,
 *  University of Bordeaux, Bordeaux INP
 *
 * @version 2.6.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/starpu/include/morse_starpu.h"
#include "runtime/starpu/include/runtime_codelet_z.h"

void MORSE_TASK_ztrasm(MORSE_option_t *options,
                       MORSE_enum storev, MORSE_enum uplo, MORSE_enum diag, int M, int N,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn)
{
    struct starpu_codelet *codelet = &cl_ztrasm;
    void (*callback)(void*) = options->profiling ? cl_ztrasm_callback : NULL;
    starpu_insert_task(
        codelet,
        STARPU_VALUE,    &storev,                    sizeof(MORSE_enum),
        STARPU_VALUE,    &uplo,                      sizeof(MORSE_enum),
        STARPU_VALUE,    &diag,                      sizeof(MORSE_enum),
        STARPU_VALUE,    &M,                         sizeof(int),
        STARPU_VALUE,    &N,                         sizeof(int),
        STARPU_R,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                       sizeof(int),
        STARPU_RW,       RTBLKADDR(B, double, Bm, Bn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
        0);
}


static void cl_ztrasm_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum storev;
    MORSE_enum uplo;
    MORSE_enum diag;
    int M;
    int N;
    MORSE_Complex64_t *A;
    int lda;
    double *work;

    A     = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    work  = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &storev, &uplo, &diag, &M, &N, &lda);
    CORE_ztrasm(storev, uplo, diag, M, N, A, lda, work);
}

/*
 * Codelet definition
 */
CODELETS_CPU(ztrasm, 2, cl_ztrasm_cpu_func)
