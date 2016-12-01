/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014, 2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_ztrtri.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

#include "runtime/starpu/include/morse_starpu.h"
#include "runtime/starpu/include/runtime_codelet_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_ztrtri(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum diag,
                       int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_ztrtri;
    void (*callback)(void*) = options->profiling ? cl_ztrtri_callback : NULL;

    if ( morse_desc_islocal( A, Am, An ) )
    {
        starpu_insert_task(
            codelet,
            STARPU_VALUE,    &uplo,              sizeof(MORSE_enum),
            STARPU_VALUE,    &diag,              sizeof(MORSE_enum),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_RW,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_VALUE,    &iinfo,             sizeof(int),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
            0);
    }
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrtri_cpu_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    MORSE_enum diag;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    int iinfo;

    int info = 0;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &N, &LDA, &iinfo);
    CORE_ztrtri(uplo, diag, N, A, LDA, &info);
}

#if defined(CHAMELEON_USE_MAGMA)
static void cl_ztrtri_cuda_func(void *descr[], void *cl_arg)
{
    MORSE_enum uplo;
    MORSE_enum diag;
    int N;
    cuDoubleComplex *A;
    int LDA;
    int iinfo;
    int info = 0;

    A = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &diag, &N, &LDA, &iinfo);
    CUDA_ztrtri( uplo, diag, N, A, LDA, &info);
    cudaThreadSynchronize();
    return;
}
#endif
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_MAGMA)
CODELETS(ztrtri, 1, cl_ztrtri_cpu_func, cl_ztrtri_cuda_func, 0)
#else
CODELETS_CPU(ztrtri, 1, cl_ztrtri_cpu_func)
#endif
