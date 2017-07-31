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
 * @file codelet_zgessq.c
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
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"
#include "chameleon/morse_tasks_z.h"

void MORSE_TASK_zgessq( const MORSE_option_t *options,
                        int m, int n,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        const MORSE_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn )
{
    struct starpu_codelet *codelet = &cl_zgessq;
    void (*callback)(void*) = options->profiling ? cl_zgessq_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_RW(SCALESUMSQ, SCALESUMSQm, SCALESUMSQn);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &m,                          sizeof(int),
        STARPU_VALUE,    &n,                          sizeof(int),
        STARPU_R,        RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                        sizeof(int),
        STARPU_RW,       RTBLKADDR(SCALESUMSQ, double, SCALESUMSQm, SCALESUMSQn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zgessq",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zgessq_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    double *SCALESUMSQ;

    A          = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    SCALESUMSQ = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda);
    CORE_zgessq( m, n, A, lda, &SCALESUMSQ[0], &SCALESUMSQ[1] );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgessq, 2, cl_zgessq_cpu_func)
