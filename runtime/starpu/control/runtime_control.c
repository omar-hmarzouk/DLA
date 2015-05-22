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
 * @file runtime_control.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "runtime/starpu/include/morse_starpu.h"

/*******************************************************************************
 * Thread rank.
 **/
int RUNTIME_rank(MORSE_context_t *morse)
{
    (void)morse;
    return starpu_worker_get_id();
}

/*******************************************************************************
 *
 **/
int RUNTIME_init_scheduler( MORSE_context_t *morse, int ncpus, int ncudas, int nthreads_per_worker)
{
    starpu_conf_t *conf = (starpu_conf_t*)(morse->schedopt);
    int hres = -1;

    conf->ncpus = ncpus;
    conf->ncuda = ncudas;
    conf->nopencl = 0;

    /* By default, use the dmdas strategy */
    if (!getenv("STARPU_SCHED"))
        conf->sched_policy_name = "dmdas";

    /* By default, enable calibration */
    if (!getenv("STARPU_CALIBRATE"))
        conf->calibrate = 1;

    /* Set scheduling to "ws" if no cuda devices used because it behaves better
     * on homognenous architecture. If the user wants to use another
     * scheduling strategy, he can set STARPU_SCHED env. var. to whatever he
     * wants
     * TODO: discuss this with ManuMathieu
     * reference, we should use "lws" strategy
     */
    if ( !getenv("STARPU_SCHED") && conf->ncuda == 0)
#if STARPU_MAJOR_VERSION > 1 || (STARPU_MAJOR_VERSION == 1 && STARPU_MINOR_VERSION >= 2)
        conf->sched_policy_name = "lws";
#else
        conf->sched_policy_name = "ws";
#endif

    if ((ncpus == -1)||(nthreads_per_worker == -1))
    {
        morse->parallel_enabled = MORSE_FALSE;

        hres = starpu_init( conf );
    }
    else {
        int worker;

        morse->parallel_enabled = MORSE_TRUE;

        for (worker = 0; worker < ncpus; worker++)
            conf->workers_bindid[worker] = (worker+1)*nthreads_per_worker - 1;

        for (worker = 0; worker < ncpus; worker++)
            conf->workers_bindid[worker + ncudas] = worker*nthreads_per_worker;

        conf->use_explicit_workers_bindid = 1;

        hres = starpu_init( conf );

        morse->nworkers = ncpus;
        morse->nthreads_per_worker = nthreads_per_worker;
    }

#if defined(CHAMELEON_USE_MPI)
    {
        int flag = 0;
        MPI_Initialized( &flag );
        starpu_mpi_init(NULL, NULL, !flag);
        MPI_Comm_rank(MPI_COMM_WORLD, &(morse->my_mpi_rank));
        MPI_Comm_size(MPI_COMM_WORLD, &(morse->mpi_comm_size));
    }
#endif

#if defined(HAVE_STARPU_FXT_PROFILING)
    starpu_fxt_stop_profiling();
#endif

#if defined(CHAMELEON_USE_CUDA)
    starpu_cublas_init();
#endif

    return hres;
}

/*******************************************************************************
 *
 */
void RUNTIME_finalize_scheduler( MORSE_context_t *morse )
{
    (void)morse;
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_shutdown();
#endif
#if defined(CHAMELEON_USE_CUDA)
    starpu_cublas_shutdown();
#endif

    starpu_shutdown();
    return;
}

/*******************************************************************************
 *  Busy-waiting barrier
 **/
void RUNTIME_barrier( MORSE_context_t *morse )
{
    (void)morse;
    starpu_task_wait_for_all();
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif
}

/*******************************************************************************
 *  To suspend the processing of new tasks by workers
 **/
void RUNTIME_pause( MORSE_context_t *morse )
{
    (void)morse;
    starpu_pause();
    return;
}

/*******************************************************************************
 *  This is the symmetrical call to RUNTIME_pause,
 *  used to resume the workers polling for new tasks.
 **/
void RUNTIME_resume( MORSE_context_t *morse )
{
    (void)morse;
    starpu_resume();
    return;
}
