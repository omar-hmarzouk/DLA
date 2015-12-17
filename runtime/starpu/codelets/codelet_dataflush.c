/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file codelet_dataflush.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Marc Sergent
 * @date 2014-02-05
 *
 **/
#include "runtime/starpu/include/morse_starpu.h"

#ifdef CHAMELEON_USE_STARPU_DATA_WONT_USE
#elif defined CHAMELEON_USE_STARPU_IDLE_PREFETCH
static void data_flush(void *handle)
{
        starpu_data_idle_prefetch_on_node(handle, STARPU_MAIN_RAM, 1);
        starpu_data_release_on_node(handle, -1);
}
#else
static void data_release(void *handle)
{
        starpu_data_release(handle);
}
#endif

void MORSE_TASK_dataflush(MORSE_option_t *options,
                          MORSE_desc_t *A, int Am, int An)
{
    (void)options;

    /*
     * We can use MORSE_Complex64_t for all precisions since it is not use to
     * compute the handle address in starpu.  We have to be careful with this if
     * something similar happen in Quark.
     */
    {
        starpu_data_handle_t *ptrtile = (starpu_data_handle_t*)(A->schedopt);
        ptrtile += ((int64_t)(A->lmt) * (int64_t)An + (int64_t)Am);

        if (*ptrtile != NULL)
        {
#if defined(CHAMELEON_USE_MPI)
            starpu_mpi_cache_flush(MPI_COMM_WORLD, *ptrtile);
#endif

            if ( A->myrank == A->get_rankof( A, Am, An ) )
            {
                /* Push data to main memory when we have time to */
#ifdef CHAMELEON_USE_STARPU_DATA_WONT_USE
                starpu_data_wont_use(*ptrtile);
#elif defined CHAMELEON_USE_STARPU_IDLE_PREFETCH
                starpu_data_acquire_on_node_cb(*ptrtile, -1, STARPU_R, data_flush, *ptrtile);
#else
                starpu_data_acquire_cb(*ptrtile, STARPU_R, data_release, *ptrtile);
#endif
            }
        }
    }
}

void MORSE_TASK_dataflush_all()
{
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_cache_flush_all_data(MPI_COMM_WORLD);
#endif
}
