/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
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
#include "chameleon_starpu.h"
#include "chameleon/morse_tasks.h"

#ifdef HAVE_STARPU_DATA_WONT_USE
#elif defined HAVE_STARPU_IDLE_PREFETCH
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

void MORSE_TASK_flush_data( const MORSE_option_t *options,
                            const MORSE_desc_t *A, int Am, int An )
{
    (void)options;

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
#ifdef HAVE_STARPU_DATA_WONT_USE
                starpu_data_wont_use(*ptrtile);
#elif defined HAVE_STARPU_IDLE_PREFETCH
                starpu_data_acquire_on_node_cb(*ptrtile, -1, STARPU_R, data_flush, *ptrtile);
#else
                starpu_data_acquire_cb(*ptrtile, STARPU_R, data_release, *ptrtile);
#endif
            }
        }
    }
}

void MORSE_TASK_flush_desc( const MORSE_option_t *options,
                            MORSE_enum uplo, const MORSE_desc_t *A )
{
    int m, n;

    switch (uplo) {
    /*
     *  MorseUpper
     */
    case MorseUpper:
        for (m = 0; m < A->mt; m++) {
            for (n = m; n < A->nt; n++) {
                MORSE_TASK_flush_data( options, A, m, n );
            }
        }
        break;
    /*
     *  MorseLower
     */
    case MorseLower:
        for (m = 0; m < A->mt; m++) {
            for (n = 0; n < chameleon_min(m+1, A->nt); n++) {
                MORSE_TASK_flush_data( options, A, m, n );
            }
        }
        break;
    /*
     *  MorseUpperLower
     */
    case MorseUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            for (n = 0; n < A->nt; n++) {
                MORSE_TASK_flush_data( options, A, m, n );
            }
        }
    }
}

void MORSE_TASK_flush_all()
{
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_cache_flush_all_data(MPI_COMM_WORLD);
#endif
}
