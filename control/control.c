/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University 
 *                          of Tennessee Research Foundation. 
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 IPB. All rights reserved. 
 *
 **/

/**
 *
 * @file control.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include "auxiliary.h"
#include "common.h"
#include "runtime.h"

/*******************************************************************************
 *
 * @ingroup Auxiliary
 *
 *  MORSE_Init - Initialize MORSE.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = MORSE_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int MORSE_Init(int cores, int gpus)
{
    return MORSE_InitPar(cores, gpus, -1);
}

/*******************************************************************************
 *
 * @ingroup Auxiliary
 *
 *  MORSE_Init_Affinity - Initialize MORSE.
 *
 *******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use (threads to launch).
 *          If cores = 0, cores = MORSE_NUM_THREADS if it is set, the
 *          system number of core otherwise.
 *
 * @param[in] coresbind
 *          Array to specify where to bind each thread.
 *          Each thread i is binded to coresbind[hwloc(i)] if hwloc is
 *          provided, or to coresbind[i] otherwise.
 *          If coresbind = NULL, coresbind = MORSE_AFF_THREADS if it
 *          is set, the identity function otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int MORSE_InitPar(int nworkers, int ncudas, int nthreads_per_worker)
{
    MORSE_context_t *morse;

    /* Create context and insert in the context map */
    morse = morse_context_create();
    if (morse == NULL) {
        morse_fatal_error("MORSE_Init", "morse_context_create() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }

#if 0    
    /* Init number of cores and topology */
    morse_topology_init();

    /* Set number of nworkers */
    if ( nworkers < 1 ) {
        morse->world_size = morse_get_numthreads();
        if ( morse->world_size == -1 ) {
            morse->world_size = 1;
            morse_warning("MORSE_Init", "Could not find the number of cores: the thread number is set to 1");
        }
    }
    else
      morse->world_size = nworkers;

    if (morse->world_size <= 0) {
        morse_fatal_error("MORSE_Init", "failed to get system size");
        return MORSE_ERR_NOT_FOUND;
    }
    nworkers = morse->world_size;
    
    /* Get the size of each NUMA node */
    morse->group_size = morse_get_numthreads_numa();
    while ( ((morse->world_size)%(morse->group_size)) != 0 ) 
        (morse->group_size)--;
#endif

#if defined(CHAMELEON_USE_MPI)
    {
      int flag = 0, provided = 0;
      MPI_Initialized( &flag );
      if ( !flag ) {
	MPI_Init_thread( NULL, NULL, MPI_THREAD_MULTIPLE, &provided );
      }
    }
#endif
#if defined(CHAMELEON_USE_MAGMA)
    magma_init();
#endif
    RUNTIME_init_scheduler( morse, nworkers, ncudas, nthreads_per_worker );
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *
 * @ingroup Auxiliary
 *
 *  MORSE_Finalize - Finalize MORSE.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int MORSE_Finalize(void)
{
    MORSE_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Finalize()", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    MORSE_TASK_dataflush_all();
    RUNTIME_finalize_scheduler( morse );
#if defined(CHAMELEON_USE_MAGMA)
    magma_finalize();
#endif
    morse_context_destroy();
#if defined(CHAMELEON_USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return MORSE_SUCCESS;
}

/*******************************************************************************
 *
 * @ingroup Auxiliary
 *
 *  MORSE_My_Mpi_Rank - Return the MPI rank of the calling process.
 *
 *******************************************************************************
 *
 *******************************************************************************
 *
 * @return
 *          \retval MPI rank
 *
 ******************************************************************************/
int MORSE_My_Mpi_Rank(void)
{
#if defined(CHAMELEON_USE_MPI)
    MORSE_context_t *morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Finalize()", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    return MORSE_MPI_RANK;
#else
    return 0;
#endif
}
