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
 * @file runtime_descriptor.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "runtime/starpu/include/morse_starpu.h"

#if defined(CHAMELEON_USE_MPI)

/* Take 24 bits for the tile id, and 7 bits for descriptor id */
#define TAG_WIDTH_MIN 20
static int tag_width = 31;
static int tag_sep   = 24;

#ifndef HAVE_STARPU_MPI_DATA_REGISTER
#define starpu_mpi_data_register( handle_, tag_, owner_ )        \
    do {                                                         \
     starpu_data_set_rank( (handle_), (owner_) );                \
     starpu_data_set_tag( (handle_), (tag_) );                   \
    } while(0)
#endif

#endif

void RUNTIME_desc_create( MORSE_desc_t *desc )
{
    int64_t lmt = desc->lmt;
    int64_t lnt = desc->lnt;
    int64_t block_ind = 0;
    starpu_data_handle_t *tiles;

    desc->occurences = 1;

    /*
     * Allocate starpu_handle_t array (handlers are initialized on the fly when
     * discovered by any algorithm to save space)
     */
    desc->schedopt = (void*)calloc(lnt*lmt,sizeof(starpu_data_handle_t));
    assert(desc->schedopt);
    tiles = (starpu_data_handle_t*)(desc->schedopt);

#if defined(CHAMELEON_USE_CUDA)
    if (desc->use_mat == 1 && desc->register_mat == 1){
        /*
         * Register allocated memory as CUDA pinned memory
         */
        {
            int64_t eltsze = MORSE_Element_Size(desc->dtyp);
            size_t size = (size_t)(desc->llm) * (size_t)(desc->lln) * eltsze;

            /* Register the matrix as pinned memory */
            if ( cudaHostRegister( desc->mat, size, cudaHostRegisterPortable ) != cudaSuccess )
            {
                morse_warning("RUNTIME_desc_create(StarPU)", "cudaHostRegister failed to register the matrix as pinned memory");
            }
        }
    }
#endif

#if defined(CHAMELEON_USE_MPI)
    /*
     * Check that we are not going over MPI tag limitations
     */
    {
        static int _tag_mpi_initialized_ = 0;
        MORSE_context_t *morse;
        int myrank = desc->myrank;

        if (!_tag_mpi_initialized_) {
            int *tag_ub = NULL;
            int ok = 0;

            MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &ok);

            if ( !ok ) {
                morse_error("RUNTIME_desc_create", "MPI_TAG_UB not known by MPI");
            }

            while ( ((uintptr_t)(1UL<<tag_width) >= (uint)(*tag_ub) )
                    && (tag_width >= TAG_WIDTH_MIN) ) {
                tag_width--;
                tag_sep--;
            }

            if ( tag_width < TAG_WIDTH_MIN ) {
                morse_error("RUNTIME_desc_create", "MPI_TAG_UB too small to identify all the data");
                return;
            }

            _tag_mpi_initialized_ = 1;
        }

        /* Check that we won't create overflow in tags used */
        if ( (lnt*lmt) > ((uintptr_t)(1UL<<tag_sep)) ) {
            morse_error("RUNTIME_desc_create", "Too many tiles in the descriptor for MPI tags");
            return;
        }
        assert(lmt*lmt<=(1<<tag_sep));

        if (desc->id >= 1UL<<(tag_width-tag_sep)) {
            morse_error("RUNTIME_desc_create", "Number of descriptor available in MPI mode out of stock");
            return;
        }
        assert( desc->id < (1UL<<(tag_width-tag_sep)) );
    }
#endif
}

void RUNTIME_desc_destroy( MORSE_desc_t *desc )
{
    desc->occurences--;

    /*
     * If this is the last descriptor using the matrix, we release the handle
     * and unregister the GPU data
     */
    if ( desc->occurences == 0 ) {
        starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
        int lmt = desc->lmt;
        int lnt = desc->lnt;
        int m, n;

        for (n = 0; n < lnt; n++)
            for (m = 0; m < lmt; m++)
            {
                if (*handle == NULL)
                {
                    handle++;
                    continue;
                }
                //printf("\nUnregister %d %d %d", MORSE_My_Mpi_Rank(), m, n);
                starpu_data_unregister(*handle);
                handle++;
            }

#if defined(CHAMELEON_USE_CUDA)
        if (desc->use_mat == 1 && desc->register_mat == 1){
            int64_t eltsze = MORSE_Element_Size(desc->dtyp);
            size_t size = (size_t)(desc->llm) * (size_t)(desc->lln) * eltsze;

            /* Unmap the pinned memory associated to the matrix */
            if (cudaHostUnregister(desc->mat) != cudaSuccess)
            {
                morse_warning("RUNTIME_desc_destroy(StarPU)",
                              "cudaHostUnregister failed to unregister the "
                              "pinned memory associated to the matrix");
            }
        }
#endif /* defined(CHAMELEON_USE_CUDA) */

        free(desc->schedopt);
    }
}

void RUNTIME_desc_init( MORSE_desc_t *desc )
{
  (void)desc;
  return;
}

void RUNTIME_desc_submatrix( MORSE_desc_t *desc )
{
    desc->occurences++;
    return;
}

/* TODO: Acquire/Release/GetonCPU need to be studied carefully and fixed
 * because we are not using them correctly */
int RUNTIME_desc_acquire( MORSE_desc_t *desc )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++)
        for (m = 0; m < lmt; m++)
        {
            if ( (*handle == NULL) ||
                 !morse_desc_islocal( desc, m, n ) )
            {
                handle++;
                continue;
            }
            starpu_data_acquire(*handle, STARPU_R);
            handle++;
        }
    return MORSE_SUCCESS;
}

int RUNTIME_desc_release( MORSE_desc_t *desc )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++)
        for (m = 0; m < lmt; m++)
        {
            if ( (*handle == NULL) ||
                 !morse_desc_islocal( desc, m, n ) )
            {
                handle++;
                continue;
            }
            starpu_data_release(*handle);
            handle++;
        }
    return MORSE_SUCCESS;
}

int RUNTIME_desc_getoncpu( MORSE_desc_t *desc )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++)
        for (m = 0; m < lmt; m++)
        {
            if ( (*handle == NULL) ||
                 !morse_desc_islocal( desc, m, n ) )
            {
                handle++;
                continue;
            }

            starpu_data_acquire(*handle, STARPU_R);
            starpu_data_release(*handle);
            handle++;
        }

    return MORSE_SUCCESS;
}

void *RUNTIME_desc_getaddr( MORSE_desc_t *desc, int m, int n )
{
    starpu_data_handle_t *ptrtile = (starpu_data_handle_t*)(desc->schedopt);
    ptrtile += ((int64_t)(desc->lmt) * (int64_t)n + (int64_t)m);

    if (*ptrtile == NULL) {
        int64_t block_ind = desc->lmt * n + m;
        int64_t eltsze = MORSE_Element_Size(desc->dtyp);
        int myrank = desc->myrank;
        int owner  = desc->get_rankof( desc, m, n );
        int tempmm = (m == desc->lmt-1) ? (desc->lm - m * desc->mb) : desc->mb;
        int tempnn = (n == desc->lnt-1) ? (desc->ln - n * desc->nb) : desc->nb;

        if ( myrank == owner ) {
            //printf("\nRegister %d %d %d", MORSE_My_Mpi_Rank(), m, n);
            starpu_matrix_data_register(ptrtile, 0,
                                        (uintptr_t)desc->get_blkaddr(desc, m, n),
                                        BLKLDD(desc, m), tempmm, tempnn, eltsze);
        }
        else {
            //printf("\nRegister dist %d %d %d", MORSE_My_Mpi_Rank(), m, n);
            starpu_matrix_data_register(ptrtile, -1,
                                        (uintptr_t) NULL,
                                        BLKLDD(desc, m), tempmm, tempnn, eltsze);
        }

#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_data_register(*ptrtile, (desc->id << tag_sep) | (block_ind), owner);
#endif /* defined(CHAMELEON_USE_MPI) */
    }

    return (void *)(*ptrtile);
}
