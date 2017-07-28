/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University of
 *                          Tennessee Research Foundation.  All rights reserved.
 * @copyright (c) 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                          Univ. Bordeaux. All rights reserved.
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
#include <unistd.h>
#include "runtime/starpu/include/morse_starpu.h"

#if defined(CHAMELEON_USE_MPI)

/* Take 24 bits for the tile id, and 7 bits for descriptor id.
 These values can be changed through the call MORSE_user_tag_size(int tag_width, int tag_sep) */
#define TAG_WIDTH_MIN 20
static int tag_width = 31;
static int tag_sep   = 24;
static int _tag_mpi_initialized_ = 0;

#ifndef HAVE_STARPU_MPI_DATA_REGISTER
#define starpu_mpi_data_register( handle_, tag_, owner_ )       \
    do {                                                        \
        starpu_data_set_rank( (handle_), (owner_) );            \
        starpu_data_set_tag( (handle_), (tag_) );               \
    } while(0)
#endif

#endif

#ifdef STARPU_MALLOC_SIMULATION_FOLDED
#define FOLDED STARPU_MALLOC_SIMULATION_FOLDED
#else
#define FOLDED 0
#endif

void RUNTIME_user_tag_size( int user_tag_width, int user_tag_sep ) {
#if defined(CHAMELEON_USE_MPI)
    if (_tag_mpi_initialized_ == 0) {
        tag_width = user_tag_width;
        tag_sep   = user_tag_sep;
    } else {
        morse_error("RUNTIME_user_tag_size",
                    "must be called before creating any Morse descriptor with MORSE_Desc_create(). The tag sizes will not be modified.");
    }
#endif
    (void)user_tag_width; (void)user_tag_sep;
}


void *RUNTIME_mat_alloc( size_t size )
{
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(CHAMELEON_USE_MPI)
    return (void*) 1;
#else
    void *mat;

    if (starpu_malloc_flags(&mat, size, STARPU_MALLOC_PINNED|FOLDED|STARPU_MALLOC_COUNT) != 0)
        return NULL;
    return mat;
#endif
}

void RUNTIME_mat_free( void *mat, size_t size )
{
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(CHAMELEON_USE_MPI)
    return;
#else
    starpu_free_flags(mat, size, STARPU_MALLOC_PINNED|FOLDED|STARPU_MALLOC_COUNT);
#endif
}

void RUNTIME_desc_create( MORSE_desc_t *desc )
{
    int64_t lmt = desc->lmt;
    int64_t lnt = desc->lnt;
    starpu_data_handle_t *tiles;
    (void)tiles;

    desc->occurences = 1;

    /*
     * Allocate starpu_handle_t array (handlers are initialized on the fly when
     * discovered by any algorithm to save space)
     */
    desc->schedopt = (void*)calloc(lnt*lmt,sizeof(starpu_data_handle_t));
    assert(desc->schedopt);
    tiles = (starpu_data_handle_t*)(desc->schedopt);

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    if (desc->use_mat == 1 && desc->register_mat == 1){
        /*
         * Register allocated memory as CUDA pinned memory
         */
        {
            int64_t eltsze = MORSE_Element_Size(desc->dtyp);
            size_t size = (size_t)(desc->llm) * (size_t)(desc->lln) * eltsze;
            cudaError_t rc;

            /* Register the matrix as pinned memory */
            rc = cudaHostRegister( desc->mat, size, cudaHostRegisterPortable );
            if ( rc != cudaSuccess )
            {
                morse_warning("RUNTIME_desc_create(StarPU): cudaHostRegister - ", cudaGetErrorString( rc ));
            }
        }
    }
#endif
    if (desc->ooc) {
        int lastmm = desc->lm - (desc->lmt-1) * desc->mb;
        int lastnn = desc->ln - (desc->lnt-1) * desc->nb;
        int64_t eltsze = MORSE_Element_Size(desc->dtyp);
        int pagesize = getpagesize();

        if ((desc->mb * desc->nb * eltsze) % pagesize != 0
            || (lastmm   * desc->nb * eltsze) % pagesize != 0
            || (desc->mb * lastnn   * eltsze) % pagesize != 0
            || (lastmm   * lastnn   * eltsze) % pagesize != 0)
        {
            morse_error("RUNTIME_desc_create", "Matrix and tile size not suitable for out-of-core: all tiles have to be multiples of 4096. Tip : choose 'n' and 'nb' as both multiples of 32.");
            return;
        }
    }

#if defined(CHAMELEON_USE_MPI)
    /*
     * Check that we are not going over MPI tag limitations
     */
    {
        if (!_tag_mpi_initialized_) {
            int *tag_ub = NULL;
            int ok = 0;

            MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &ok);

            if ( !ok ) {
                morse_error("RUNTIME_desc_create", "MPI_TAG_UB not known by MPI");
            }

            while ( ((uintptr_t)((1UL<<tag_width) - 1) > (uintptr_t)(*tag_ub) )
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
        if ( ((uintptr_t)(lnt*lmt)) > ((uintptr_t)(1UL<<tag_sep)) ) {
            morse_fatal_error("RUNTIME_desc_create", "Too many tiles in the descriptor for MPI tags");
            return;
        }
        assert(lmt*lmt<=(1<<tag_sep));

        if ( ((uintptr_t)desc->id) >= (uintptr_t)(1UL<<(tag_width-tag_sep)) ) {
            morse_fatal_error("RUNTIME_desc_create", "Number of descriptor available in MPI mode out of stock");
            return;
        }
        assert( ((uintptr_t)desc->id) < (uintptr_t)(1UL<<(tag_width-tag_sep)) );
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

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
        if (desc->use_mat == 1 && desc->register_mat == 1){
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

/**
 * For older revision of StarPU, STARPU_MAIN_RAM is not defined
 */
#ifndef STARPU_MAIN_RAM
#define STARPU_MAIN_RAM 0
#endif

int RUNTIME_desc_getoncpu( MORSE_desc_t *desc )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    if (desc->ooc)
        /* May not even fit */
        return MORSE_SUCCESS;

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

void *RUNTIME_desc_getaddr( const MORSE_desc_t *desc, int m, int n )
{
    int64_t im = m + (desc->i / desc->mb);
    int64_t jn = n + (desc->j / desc->nb);

    starpu_data_handle_t *ptrtile = desc->schedopt;
    ptrtile += ((int64_t)desc->lmt) * jn + im;

    if (*ptrtile == NULL) {
        int home_node = -1;
        void *user_ptr = NULL;
        int myrank = desc->myrank;
        int owner  = desc->get_rankof( desc, m, n );
        int64_t eltsze = MORSE_Element_Size(desc->dtyp);
        int tempmm = (im == desc->lmt-1) ? (desc->lm - im * desc->mb) : desc->mb;
        int tempnn = (jn == desc->lnt-1) ? (desc->ln - jn * desc->nb) : desc->nb;

        if ( myrank == owner ) {
            user_ptr = desc->get_blkaddr(desc, m, n);
            if ( user_ptr != NULL ) {
                home_node = STARPU_MAIN_RAM;
            }
        }

        starpu_matrix_data_register(ptrtile, home_node, (uintptr_t) user_ptr,
                                    BLKLDD(desc, im),
                                    tempmm, tempnn, eltsze);

#ifdef HAVE_STARPU_DATA_SET_COORDINATES
        starpu_data_set_coordinates(*ptrtile, 2, m, n);
#endif

#if defined(CHAMELEON_USE_MPI)
        {
            int64_t block_ind = desc->lmt * jn + im;
            starpu_mpi_data_register(*ptrtile, (desc->id << tag_sep) | (block_ind), owner);
        }
#endif /* defined(CHAMELEON_USE_MPI) */
    }

    return *ptrtile;
}
