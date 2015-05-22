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
 * @file descriptor.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

/**
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users
 *
 */

#include <stdlib.h>
#include <assert.h>
#include "control/common.h"
#include "control/descriptor.h"
#include "runtime.h"

static int nbdesc = 0;

/*******************************************************************************
 *  Internal static descriptor initializer
 **/
MORSE_desc_t morse_desc_init(MORSE_enum dtyp, int mb, int nb, int bsiz,
                             int lm, int ln, int i, int j,
                             int m,  int n,  int p, int q)
{
    MORSE_desc_t desc;
    desc.get_blkaddr = morse_getaddr_ccrb;
    desc.get_blkldd  = morse_getblkldd_ccrb;
    desc.get_rankof  = morse_getrankof_2d;
    // Matrix properties
    desc.dtyp = dtyp;
    // seems useless
    //desc.styp = MorseCCRB;
    desc.mb = mb;
    desc.nb = nb;
    desc.bsiz = bsiz;
    // Large matrix parameters
    desc.lm = lm;
    desc.ln = ln;
    // Large matrix derived parameters
    desc.lmt = (lm%mb==0) ? (lm/mb) : (lm/mb+1);
    desc.lnt = (ln%nb==0) ? (ln/nb) : (ln/nb+1);
    // Submatrix parameters
    desc.i = i;
    desc.j = j;
    desc.m = m;
    desc.n = n;
    // Submatrix derived parameters
    desc.mt = (m == 0) ? 0 : (i+m-1)/mb - i/mb + 1;
    desc.nt = (n == 0) ? 0 : (j+n-1)/nb - j/nb + 1;

    desc.id = nbdesc; nbdesc++;
    desc.occurences = 0;
    desc.use_mat = 1;
    desc.alloc_mat = 1;
    desc.register_mat = 1;

#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &(desc.myrank) );
#else
    desc.myrank = 0;
#endif

    // Grid size
    desc.p = p;
    desc.q = q;

    // Local dimensions in tiles
    if ( desc.myrank < (p*q) ) {
        desc.llmt = (desc.lmt + p - 1) / p;
        desc.llnt = (desc.lnt + q - 1) / q;

        // Local dimensions
        if ( ((desc.lmt-1) % p) == (desc.myrank / q) ) {
            desc.llm  = ( desc.llmt - 1 ) * mb + ((lm%mb==0) ? mb : (lm%mb));
        } else {
            desc.llm  =  desc.llmt * mb;
        }

        if ( ((desc.lnt-1) % q) == (desc.myrank % q) ) {
            desc.lln  = ( desc.llnt - 1 ) * nb + ((ln%nb==0) ? nb : (ln%nb));
        } else {
            desc.lln  =  desc.llnt * nb;
        }

        desc.llm1 = (desc.llm/mb);
        desc.lln1 = (desc.lln/nb);
    } else {
      desc.llmt = 0;
      desc.llnt = 0;
      desc.llm  = 0;
      desc.lln  = 0;
      desc.llm1 = 0;
      desc.lln1 = 0;
    }

    // Matrix address
    desc.mat = NULL;
    desc.A21 = (size_t)(desc.llm - desc.llm%mb)*(size_t)(desc.lln - desc.lln%nb);
    desc.A12 = (size_t)(           desc.llm%mb)*(size_t)(desc.lln - desc.lln%nb) + desc.A21;
    desc.A22 = (size_t)(desc.llm - desc.llm%mb)*(size_t)(           desc.lln%nb) + desc.A12;

    RUNTIME_desc_init( &desc );

    return desc;
}

/*******************************************************************************
 *  Static descriptor initializer for users
 **/
MORSE_desc_t morse_desc_init_user(MORSE_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j,
                                  int m,  int n,  int p, int q,
                                  void* (*get_blkaddr)( const MORSE_desc_t*, int, int ),
                                  int (*get_blkldd)( const MORSE_desc_t*, int ),
                                  int (*get_rankof)( const MORSE_desc_t*, int, int ))
{
    MORSE_desc_t desc;
    desc.get_blkaddr = get_blkaddr;
    desc.get_blkldd  = get_blkldd;
    desc.get_rankof  = get_rankof;
    // Matrix properties
    desc.dtyp = dtyp;
    // seems useless
    //desc.styp = MorseCCRB
    desc.mb = mb;
    desc.nb = nb;
    desc.bsiz = bsiz;
    // Large matrix parameters
    desc.lm = lm;
    desc.ln = ln;
    // Large matrix derived parameters
    desc.lmt = (lm%mb==0) ? (lm/mb) : (lm/mb+1);
    desc.lnt = (ln%nb==0) ? (ln/nb) : (ln/nb+1);
    // Submatrix parameters
    desc.i = i;
    desc.j = j;
    desc.m = m;
    desc.n = n;
    // Submatrix derived parameters
    desc.mt = (m == 0) ? 0 : (i+m-1)/mb - i/mb + 1;
    desc.nt = (n == 0) ? 0 : (j+n-1)/nb - j/nb + 1;

    desc.id = nbdesc; nbdesc++;
    desc.occurences = 0;
    desc.use_mat = 1;
    desc.alloc_mat = 1;
    desc.register_mat = 1;

#if defined(CHAMELEON_USE_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &(desc.myrank) );
#else
    desc.myrank = 0;
#endif

    // Grid size
    desc.p = p;
    desc.q = q;

    // Local dimensions in tiles
    if ( desc.myrank < (p*q) ) {
        desc.llmt = (desc.lmt + p - 1) / p;
        desc.llnt = (desc.lnt + q - 1) / q;

        // Local dimensions
        if ( ((desc.lmt-1) % p) == (desc.myrank / q) ) {
            desc.llm  = ( desc.llmt - 1 ) * mb + ((lm%mb==0) ? mb : (lm%mb));
        } else {
            desc.llm  =  desc.llmt * mb;
        }

        if ( ((desc.lnt-1) % q) == (desc.myrank % q) ) {
            desc.lln  = ( desc.llnt - 1 ) * nb + ((ln%nb==0) ? nb : (ln%nb));
        } else {
            desc.lln  =  desc.llnt * nb;
        }

        desc.llm1 = (desc.llm/mb);
        desc.lln1 = (desc.lln/nb);
    } else {
      desc.llmt = 0;
      desc.llnt = 0;
      desc.llm  = 0;
      desc.lln  = 0;
      desc.llm1 = 0;
      desc.lln1 = 0;
    }

    // Matrix address
    desc.mat = NULL;
    desc.A21 = (size_t)(desc.llm - desc.llm%mb)*(size_t)(desc.lln - desc.lln%nb);
    desc.A12 = (size_t)(           desc.llm%mb)*(size_t)(desc.lln - desc.lln%nb) + desc.A21;
    desc.A22 = (size_t)(desc.llm - desc.llm%mb)*(size_t)(           desc.lln%nb) + desc.A12;

    RUNTIME_desc_init( &desc );

    return desc;
}

/*******************************************************************************
 *  Internal static descriptor initializer for submatrices
 **/
MORSE_desc_t* morse_desc_submatrix(MORSE_desc_t *descA, int i, int j, int m, int n )
{
    MORSE_desc_t *descB = malloc(sizeof(MORSE_desc_t));
    int mb, nb;

    if ( (descA->i + i + m) > descA->lm ) {
        morse_error("morse_desc_submatrix", "The number of rows (i+m) of the submatrix doesn't fit in the parent matrix");
        assert((descA->i + i + m) > descA->lm);
    }
    if ( (descA->j + j + n) > descA->ln ) {
        morse_error("morse_desc_submatrix", "The number of rows (j+n) of the submatrix doesn't fit in the parent matrix");
        assert((descA->j + j + n) > descA->ln);
    }

    memcpy( descB, descA, sizeof(MORSE_desc_t) );
    mb = descA->mb;
    nb = descA->nb;
    // Submatrix parameters
    descB->i = descA->i + i;
    descB->j = descA->j + j;
    descB->m = m;
    descB->n = n;
    // Submatrix derived parameters
    descB->mt = (m == 0) ? 0 : (descB->i+m-1)/mb - descB->i/mb + 1;
    descB->nt = (n == 0) ? 0 : (descB->j+n-1)/nb - descB->j/nb + 1;

    RUNTIME_desc_submatrix( descB );

    return descB;
}

/*******************************************************************************
 *  Check for descriptor correctness
 **/
int morse_desc_check(MORSE_desc_t *desc)
{
    if (desc == NULL) {
        morse_error("morse_desc_check", "NULL descriptor");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (desc->mat == NULL && desc->use_mat == 1) {
        morse_error("morse_desc_check", "NULL matrix pointer");
        return MORSE_ERR_UNALLOCATED;
    }
    if (desc->dtyp != MorseRealFloat &&
        desc->dtyp != MorseRealDouble &&
        desc->dtyp != MorseComplexFloat &&
        desc->dtyp != MorseComplexDouble  ) {
        morse_error("morse_desc_check", "invalid matrix type");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    if (desc->mb <= 0 || desc->nb <= 0) {
        morse_error("morse_desc_check", "negative tile dimension");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    if (desc->bsiz < desc->mb*desc->nb) {
        morse_error("morse_desc_check", "tile memory size smaller than the product of dimensions");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    if (desc->lm <= 0 || desc->ln <= 0) {
        morse_error("morse_desc_check", "negative matrix dimension");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    if ((desc->lm < desc->m) || (desc->ln < desc->n)) {
        morse_error("morse_desc_check", "matrix dimensions larger than leading dimensions");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    if ((desc->i > 0 && desc->i >= desc->lm) || (desc->j > 0 && desc->j >= desc->ln)) {
        morse_error("morse_desc_check", "beginning of the matrix out of scope");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    if (desc->i+desc->m > desc->lm || desc->j+desc->n > desc->ln) {
        morse_error("morse_desc_check", "submatrix out of scope");
        return MORSE_ERR_ILLEGAL_VALUE;
    }
    return MORSE_SUCCESS;
}

/*******************************************************************************
 *
 **/
int morse_desc_mat_alloc( MORSE_desc_t *desc )
{

    size_t size = (size_t)(desc->llm) * (size_t)(desc->lln)
      * (size_t)MORSE_Element_Size(desc->dtyp);
    if ((desc->mat = malloc(size)) == NULL) {
        morse_error("morse_desc_mat_alloc", "malloc() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }

    RUNTIME_desc_create(desc);

    return MORSE_SUCCESS;
}

/*******************************************************************************
 *
 **/
int morse_desc_mat_free( MORSE_desc_t *desc )
{

    RUNTIME_desc_destroy( desc );

    if (desc->mat != NULL  &&
        desc->use_mat == 1 &&
        desc->alloc_mat == 1) {
#ifndef CHAMELEON_SIMULATION
        free(desc->mat);
#endif
        desc->mat = NULL;
    }
    return MORSE_SUCCESS;
}

/** ***************************************************************************
 *
 * @ingroup Descriptor
 *
 *  MORSE_Desc_Create - Create matrix descriptor.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg MorseRealFloat:     single precision real (S),
 *          @arg MorseRealDouble:    double precision real (D),
 *          @arg MorseComplexFloat:  single precision complex (C),
 *          @arg MorseComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in bytes including padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Desc_Create(MORSE_desc_t **desc, void *mat, MORSE_enum dtyp, int mb, int nb, int bsiz,
                      int lm, int ln, int i, int j, int m, int n, int p, int q)
{
    MORSE_context_t *morse;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Desc_Create", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Allocate memory and initialize the descriptor */
    *desc = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    if (*desc == NULL) {
        morse_error("MORSE_Desc_Create", "malloc() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }
    **desc = morse_desc_init(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q);

    if (mat == NULL) {

#ifdef CHAMELEON_SIMULATION
        (**desc).mat = (void*) 1;
#else
        /* TODO: a call to morse_desc_mat_alloc should be made, but require to
        move the call to RUNTIME_desc_create within the function */
        size_t size = (size_t)((*desc)->llm) * (size_t)((*desc)->lln)
            * (size_t)MORSE_Element_Size((*desc)->dtyp);

        if (((**desc).mat = malloc(size)) == NULL) {
            morse_error("MORSE_Desc_Create", "malloc() failed");
            return MORSE_ERR_OUT_OF_RESOURCES;
        }
#endif

    } else {
        (**desc).mat = mat;
        /* memory of the matrix is handle by users */
        (**desc).alloc_mat = 0;
    }

    /* Create scheduler structure like registering data */
    RUNTIME_desc_create( *desc );

    status = morse_desc_check(*desc);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_Desc_Create", "invalid descriptor");
        return status;
    }

    return MORSE_SUCCESS;
}

/** ***************************************************************************
 *
 * @ingroup Descriptor
 *
 *  MORSE_Desc_Create_User - Create matrix descriptor for general applications.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg MorseRealFloat:     single precision real (S),
 *          @arg MorseRealDouble:    double precision real (D),
 *          @arg MorseComplexFloat:  single precision complex (C),
 *          @arg MorseComplexDouble: double precision complex (Z).
 *
 * @param[in] nb
 *          Number of rows and columns in a tile.
 *
 * @param[in] m
 *          Number of rows of the entire matrix.
 *
 * @param[in] n
 *          Number of columns of the entire matrix.
 *
  * @param[in] p
 *          2d-block cyclic partitioning, number of tiles in rows.
 *
 * @param[in] q
 *          2d-block cyclic partitioning, number of tiles in columns.
 *
 * @param[in] (*get_blkaddr)( const MORSE_desc_t *A, int m, int n)
 *          A function which return the address of the data corresponding to
 *          the tile A(m,n).
 *
 * @param[in] (*get_blkldd)( const MORSE_desc_t *A, int m)
 *          A function that return the leading dimension of the tile A(m,*).
 *
 * @param[in] (*get_rankof)( const MORSE_desc_t *A, int m, int n)
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Desc_Create_User(MORSE_desc_t **desc, void *mat, MORSE_enum dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q,
                           void* (*get_blkaddr)( const MORSE_desc_t*, int, int ),
                           int (*get_blkldd)( const MORSE_desc_t*, int ),
                           int (*get_rankof)( const MORSE_desc_t*, int, int ))
{
    MORSE_context_t *morse;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Desc_Create", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    /* Allocate memory and initialize the descriptor */
    *desc = (MORSE_desc_t*)malloc(sizeof(MORSE_desc_t));
    if (*desc == NULL) {
        morse_error("MORSE_Desc_Create", "malloc() failed");
        return MORSE_ERR_OUT_OF_RESOURCES;
    }
    **desc = morse_desc_init_user(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q,
        get_blkaddr, get_blkldd, get_rankof);


    /* if the user gives a pointer to the overall data (tiles) we can use it */
    (**desc).use_mat = (mat == NULL) ? 0 : 1;

    /* memory of the matrix is handle by users */
    (**desc).alloc_mat = 0;

    /* users data can have multiple forms: let him register tiles */
    (**desc).register_mat = 0;

    (**desc).mat = mat;

    /* Create scheduler structure like registering data */
    RUNTIME_desc_create( *desc );

    status = morse_desc_check(*desc);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_Desc_Create", "invalid descriptor");
        return status;
    }

    return MORSE_SUCCESS;
}

/** ***************************************************************************
 *
 * @ingroup Descriptor
 *
 *  MORSE_Desc_Destroy - Destroys matrix descriptor.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Desc_Destroy(MORSE_desc_t **desc)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_error("MORSE_Desc_Destroy", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    if (*desc == NULL) {
        morse_error("MORSE_Desc_Destroy", "attempting to destroy a NULL descriptor");
        return MORSE_ERR_UNALLOCATED;
    }

    morse_desc_mat_free( *desc );
    free(*desc);
    *desc = NULL;
    return MORSE_SUCCESS;
}

/** ***************************************************************************
 *
 * @ingroup Descriptor
 *
 *  MORSE_Desc_Acquire - Ensures that all data of the descriptor are
 *  up-to-date.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Desc_Acquire (MORSE_desc_t  *desc) {
  return RUNTIME_desc_acquire( desc );
}

/** ***************************************************************************
 *
 * @ingroup Descriptor
 *
 *  MORSE_Desc_Release - Release the data of the descriptor acquired by the
 *  application. Should be called if MORSE_Desc_Acquire has been called on the
 *  descriptor and if you do not need to access to its data anymore.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Desc_Release (MORSE_desc_t  *desc) {
  return RUNTIME_desc_release( desc );
}

/** ***************************************************************************
 *
 * @ingroup Descriptor
 *
 *  MORSE_Desc_Getoncpu - Apply an acquire and a release on the data of the
 *  descriptors so that its values are up-to-date in the main memory.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *****************************************************************************/
int MORSE_Desc_Getoncpu(MORSE_desc_t  *desc) {
  return RUNTIME_desc_getoncpu( desc );
}
