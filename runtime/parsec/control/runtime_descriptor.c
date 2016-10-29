/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
#include <stdlib.h>
#include "runtime/parsec/include/morse_parsec.h"
#include <dague/data.h>

void RUNTIME_user_tag_size(int user_tag_width, int user_tag_sep) {
  (void)user_tag_width;
  (void)user_tag_sep;
}

void *RUNTIME_mat_alloc( size_t size)
{
    return malloc(size);
}

void RUNTIME_mat_free( void *mat, size_t size)
{
    (void)size;
    free(mat);
    return;
}

struct morse_parsec_desc_s {
    dague_ddesc_t  super;
    MORSE_desc_t  *desc;
    dague_data_t **data_map;
};

static void
morse_parsec_key_to_coordinates(dague_ddesc_t *ddesc, dague_data_key_t key,
                                int *m, int *n)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)ddesc;
    MORSE_desc_t *mdesc = pdesc->desc;
    int _m, _n;

    _m = key % mdesc->lmt;
    _n = key / mdesc->lmt;
    *m = _m - mdesc->i / mdesc->mb;
    *n = _n - mdesc->j / mdesc->nb;
}

static dague_data_key_t
morse_parsec_data_key(dague_ddesc_t *ddesc, ...)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)ddesc;
    MORSE_desc_t *mdesc = pdesc->desc;
    va_list ap;
    int m, n;

    /* Get coordinates */
    va_start(ap, ddesc);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Offset by (i,j) to translate (m,n) in the global matrix */
    m += mdesc->i / mdesc->mb;
    n += mdesc->j / mdesc->nb;

    return ((n * mdesc->lmt) + m);
}

static uint32_t
morse_parsec_rank_of(dague_ddesc_t *ddesc, ...)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)ddesc;
    MORSE_desc_t *mdesc = pdesc->desc;
    va_list ap;
    int m, n;

    /* Get coordinates */
    va_start(ap, ddesc);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Offset by (i,j) to translate (m,n) in the global matrix */
    m += mdesc->i / mdesc->mb;
    n += mdesc->j / mdesc->nb;

    return mdesc->get_rankof( mdesc, m, n );
}

static uint32_t
morse_parsec_rank_of_key(dague_ddesc_t *ddesc, dague_data_key_t key)
{
    int m, n;
    morse_parsec_key_to_coordinates(ddesc, key, &m, &n);
    return morse_parsec_rank_of(ddesc, m, n);
}

static int32_t
morse_parsec_vpid_of(dague_ddesc_t *ddesc, ...)
{
    return 0;
}

static int32_t
morse_parsec_vpid_of_key(dague_ddesc_t *ddesc, dague_data_key_t key)
{
    int m, n;
    morse_parsec_key_to_coordinates(ddesc, key, &m, &n);
    return morse_parsec_vpid_of(ddesc, m, n);
}

static dague_data_t*
morse_parsec_data_of(dague_ddesc_t *ddesc, ...)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)ddesc;
    MORSE_desc_t *mdesc = pdesc->desc;
    va_list ap;
    int m, n;

    /* Get coordinates */
    va_start(ap, ddesc);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Offset by (i,j) to translate (m,n) in the global matrix */
    m += mdesc->i / mdesc->mb;
    n += mdesc->j / mdesc->nb;

#if defined(CHAMELEON_USE_MPI)
    /* TODO: change displacement in data_map when in distributed */
    assert( mdesc->nodes == 1 );
#endif
    return dague_data_create( pdesc->data_map + n * mdesc->lmt + m, ddesc,
                              morse_parsec_data_key( ddesc, m, n ),
                              mdesc->get_blkaddr( mdesc, m, n ),
                              mdesc->bsiz * MORSE_Element_Size(mdesc->dtyp) );
}

static dague_data_t*
morse_parsec_data_of_key(dague_ddesc_t *ddesc, dague_data_key_t key)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)ddesc;
    MORSE_desc_t *mdesc = pdesc->desc;
    int m, n;
    morse_parsec_key_to_coordinates(ddesc, key, &m, &n);

#if defined(CHAMELEON_USE_MPI)
    /* TODO: change displacement in data_map when in distributed */
    assert( mdesc->nodes == 1 );
#endif
    return dague_data_create( pdesc->data_map + key, ddesc, key,
                              mdesc->get_blkaddr( mdesc, m, n ),
                              mdesc->bsiz * MORSE_Element_Size(mdesc->dtyp) );
}

#ifdef DAGUE_PROF_TRACE
static int
morse_parsec_key_to_string(dague_ddesc_t *ddesc, dague_data_key_t key, char * buffer, uint32_t buffer_size)
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)ddesc;
    MORSE_desc_t *mdesc = pdesc->desc;
    int m, n, res;
    morse_parsec_key_to_coordinates(ddesc, key, &m, &n);
    res = snprintf(buffer, buffer_size, "(%d, %d)", m, n);
    if (res < 0)
    {
        printf("error in key_to_string for tile (%u, %u) key: %u\n", m, n, datakey);
    }
    return res;
}
#endif

void RUNTIME_desc_init( MORSE_desc_t *mdesc )
{
    (void)mdesc;
    return;
}

void RUNTIME_desc_create( MORSE_desc_t *mdesc )
{
    dague_ddesc_t       *ddesc;
    morse_parsec_desc_t *pdesc;
    int comm_size;

    pdesc = malloc( sizeof(morse_parsec_desc_t) );
    ddesc = (dague_ddesc_t*)pdesc;

    /* Super setup */
    RUNTIME_comm_size( &comm_size );
    ddesc->nodes  = comm_size;
    ddesc->myrank = mdesc->myrank;

    ddesc->data_key    = morse_parsec_data_key;
    ddesc->rank_of     = morse_parsec_rank_of;
    ddesc->rank_of_key = morse_parsec_rank_of_key;
    ddesc->data_of     = morse_parsec_data_of;
    ddesc->data_of_key = morse_parsec_data_of_key;
    ddesc->vpid_of     = morse_parsec_vpid_of;
    ddesc->vpid_of_key = morse_parsec_vpid_of_key;
#if defined(DAGUE_PROF_TRACE)
    {
        int rc;
        ddesc->key_to_string = morse_parsec_key_to_string;
        ddesc->key           = NULL;
        rc = asprintf(&(ddesc->key_dim), "(%d, %d)", mdesc->lmt, mdesc->lnt);
        (void)rc;
    }
#endif
    ddesc->memory_registration_status = MEMORY_STATUS_UNREGISTERED;

    pdesc->data_map = calloc( mdesc->lmt * mdesc->lnt, sizeof(dague_data_t*) );

    /* Double linking */
    pdesc->desc     = mdesc;
    mdesc->schedopt = pdesc;

    /* /\* Overwrite the leading dimensions to store the padding *\/ */
    /* mdesc->llm = mdesc->mb * mdesc->lmt; */
    /* mdesc->lln = mdesc->nb * mdesc->lnt; */
    return;
}

void RUNTIME_desc_destroy( MORSE_desc_t *mdesc )
{
    morse_parsec_desc_t *pdesc = (morse_parsec_desc_t*)(mdesc->schedopt);

    if ( pdesc->data_map != NULL ) {
        dague_data_t **data = pdesc->data_map;
        int nb_local_tiles = mdesc->lmt * mdesc->lnt;
        int i;

        for(i=0; i<nb_local_tiles; i++, data++)
        {
            dague_data_destroy( *data );
        }

        free( pdesc->data_map );
        pdesc->data_map = NULL;
    }
    free(pdesc);
    return;
}

void RUNTIME_desc_submatrix( MORSE_desc_t *desc )
{
    (void)desc;
    return;
}

int RUNTIME_desc_acquire( MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

int RUNTIME_desc_release( MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

int RUNTIME_desc_getoncpu( MORSE_desc_t *desc )
{
    (void)desc;
    return MORSE_SUCCESS;
}

void *RUNTIME_desc_getaddr( MORSE_desc_t *desc, int m, int n )
{
    assert(0); /* This should not be called because we also need the handle to match the address we need. */
    return desc->get_blkaddr( desc, m, n );
}
