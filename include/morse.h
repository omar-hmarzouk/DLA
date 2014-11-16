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
 *  @file morse.h
 *
 *  MORSE main header
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 1.1.0
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @author Cedric Castagnede
 *  @date 2012-09-15
 *
 **/
#ifndef _MORSE_H_
#define _MORSE_H_

#define MORSE_VERSION_MAJOR 1
#define MORSE_VERSION_MINOR 1
#define MORSE_VERSION_MICRO 0

#define CHAMELEON_VERSION_MAJOR 1
#define CHAMELEON_VERSION_MINOR 1
#define CHAMELEON_VERSION_MICRO 0


/** ****************************************************************************
 * MORSE types and constants
 **/
//#include "morse_mangling.h"
#include "morse_types.h"
#include "morse_struct.h"
#include "morse_constants.h"


/** ****************************************************************************
 * RUNTIME Functions
 */
#include "runtime.h"


/** ****************************************************************************
 * For Simulation mode
 */
#include "morse_simulate.h"


/** ****************************************************************************
 * Set of routines which can be useful fo users
 */
#include "context.h"
#include "descriptor.h"


/** ****************************************************************************
 * MORSE Functions
 */
#ifdef __cplusplus
extern "C" {
#endif

/* Auxiliary */
int MORSE_Version        (int *ver_major, int *ver_minor, int *ver_micro);
int MORSE_Init           (int nworkers, int ncudas);
int MORSE_InitPar        (int nworkers, int ncudas, int nthreads_per_worker);
int MORSE_Finalize       (void);
int MORSE_My_Mpi_Rank    (void);
int MORSE_Lapack_to_Tile (void *Af77, int LDA, MORSE_desc_t *A);
int MORSE_Tile_to_Lapack (MORSE_desc_t *A, void *Af77, int LDA);

/* Descriptor */
int MORSE_Desc_Create  (MORSE_desc_t **desc, void *mat, MORSE_enum dtyp,
                        int mb, int nb, int bsiz, int lm, int ln,
                        int i, int j, int m, int n, int p, int q);
int MORSE_Desc_Create_User(MORSE_desc_t **desc, void *mat, MORSE_enum dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q,
                           void* (*get_blkaddr)( const MORSE_desc_t*, int, int ),
                           int (*get_blkldd)( const MORSE_desc_t*, int ),
                           int (*get_rankof)( const MORSE_desc_t*, int, int ));
int MORSE_Desc_Destroy (MORSE_desc_t **desc);
int MORSE_Desc_Acquire (MORSE_desc_t  *desc);
int MORSE_Desc_Release (MORSE_desc_t  *desc);
int MORSE_Desc_Getoncpu(MORSE_desc_t  *desc);

/* Workspaces */
int MORSE_Dealloc_Workspace (MORSE_desc_t **desc);

/* Options */
int MORSE_Enable  (MORSE_enum option);
int MORSE_Disable (MORSE_enum option);
int MORSE_Set     (MORSE_enum param, int  value);
int MORSE_Get     (MORSE_enum param, int *value);

/* Sequences */
int MORSE_Sequence_Create  (MORSE_sequence_t **sequence);
int MORSE_Sequence_Destroy (MORSE_sequence_t *sequence);
int MORSE_Sequence_Wait    (MORSE_sequence_t *sequence);

#ifdef __cplusplus
}
#endif

#include "morse_z.h"
#include "morse_c.h"
#include "morse_d.h"
#include "morse_s.h"
#include "morse_zc.h"
#include "morse_ds.h"

#endif
