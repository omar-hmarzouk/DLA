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
 * @file common.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 **/

/*******************************************************************************
 *  MORSE facilities of interest to both MORSE core developer
 *  and also of interest to MORSE community contributor.
 **/
#ifndef _MORSE_COMMON_H_
#define _MORSE_COMMON_H_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif


/** ****************************************************************************
 * Implementation headers
 **/
#if defined(CHAMELEON_USE_CUDA)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas_v2.h>
#else
#include <cublas.h>
#endif
#endif

#if defined(CHAMELEON_USE_OPENCL)
#include <OpenCL/cl.h>
#endif

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif


/** ****************************************************************************
 * Linear Algebra headers
 **/
#if defined(CHAMELEON_USE_MAGMA)
#include <magma.h>
#endif


/** ****************************************************************************
 *  Line to avoid conflict with magma, because, we don't know why
 *  but lapacke provide a wrong interface of lapack in fortran
 **/
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif
#include "coreblas/include/lapacke.h"
#include "coreblas/include/coreblas.h"

#include "morse.h"

#include "control/global.h"
#include "control/auxiliary.h"
#include "control/tile.h"
#include "control/async.h"
#include "control/bulge.h"

/** ****************************************************************************
 *  Determine FORTRAN names
 **/
#if defined(ADD_)
#define MORSE_FNAME(lcname, UCNAME)        morse_##lcname##_
#define MORSE_TILE_FNAME(lcname, UCNAME)   morse_##lcname##_tile_
#define MORSE_ASYNC_FNAME(lcname, UCNAME)  morse_##lcname##_tile_async_
#define MORSE_WS_FNAME(lcname, UCNAME)     morse_alloc_workspace_##lcname##_
#define MORSE_WST_FNAME(lcname, UCNAME)    morse_alloc_workspace_##lcname##_tile_
#elif defined(NOCHANGE)
#define MORSE_FNAME(lcname, UCNAME)        morse_##lcname
#define MORSE_TILE_FNAME(lcname, UCNAME)   morse_##lcname##_tile
#define MORSE_ASYNC_FNAME(lcname, UCNAME)  morse_##lcname##_tile_async
#define MORSE_WS_FNAME(lcname, UCNAME)     morse_alloc_workspace_##lcname
#define MORSE_WST_FNAME(lcname, UCNAME)    morse_alloc_workspace_##lcname##_tile
#elif defined(UPCASE)
#define MORSE_FNAME(lcname, UCNAME)        MORSE_##UCNAME
#define MORSE_TILE_FNAME(lcname, UCNAME)   MORSE_##UCNAME##_TILE
#define MORSE_ASYNC_FNAME(lcname, UCNAME)  MORSE_##UCNAME##_TILE_ASYNC
#define MORSE_WS_FNAME(lcname, UCNAME)     MORSE_ALLOC_WORKSPACE_##UCNAME
#define MORSE_WST_FNAME(lcname, UCNAME)    MORSE_ALLOC_WORKSPACE_##UCNAME##_TILE
#endif


/*******************************************************************************
 *  Global shortcuts
 **/
#define MORSE_RANK        morse_rank(morse)
#define MORSE_SIZE        morse->world_size
#define MORSE_GRPSIZE     morse->group_size
#define MORSE_NB          morse->nb
#define MORSE_IB          morse->ib
#define MORSE_NBNBSIZE    morse->nbnbsize
#define MORSE_IBNBSIZE    morse->ibnbsize
#define MORSE_SCHEDULING  morse->scheduling
#define MORSE_RHBLK       morse->rhblock
#define MORSE_TRANSLATION morse->translation
#define MORSE_PARALLEL    morse->parallel_enabled
#define MORSE_PROFILING   morse->profiling_enabled
#if defined(CHAMELEON_USE_MPI)
#define MORSE_MPI_RANK    morse->my_mpi_rank
#define MORSE_MPI_SIZE    morse->mpi_comm_size
#endif

/*******************************************************************************
 *  Activate copy of diagonal tile (StarPU only) for some tile algorithms (pz)
 **/
#if defined(CHAMELEON_SCHED_STARPU)
#define CHAMELEON_COPY_DIAG
#endif

/*******************************************************************************
 *  IPT internal define
 **/
#define MorseIPT_NoDep   0
#define MorseIPT_Panel   1
#define MorseIPT_All     2


/*******************************************************************************
 *  Global array of LAPACK constants
 **/
extern char *morse_lapack_constants[];
#define morse_lapack_const(morse_const) morse_lapack_constants[morse_const][0]

#ifdef __cplusplus
extern "C" {
#endif

#include "control/compute_s.h"
#include "control/compute_d.h"
#define COMPLEX
#include "control/compute_c.h"
#include "control/compute_z.h"
#undef COMPLEX

/*
void morse_pdlag2s(MORSE_context_t *morse);
void morse_pzlag2c(MORSE_context_t *morse);
void morse_pslag2d(MORSE_context_t *morse);
void morse_pclag2z(MORSE_context_t *morse);
*/

#ifdef __cplusplus
}
#endif

#endif
