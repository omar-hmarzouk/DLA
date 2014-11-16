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
 * @file morse_starpu.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 1.1.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 **/

/******************************************************************************/

/*
 *  MORSE facilities of interest to both src and magmablas directories
 **/
#ifndef _MORSE_STARPU_H_
#define _MORSE_STARPU_H_

#if defined(MAGMAMORSE_USE_MPI)
#include <starpu_mpi.h>
#else
#include <starpu.h>
#endif

#include <starpu_profiling.h>

#if defined(MAGMAMORSE_USE_CUDA)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>
#endif

#include "common.h"
#include "codelets.h"
#include "profiling.h"
#include "codelet_profile.h"
#include "workspace.h"

typedef struct starpu_conf starpu_conf_t;

/******************************************************************************/

/*
 * MPI Redefinitions
 */
#if defined(MAGMAMORSE_USE_MPI)
#undef STARPU_REDUX
#define starpu_insert_task(...) starpu_mpi_insert_task(MPI_COMM_WORLD, __VA_ARGS__)
#endif

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)RUNTIME_desc_getaddr( desc, m, n ) )

void RUNTIME_set_reduction_methods(starpu_data_handle_t handle, MORSE_enum dtyp);

#endif /* _MORSE_STARPU_H_ */
