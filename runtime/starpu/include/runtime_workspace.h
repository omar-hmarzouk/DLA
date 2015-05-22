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
 * @file runtime_workspace.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Cedric Augonnet
 * @date 2011-06-01
 *
 **/

#ifndef _MORSE_STARPU_WORKSPACE_H_
#define _MORSE_STARPU_WORKSPACE_H_

/* 
 * Allocate workspace in host memory: CPU for any worker 
 * or allocate workspace in worker's memory: main memory for cpu workers,
 * and embedded memory for CUDA devices. 
 */
#define MORSE_HOST_MEM    0
#define MORSE_WORKER_MEM  1

struct morse_starpu_ws_s {
    size_t size;
    int    memory_location;
    int    which_workers;
    void  *workspaces[STARPU_NMAXWORKERS];
};

typedef struct morse_starpu_ws_s MORSE_starpu_ws_t;

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. MORSE_CUDA|MORSE_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in worker's memory (MORSE_HOST_MEM or MORSE_WORKER_MEM). This function
 * returns 0 upon successful completion. 
 */
int   RUNTIME_starpu_ws_alloc   ( MORSE_starpu_ws_t **workspace, size_t size, int which_workers, int memory_location);
int   RUNTIME_starpu_ws_free    ( MORSE_starpu_ws_t  *workspace);
void *RUNTIME_starpu_ws_getlocal( MORSE_starpu_ws_t  *workspace);

#endif /* _MORSE_STARPU_WORKSPACE_H_ */
