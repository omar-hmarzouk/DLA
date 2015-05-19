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
 * @file runtime.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 **/
#ifndef _RUNTIME_H_
#define _RUNTIME_H_

#include "morse_struct.h"

/*******************************************************************************
 * RUNTIME Async
 **/
int   RUNTIME_sequence_create  (MORSE_context_t*, MORSE_sequence_t*);
int   RUNTIME_sequence_destroy (MORSE_context_t*, MORSE_sequence_t*);
int   RUNTIME_sequence_wait    (MORSE_context_t*, MORSE_sequence_t*);
void  RUNTIME_sequence_flush   (void* quark, MORSE_sequence_t*, MORSE_request_t*, int);

/*******************************************************************************
 * RUNTIME Context
 **/
void  RUNTIME_context_create  (MORSE_context_t*);
void  RUNTIME_context_destroy (MORSE_context_t*);
void  RUNTIME_enable          (MORSE_enum);
void  RUNTIME_disable         (MORSE_enum);

/*******************************************************************************
 * RUNTIME Control
 **/
int   RUNTIME_rank               (MORSE_context_t*);
int   RUNTIME_init_scheduler     (MORSE_context_t*, int, int, int);
void  RUNTIME_finalize_scheduler (MORSE_context_t*);
void  RUNTIME_barrier            (MORSE_context_t*);
void  RUNTIME_pause              (MORSE_context_t*);
void  RUNTIME_resume             (MORSE_context_t*);

/*******************************************************************************
 * RUNTIME Descriptor
 **/
void  RUNTIME_desc_init      (MORSE_desc_t*);
void  RUNTIME_desc_create    (MORSE_desc_t*);
void  RUNTIME_desc_destroy   (MORSE_desc_t*);
void  RUNTIME_desc_submatrix (MORSE_desc_t*);
void* RUNTIME_desc_getaddr   (MORSE_desc_t*, int, int);
/* Acquire in main memory an up-to-date copy of the data described by the descriptor for read-write access. */
int   RUNTIME_desc_acquire   (MORSE_desc_t*);
/* Release the data described by the descriptor to be used by the StarPU tasks again. */
int   RUNTIME_desc_release   (MORSE_desc_t*);
int   RUNTIME_desc_getoncpu  (MORSE_desc_t*);

/*******************************************************************************
 * RUNTIME Options
 **/
void  RUNTIME_options_init     (MORSE_option_t*, MORSE_context_t*, MORSE_sequence_t*, MORSE_request_t*);
void  RUNTIME_options_finalize (MORSE_option_t*, MORSE_context_t *);
int   RUNTIME_options_ws_alloc (MORSE_option_t*, size_t, size_t);
int   RUNTIME_options_ws_free  (MORSE_option_t*);
/* int   RUNTIME_options_ws_gethost   (MORSE_option_t*); */
/* int   RUNTIME_options_ws_getdevice (MORSE_option_t*); */

/*******************************************************************************
 * RUNTIME Locality
 **/
void RUNTIME_zlocality_allrestore ();
void RUNTIME_clocality_allrestore ();
void RUNTIME_dlocality_allrestore ();
void RUNTIME_slocality_allrestore ();

void RUNTIME_zlocality_allrestrict(uint32_t);
void RUNTIME_zlocality_onerestrict(MORSE_kernel_t, uint32_t);
void RUNTIME_zlocality_onerestore (MORSE_kernel_t);

void RUNTIME_clocality_allrestrict(uint32_t);
void RUNTIME_clocality_onerestrict(MORSE_kernel_t, uint32_t);
void RUNTIME_clocality_onerestore (MORSE_kernel_t);

void RUNTIME_dlocality_allrestrict(uint32_t);
void RUNTIME_dlocality_onerestrict(MORSE_kernel_t, uint32_t);
void RUNTIME_dlocality_onerestore (MORSE_kernel_t);

void RUNTIME_slocality_allrestrict(uint32_t);
void RUNTIME_slocality_onerestrict(MORSE_kernel_t, uint32_t);
void RUNTIME_slocality_onerestore (MORSE_kernel_t);

/*******************************************************************************
 * RUNTIME Profiling
 **/
void  RUNTIME_schedprofile_display ();
void  RUNTIME_kernelprofile_display();
double RUNTIME_get_time();

void RUNTIME_start_profiling();
void RUNTIME_stop_profiling();

#if defined(PRECISION_z)
void RUNTIME_zdisplay_allprofile ();
void RUNTIME_zdisplay_oneprofile (MORSE_kernel_t);
#endif

#if defined(PRECISION_c)
void RUNTIME_cdisplay_allprofile ();
void RUNTIME_cdisplay_oneprofile (MORSE_kernel_t);
#endif

#if defined(PRECISION_d)
void RUNTIME_ddisplay_allprofile ();
void RUNTIME_ddisplay_oneprofile (MORSE_kernel_t);
#endif

#if defined(PRECISION_s)
void RUNTIME_sdisplay_allprofile ();
void RUNTIME_sdisplay_oneprofile (MORSE_kernel_t);
#endif

/*******************************************************************************
 * RUNTIME Kernels
 **/
#include "runtime_z.h"
#include "runtime_d.h"
#include "runtime_c.h"
#include "runtime_s.h"
#include "runtime_zc.h"
#include "runtime_ds.h"

void MORSE_TASK_ztile_zero(MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_dtile_zero(MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_ctile_zero(MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           MORSE_desc_t *A, int Am, int An, int lda);
void MORSE_TASK_stile_zero(MORSE_option_t *options,
                           int X1, int X2, int Y1, int Y2,
                           MORSE_desc_t *A, int Am, int An, int lda);

/*
 * Mark a data as unused after this call
 */
void MORSE_TASK_dataflush(MORSE_option_t *options,
                          MORSE_desc_t *A, int Am, int An);
void MORSE_TASK_dataflush_all();

#endif
