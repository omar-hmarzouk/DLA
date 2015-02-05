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
 * @file morse_quark.h
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.1.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 **/

/*******************************************************************************
 *  MAGMA facilities of interest to both src and magmablas directories
 **/
#ifndef _MORSE_QUARK_H_
#define _MORSE_QUARK_H_

#include <quark.h>
#include "coreblas/include/coreblas.h"
#include "runtime/quark/include/quark_blas.h"
#include "runtime/quark/include/core_blas_dag.h"

#include "control/common.h"

typedef struct quark_option_s {
    Quark_Task_Flags flags;
    Quark *quark;
} quark_option_t;

/*
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (type*)RUNTIME_desc_getaddr( desc, m, n ) )

#endif /* _MORSE_QUARK_H_ */
