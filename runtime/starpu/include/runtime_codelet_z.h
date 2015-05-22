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
 * @file runtime_codelet_z.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 * @precisions normal z -> c d s
 *
 **/

#ifndef _CODELETS_Z_H_
#define _CODELETS_Z_H_

#include <stdio.h>
#include "runtime/starpu/include/runtime_codelets.h"

/*
 * Management functions
 */
ZCODELETS_HEADER(tile_zero)

/*
 * BLAS 1 functions
 */
ZCODELETS_HEADER(axpy)

/*
 * BLAS 3 functions
 */
ZCODELETS_HEADER(gemm)
ZCODELETS_HEADER(hemm)
ZCODELETS_HEADER(her2k)
ZCODELETS_HEADER(herk)
ZCODELETS_HEADER(symm)
ZCODELETS_HEADER(syr2k)
ZCODELETS_HEADER(syrk)
ZCODELETS_HEADER(trmm)
ZCODELETS_HEADER(trsm)

/*
 * LAPACK functions
 */
ZCODELETS_HEADER(gelqt)
ZCODELETS_HEADER(geqrt)
ZCODELETS_HEADER(gessm)
ZCODELETS_HEADER(gessq)
ZCODELETS_HEADER(getrf)
ZCODELETS_HEADER(getrf_incpiv)
ZCODELETS_HEADER(getrf_nopiv)
ZCODELETS_HEADER(lauum)
ZCODELETS_HEADER(potrf)
ZCODELETS_HEADER(ssssm)
ZCODELETS_HEADER(syssq)
ZCODELETS_HEADER(trasm)
ZCODELETS_HEADER(trssq)
ZCODELETS_HEADER(trtri)
ZCODELETS_HEADER(tslqt)
ZCODELETS_HEADER(tsmlq)
ZCODELETS_HEADER(tsmqr)
ZCODELETS_HEADER(tsqrt)
ZCODELETS_HEADER(tstrf)
ZCODELETS_HEADER(ttlqt)
ZCODELETS_HEADER(ttmlq)
ZCODELETS_HEADER(ttmqr)
ZCODELETS_HEADER(ttqrt)
ZCODELETS_HEADER(unmlq)
ZCODELETS_HEADER(unmqr)

/*
 * Auxiliary functions
 */
ZCODELETS_HEADER(geadd)
ZCODELETS_HEADER(lacpy)
ZCODELETS_HEADER(lange)
ZCODELETS_HEADER(lange_max)
ZCODELETS_HEADER(lansy)
ZCODELETS_HEADER(lantr)
ZCODELETS_HEADER(laset)
ZCODELETS_HEADER(laset2)
ZCODELETS_HEADER(plssq)
ZCODELETS_HEADER(plssq2)

/*
 * MIXED PRECISION functions
 */
ZCODELETS_HEADER(lag2c)

/*
 * DZ functions
 */
ZCODELETS_HEADER(asum)

/*
 * CPU only functions
 */
ZCODELETS_HEADER(plrnt)

#if defined(PRECISION_z) || defined(PRECISION_c)
ZCODELETS_HEADER(hessq)
ZCODELETS_HEADER(lanhe)
ZCODELETS_HEADER(plghe)
ZCODELETS_HEADER(sytrf_nopiv)
#endif
ZCODELETS_HEADER(plgsy)

#endif /* _CODELETS_Z_H_ */
