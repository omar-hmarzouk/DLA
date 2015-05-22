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
 *  @file codelet_zpotrf.c
 *
 *  MAGMA codelets kernel
 *  MAGMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 *  @version 2.3.1
 *  @author Mathieu Faverge
 *  @author Cedric Augonnet
 *  @date 2011-06-01
 *  @precisions normal z -> c d s
 *
 **/
#include "runtime/starpu/include/morse_starpu.h"
#include "runtime/starpu/include/runtime_codelet_z.h"

CHAMELEON_CL_CB(zasum,         starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      M*N);
CHAMELEON_CL_CB(zaxpy,         starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[1]), 0,                                      M);
CHAMELEON_CL_CB(zgeadd,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      M*N);
CHAMELEON_CL_CB(zgelqt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      (4./3.)*M*N*K);
CHAMELEON_CL_CB(zgemm,         starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_ny(task->handles[2]), starpu_matrix_get_ny(task->handles[0]),     2. *M*N*K); /* If A^t, computation is wrong */
CHAMELEON_CL_CB(zgeqrt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      (4./3.)*M*M*N);
CHAMELEON_CL_CB(zgessm,        starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_nx(task->handles[2]),     2. *M*N*K);
CHAMELEON_CL_CB(zgessq,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), 0,                                      4.*M*N);
CHAMELEON_CL_CB(zgetrf,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (2./3.)*M*N*K);
CHAMELEON_CL_CB(zgetrf_incpiv, starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (2./3.)*M*N*K);
CHAMELEON_CL_CB(zgetrf_nopiv,  starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (2./3.)*M*N*K);
#if defined(PRECISION_z) || defined(PRECISION_c)
CHAMELEON_CL_CB(zhemm,         starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_ny(task->handles[2]), 0,                                          2.*M*M *N);
CHAMELEON_CL_CB(zher2k,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                     ( 1.+2.*M*N)*M);
CHAMELEON_CL_CB(zherk,         starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                     ( 1.+   M)*M*N);
#endif
CHAMELEON_CL_CB(zlacpy,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zlange,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zlaset,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zlaset2,       starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zlauum,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (1./3.)*M* M*M);
#if defined(PRECISION_z) || defined(PRECISION_c)
CHAMELEON_CL_CB(zplghe,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zsytrf_nopiv,  starpu_matrix_get_nx(task->handles[0]), 0, 0,                                                                           (1./3.)*M* M*M);
#endif
CHAMELEON_CL_CB(zplgsy,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zplrnt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N);
CHAMELEON_CL_CB(zplssq,                                             1,                                      1, 0,                                                4);
CHAMELEON_CL_CB(zplssq2,                                            1,                                      1, 0,                                                1);
CHAMELEON_CL_CB(zpotrf,        starpu_matrix_get_nx(task->handles[0]), 0, 0,                                                                           (1./3.)*M* M*M);
CHAMELEON_CL_CB(zssssm,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), M*M*(2.*M+starpu_matrix_get_nx(task->handles[2])));
CHAMELEON_CL_CB(zsymm,         starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_ny(task->handles[2]), 0,                                           2.*M*M *N);
CHAMELEON_CL_CB(zsyr2k,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      ( 1.+2.*M*N)*M);
CHAMELEON_CL_CB(zsyrk,         starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      ( 1.+   M)*M*N);
CHAMELEON_CL_CB(ztrasm,       starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                               0.5*M*(M+1));
CHAMELEON_CL_CB(ztrmm,         starpu_matrix_get_nx(task->handles[1]), starpu_matrix_get_ny(task->handles[1]), 0,                                               M*M*N);
CHAMELEON_CL_CB(ztrsm,         starpu_matrix_get_nx(task->handles[1]), starpu_matrix_get_ny(task->handles[1]), 0,                                               M*M*N);
CHAMELEON_CL_CB(ztrtri,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (1./3.)*M *M*M);
CHAMELEON_CL_CB(ztslqt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),     2. *M* M*M);
CHAMELEON_CL_CB(ztsmlq,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (4.0*M+starpu_matrix_get_nx(task->handles[3]))*M*M);
CHAMELEON_CL_CB(ztsmqr,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (4.0*M+starpu_matrix_get_nx(task->handles[3]))*M*M);
CHAMELEON_CL_CB(ztsqrt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),     2. *M* M*M);
CHAMELEON_CL_CB(ztstrf,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),         M* M*M);
CHAMELEON_CL_CB(zttlqt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),     1. *M* M*M);
CHAMELEON_CL_CB(zttmlq,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (2.0*M+starpu_matrix_get_nx(task->handles[3]))*M*M);
CHAMELEON_CL_CB(zttmqr,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), (2.0*M+starpu_matrix_get_nx(task->handles[3]))*M*M);
CHAMELEON_CL_CB(zttqrt,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),     1. *M* M*M);
CHAMELEON_CL_CB(zunmlq,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),     2. *M* M*M);
CHAMELEON_CL_CB(zunmqr,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_nx(task->handles[0]),     2. *M* M*M);
