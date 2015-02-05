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
 * @file codelet_zgetrf.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "coreblas/include/lapacke.h"
#include "runtime/quark/include/morse_quark.h"

void MORSE_TASK_zgetrf(MORSE_option_t *options,
                       int m, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       int *IPIV,
                       MORSE_bool check_info, int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_GETRF;
    QUARK_Insert_Task(opt->quark, CORE_zgetrf_quark, (Quark_Task_Flags*)opt,
        sizeof(int),                        &m,             VALUE,
        sizeof(int),                        &n,             VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                     INOUT | LOCALITY,
        sizeof(int),                        &lda,           VALUE,
        sizeof(int)*nb,                      IPIV,                  OUTPUT,
        sizeof(MORSE_sequence_t*),           &(options->sequence),      VALUE,
        sizeof(MORSE_request_t*),            &(options->request),       VALUE,
        sizeof(MORSE_bool),                &check_info,    VALUE,
        sizeof(int),                        &iinfo,         VALUE,
        0);
}


void CORE_zgetrf_quark(Quark *quark)
{
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int *IPIV;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    MORSE_bool check_info;
    int iinfo;
    int info;

    quark_unpack_args_9(quark, m, n, A, lda, IPIV, sequence, request, check_info, iinfo);
    info = LAPACKE_zgetrf_work(LAPACK_COL_MAJOR, m, n, A, lda, IPIV );
    if (info != MORSE_SUCCESS && check_info)
        RUNTIME_sequence_flush(quark, sequence, request, iinfo+info);
}
