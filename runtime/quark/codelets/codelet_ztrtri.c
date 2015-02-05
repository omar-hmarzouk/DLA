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
 * @file codelet_ztrtri.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "coreblas/include/lapacke.h"
#include "runtime/quark/include/morse_quark.h"

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/

void MORSE_TASK_ztrtri(MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum diag,
                       int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    QUARK_Insert_Task(
        opt->quark, CORE_ztrtri_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),                &uplo,      VALUE,
        sizeof(MORSE_enum),                &diag,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(MORSE_Complex64_t)*nb*nb,    RTBLKADDR(A, MORSE_Complex64_t, Am, An),                 INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(MORSE_sequence_t*),           &(options->sequence),  VALUE,
        sizeof(MORSE_request_t*),            &(options->request),   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}


void CORE_ztrtri_quark(Quark *quark)
{
    MORSE_enum uplo;
    MORSE_enum diag;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    MORSE_sequence_t *sequence;
    MORSE_request_t *request;
    int iinfo;

    int info;

    quark_unpack_args_8(quark, uplo, diag, N, A, LDA, sequence, request, iinfo);
    info = LAPACKE_ztrtri_work(
        LAPACK_COL_MAJOR,
        morse_lapack_const(uplo), morse_lapack_const(diag),
        N, A, LDA);
    if ((sequence->status == MORSE_SUCCESS) && (info > 0))
        RUNTIME_sequence_flush(quark, sequence, request, iinfo + info);
}
