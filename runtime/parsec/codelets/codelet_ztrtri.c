/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 * @file codelet_ztrtri.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

static int
CORE_ztrtri_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    MORSE_enum *uplo;
    MORSE_enum *diag;
    int *N;
    MORSE_Complex64_t *A;
    int *LDA;
    int *iinfo;
    int info;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &uplo,
                          UNPACK_VALUE, &diag,
                          UNPACK_VALUE, &N,
                          UNPACK_DATA,  &A,
                          UNPACK_VALUE, &LDA,
                          UNPACK_VALUE, &iinfo
                        );


    CORE_ztrtri(*uplo, *diag, *N, A, *LDA, &info);

    return 0;
}

void MORSE_TASK_ztrtri(const MORSE_option_t *options,
                       MORSE_enum uplo, MORSE_enum diag,
                       int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,          CORE_ztrtri_parsec,    "trtri",
                            sizeof(MORSE_enum),         &uplo,                  VALUE,
                            sizeof(MORSE_enum),         &diag,                  VALUE,
                            sizeof(int),                &n,                     VALUE,
                            PASSED_BY_REF,              RTBLKADDR( A, MORSE_Complex64_t, Am, An ),   INOUT | REGION_FULL,
                            sizeof(int),                &lda,                   VALUE,
                            sizeof(int),                &iinfo,                 VALUE,
                            0);
}
