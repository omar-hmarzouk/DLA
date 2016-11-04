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
 * @file codelet_zplghe.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @precisions normal z -> c
 *
 **/
#include "runtime/parsec/include/morse_parsec.h"

static int
CORE_zplghe_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    double *bump;
    int *m;
    int *n;
    MORSE_Complex64_t *A;
    int *lda;
    int *bigM;
    int *m0;
    int *n0;
    unsigned long long int *seed;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &bump,
                          UNPACK_VALUE, &m,
                          UNPACK_VALUE, &n,
                          UNPACK_DATA,  &A,
                          UNPACK_VALUE, &lda,
                          UNPACK_VALUE, &bigM,
                          UNPACK_VALUE, &m0,
                          UNPACK_VALUE, &n0,
                          UNPACK_VALUE, &seed
                        );


    CORE_zplghe( *bump, *m, *n, A, *lda, *bigM, *m0, *n0, *seed );

    return 0;
}

void MORSE_TASK_zplghe( const MORSE_option_t *options,
                        double bump, int m, int n, MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,  CORE_zplghe_parsec,            "zplghe",
                             sizeof(double),    &bump,                          VALUE,
                             sizeof(int),       &m,                             VALUE,
                             sizeof(int),       &n,                             VALUE,
                             PASSED_BY_REF,     RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     OUTPUT | REGION_FULL,
                             sizeof(int),       &lda,                           VALUE,
                             sizeof(int),       &bigM,                          VALUE,
                             sizeof(int),       &m0,                            VALUE,
                             sizeof(int),       &n0,                            VALUE,
                             sizeof(unsigned long long int),       &seed,       VALUE,
                             0);
}
