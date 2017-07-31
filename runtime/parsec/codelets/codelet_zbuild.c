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
 * @file codelet_zbuild.c
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Reazul Hoque
 * @author Guillaume Sylvand
 * @date 2016-09-05
 * @precisions normal z -> c d s
 *
 **/
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"

static inline int
CORE_zbuild_parsec(dague_execution_unit_t    *context,
                   dague_execution_context_t *this_task)
{
    MORSE_Complex64_t *A;
    int lda;
    void *user_data;
    void (*user_build_callback)(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data) ;
    int row_min, row_max, col_min, col_max;

    dague_dtd_unpack_args(
        this_task,
        UNPACK_VALUE, &row_min,
        UNPACK_VALUE, &row_max,
        UNPACK_VALUE, &col_min,
        UNPACK_VALUE, &col_max,
        UNPACK_DATA,  &A,
        UNPACK_VALUE, &lda,
        UNPACK_VALUE, &user_data,
        UNPACK_VALUE, &user_build_callback );

    user_build_callback(row_min, row_max, col_min, col_max, A, lda, user_data);

    return 0;
}

void MORSE_TASK_zbuild( const MORSE_option_t *options,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        void *user_data, void* user_build_callback )
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);
    int row_min, row_max, col_min, col_max;
    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;

    dague_insert_task(
        DAGUE_dtd_handle, CORE_zbuild_parsec, "zbuild",
        sizeof(int),   &row_min,                          VALUE,
        sizeof(int),   &row_max,                          VALUE,
        sizeof(int),   &col_min,                          VALUE,
        sizeof(int),   &col_max,                          VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     OUTPUT | REGION_FULL,
        sizeof(int),   &lda,                              VALUE,
        sizeof(void*), &user_data,                        VALUE,
        sizeof(void*), &user_build_callback,              VALUE,
        0);
}
