/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 * @file codelet_zlag2c.c
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

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
static int
CORE_zlag2c_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *m;
    int *n;
    dague_data_copy_t *gA;
    int *lda;
    dague_data_copy_t *gB;
    int *ldb;
    int info;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &m,
                          UNPACK_VALUE, &n,
                          UNPACK_DATA,  &gA,
                          UNPACK_VALUE, &lda,
                          UNPACK_DATA,  &gB,
                          UNPACK_VALUE, &ldb,
                          );

    void *A = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gA);
    void *B = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gB);

    CORE_zlag2c( *m, *n, A, *lda, B, *ldb);

    return 0;
}

void MORSE_TASK_zlag2c(MORSE_option_t *options,
                       int m, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_zlag2c_parsec,               "lag2c",
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT | REGION_FULL,
        sizeof(int),                        &lda,       VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ),     OUTPUT | REGION_FULL,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 **/
static int
CORE_clag2z_parsec(dague_execution_unit_t *context, dague_execution_context_t *this_task)
{
    int *m;
    int *n;
    dague_data_copy_t *gA;
    int *lda;
    dague_data_copy_t *gB;
    int *ldb;

    dague_dtd_unpack_args(this_task,
                          UNPACK_VALUE, &m,
                          UNPACK_VALUE, &n,
                          UNPACK_DATA,  &gA,
                          UNPACK_VALUE, &lda,
                          UNPACK_DATA,  &gB,
                          UNPACK_VALUE, &ldb,
                          );

    void *A = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gA);
    void *B = DAGUE_DATA_COPY_GET_PTR((dague_data_copy_t *)gB);

    CORE_clag2z( *m, *n, A, *lda, B, *ldb );

    return 0;
}

void MORSE_TASK_clag2z(MORSE_option_t *options,
                       int m, int n, int nb,
                       MORSE_desc_t *A, int Am, int An, int lda,
                       MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    dague_dtd_handle_t* DAGUE_dtd_handle = (dague_dtd_handle_t *)(options->sequence->schedopt);

    insert_task_generic_fptr(DAGUE_dtd_handle,      CORE_clag2z_parsec,               "lag2z",
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ),     INPUT | REGION_FULL,
        sizeof(int),                        &lda,       VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ),     INOUT | REGION_FULL,
        sizeof(int),                        &ldb,       VALUE,
        0);
}
