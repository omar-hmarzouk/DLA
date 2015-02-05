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
 * @file pzlauum.c
 *
 *  MORSE auxiliary routines
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
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
/***************************************************************************//**
 *  Parallel UU' or L'L operation - dynamic scheduling
 **/
void morse_pzlauum(MORSE_enum uplo, MORSE_desc_t *A,
                          MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int ldam;
    int tempkm, tempmm, tempnn;

    MORSE_Complex64_t zone = (MORSE_Complex64_t)1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);
    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (m = 0; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            for(n = 0; n < m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseConjTrans,
                    tempnn, tempmm, A->mb,
                    1.0, A(m, n), ldam,
                    1.0, A(n, n), A->mb);

                for(k = n+1; k < m; k++) {
                    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        tempkm, tempnn, tempmm, A->mb,
                        zone, A(m, k), ldam,
                              A(m, n), ldam,
                        zone, A(k, n), A->mb);
                }
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztrmm(
                    &options,
                    MorseLeft, uplo, MorseConjTrans, MorseNonUnit,
                    tempmm, tempnn, A->mb,
                    zone, A(m, m), ldam,
                          A(m, n), ldam);
            }
            MORSE_TASK_zlauum(
                &options,
                uplo,
                tempmm,
                A->mb, A(m, m), ldam);
        }
    }
    /*
     *  MorseUpper
     */
    else {
        for (m = 0; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            for (n = 0; n < m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseNoTrans,
                    tempnn, tempmm, A->mb,
                    1.0, A(n, m), A->mb,
                    1.0, A(n, n), A->mb);

                for (k = n+1; k < m; k++){
                    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        tempnn, tempkm, tempmm, A->mb,
                        zone, A(n, m), A->mb,
                              A(k, m), A->mb,
                        zone, A(n, k), A->mb);
                }
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztrmm(
                    &options,
                    MorseRight, uplo, MorseConjTrans, MorseNonUnit,
                    tempnn, tempmm, A->mb,
                    zone, A(m, m), ldam,
                          A(n, m), A->mb);
            }
            MORSE_TASK_zlauum(
                &options,
                uplo,
                tempmm,
                A->mb, A(m, m), ldam);
        }
    }
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
}
