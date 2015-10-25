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
 * @file pzpotrimm.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void morse_pzpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C,
                     MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int lda, ldb, ldc;
    int ldak, ldam, ldan;
    int tempkm, tempmm, tempnn, tempkn;

    MORSE_Complex64_t alpha  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t beta  = (MORSE_Complex64_t) 0.0;
    MORSE_Complex64_t zbeta;
    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;


    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

#ifdef CHAMELEON_USE_MAGMA
    {
      int nb = magma_get_zpotrf_nb(A->nb);
      RUNTIME_options_ws_alloc( &options, nb*nb, 0 );
    }
#endif

    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);

            MORSE_TASK_zpotrf(
                &options,
                MorseLower, tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, MorseLower, MorseConjTrans, MorseNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                ldan = BLKLDD(A, n);
                MORSE_TASK_zherk(
                    &options,
                    MorseLower, MorseNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k), ldan,
                     1.0, A(n, n), ldan);

                for (m = n+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k), ldam,
                               A(n, k), ldan,
                        zone,  A(m, n), ldam);
                }
                MORSE_TASK_dataflush( &options, A(n, k) );
            }
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->nt; k++) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, MorseNonUnit,
                    tempmm, tempkn, A->mb,
                    mzone, A(k, k), ldak,
                           A(m, k), ldam);
            }
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                ldam = BLKLDD(A, m);
                for (n = 0; n < k; n++) {
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        tempmm, A->nb, tempkn, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                MORSE_TASK_dataflush( &options, A(m, k) );
            }
            for (n = 0; n < k; n++) {
                MORSE_TASK_dataflush( &options, A(k, n) );
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, MorseNonUnit,
                    tempkn, A->nb, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );
            MORSE_TASK_ztrtri(
                &options,
                uplo, MorseNonUnit,
                tempkn, A->mb,
                A(k, k), ldak, A->nb*k);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for(n = 0; n < k; n++) {
                ldan = BLKLDD(A, n);
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseConjTrans,
                    A->mb, tempkm, A->mb,
                    1.0, A(k, n), ldak,
                    1.0, A(n, n), ldan);

                for(m = n+1; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        A->mb, A->nb, tempkm, A->mb,
                        1.0, A(k, m), ldak,
                             A(k, n), ldak,
                        1.0, A(m, n), ldam);
                }
            }
            for (n = 0; n < k; n++) {
                MORSE_TASK_dataflush( &options, A(k, n) );
                MORSE_TASK_ztrmm(
                    &options,
                    MorseLeft, uplo, MorseConjTrans, MorseNonUnit,
                    tempkm, A->nb, A->mb,
                    1.0, A(k, k), ldak,
                         A(k, n), ldak);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );
            MORSE_TASK_zlauum(
                &options,
                uplo, tempkm, A->mb,
                A(k, k), ldak);
        }
        /*
         *  ZSYMM
         */
        for (m = 0; m < C->mt; m++) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
            ldc = BLKLDD(C, m);
            for (n = 0; n < C->nt; n++) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                lda = BLKLDD(A, n);
                ldb = BLKLDD(B, m);
                for (k = 0; k < C->nt; k++) {
                    tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
                    ldak = BLKLDD(A, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < n) {
                       MORSE_TASK_zgemm(
                           &options,
                           MorseNoTrans, MorseTrans,
                           tempmm, tempnn, tempkn, A->mb,
                           alpha, B(m, k), ldb,  /* ldb * K */
                           A(n, k), lda,  /* lda * K */
                           zbeta, C(m, n), ldc); /* ldc * Y */
                    }
                    else {
                        if (k == n) {
                           MORSE_TASK_zsymm(
                               &options,
                               MorseRight, uplo,
                               tempmm, tempnn, A->mb,
                               alpha, A(k, k), ldak, /* ldak * Y */
                               B(m, k), ldb,  /* ldb  * Y */
                               zbeta, C(m, n), ldc); /* ldc  * Y */
                        }
                        else {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldb,  /* ldb  * K */
                                A(k, n), ldak, /* ldak * Y */
                                zbeta, C(m, n), ldc); /* ldc  * Y */
                        }
                    }
                }
            }
        }
    }
    /*
     *  MorseUpper
     */
    else {
        /*
         *  ZPOTRF
         */
        for (k = 0; k < A->nt; k++) {
            tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);
            MORSE_TASK_zpotrf(
                &options,
                MorseUpper,
                tempkm, A->mb,
                A(k, k), ldak, A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, MorseUpper, MorseConjTrans, MorseNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k), ldak,
                          A(k, n), ldak);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
                ldam = BLKLDD(A, m);

                MORSE_TASK_zherk(
                    &options,
                    MorseUpper, MorseConjTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam);

                for (n = m+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                    MORSE_TASK_zgemm(
                        &options,
                        MorseConjTrans, MorseNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m), ldak,
                               A(k, n), ldak,
                        zone,  A(m, n), ldam);
                }
                MORSE_TASK_dataflush( &options, A(k, m) );
            }
        }
        /*
         *  ZTRTRI
         */
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            ldak = BLKLDD(A, k);
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                MORSE_TASK_ztrsm(
                    &options,
                    MorseLeft, uplo, MorseNoTrans, MorseNonUnit,
                    tempkm, tempnn, A->mb,
                    mzone, A(k, k), ldak,
                           A(k, n), ldak);
            }
            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                for (m = 0; m < k; m++) {
                    ldam = BLKLDD(A, m);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseNoTrans,
                        A->mb, tempnn, tempkm, A->mb,
                        zone, A(m, k), ldam,
                              A(k, n), ldak,
                        zone, A(m, n), ldam);
                }
                MORSE_TASK_dataflush( &options, A(k, n) );
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                MORSE_TASK_dataflush( &options, A(m, k) );
                MORSE_TASK_ztrsm(
                    &options,
                    MorseRight, uplo, MorseNoTrans, MorseNonUnit,
                    A->mb, tempkm, A->mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );
            MORSE_TASK_ztrtri(
                &options,
                uplo, MorseNonUnit,
                tempkm, A->mb,
                A(k, k), ldak, A->mb*k);
        }
        /*
         *  ZLAUUM
         */
        for (k = 0; k < A->mt; k++) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
            ldak = BLKLDD(A, k);

            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                MORSE_TASK_zherk(
                    &options,
                    uplo, MorseNoTrans,
                    A->mb, tempkn, A->mb,
                    1.0, A(m, k), ldam,
                    1.0, A(m, m), ldam);

                for (n = m+1; n < k; n++){
                    ldan = BLKLDD(A, n);
                    MORSE_TASK_zgemm(
                        &options,
                        MorseNoTrans, MorseConjTrans,
                        A->mb, A->nb, tempkn, A->mb,
                        1.0, A(m, k), ldam,
                             A(n, k), ldan,
                        1.0, A(m, n), ldam);
                }
            }
            for (m = 0; m < k; m++) {
                ldam = BLKLDD(A, m);
                MORSE_TASK_dataflush( &options, A(m, k) );
                MORSE_TASK_ztrmm(
                    &options,
                    MorseRight, uplo, MorseConjTrans, MorseNonUnit,
                    A->mb, tempkn, A->mb,
                    1.0, A(k, k), ldak,
                         A(m, k), ldam);
            }
            MORSE_TASK_dataflush( &options, A(k, k) );
            MORSE_TASK_zlauum(
                &options,
                uplo, tempkn, A->mb,
                A(k, k), ldak);
        }
        /*
         *  ZSYMM
         */
        for (m = 0; m < C->mt; m++) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
            ldc = BLKLDD(C, m);
            for (n = 0; n < C->nt; n++) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                lda = BLKLDD(A, n);
                ldb = BLKLDD(B, m);
                for (k = 0; k < C->nt; k++) {
                    tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
                    ldak = BLKLDD(A, k);
                    zbeta = k == 0 ? beta : zone;
                    if (k < n) {
                        MORSE_TASK_zgemm(
                            &options,
                            MorseNoTrans, MorseNoTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, B(m, k), ldb,  /* ldb  * K */
                            A(k, n), ldak, /* ldak * Y */
                            zbeta, C(m, n), ldc); /* ldc  * Y */
                    }
                    else {
                        if (k == n) {
                            MORSE_TASK_zsymm(
                                &options,
                                MorseRight, uplo,
                                tempmm, tempnn, A->mb,
                                alpha, A(k, k), ldak, /* ldak * Y */
                                B(m, k), ldb,  /* ldb  * Y */
                                zbeta, C(m, n), ldc); /* ldc  * Y */
                        }
                        else {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseTrans,
                                tempmm, tempnn, tempkn, A->mb,
                                alpha, B(m, k), ldb,  /* ldb * K */
                                A(n, k), lda,  /* lda * K */
                                zbeta, C(m, n), ldc); /* ldc * Y */
                        }
                    }
                }
            }
        }
    }

#ifdef CHAMELEON_USE_MAGMA
    RUNTIME_options_ws_free(&options);
#endif
    RUNTIME_options_finalize(&options, morse);
}
