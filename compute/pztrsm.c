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
 * @file pztrsm.c
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
/***************************************************************************//**
 *  Parallel tile triangular solve - dynamic scheduling
 **/
void morse_pztrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum trans, MORSE_enum diag,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A, MORSE_desc_t *B,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n;
    int lda, ldan, ldb;
    int tempkm, tempkn, tempmm, tempnn;

    MORSE_Complex64_t zone       = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone      = (MORSE_Complex64_t)-1.0;
    MORSE_Complex64_t minvalpha  = (MORSE_Complex64_t)-1.0 / alpha;
    MORSE_Complex64_t lalpha;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);
    /*
     *  MorseLeft / MorseUpper / MorseNoTrans
     */
    if (side == MorseLeft) {
        if (uplo == MorseUpper) {
            if (trans == MorseNoTrans) {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == 0 ? B->m-(B->mt-1)*B->mb : B->mb;
                    lda = BLKLDD(A, B->mt-1-k);
                    ldb = BLKLDD(B, B->mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(B->mt-1-k, B->mt-1-k), lda,  /* lda * tempkm */
                                    B(B->mt-1-k,        n), ldb); /* ldb * tempnn */
                    }
                    MORSE_TASK_dataflush( &options, A(B->mt-1-k, B->mt-1-k) );

                    for (m = k+1; m < B->mt; m++) {
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                B->mb, tempnn, tempkm, A->mb,
                                mzone,  A(B->mt-1-m, B->mt-1-k), A->mb,
                                        B(B->mt-1-k, n       ), ldb,
                                lalpha, B(B->mt-1-m, n       ), B->mb);
                        }
                        MORSE_TASK_dataflush( &options, A(B->mt-1-m, B->mt-1-k) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, B(B->mt-1-k, n) );
                    }
                }
            }
            /*
             *  MorseLeft / MorseUpper / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                    lda = BLKLDD(A, k);
                    ldb = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(k, k), lda,
                                    B(k, n), ldb);
                    }
                    MORSE_TASK_dataflush( &options, A(k, k) );

                    for (m = k+1; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldb = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                mzone,  A(k, m), A->mb,
                                        B(k, n), B->mb,
                                lalpha, B(m, n), ldb);
                        }
                        MORSE_TASK_dataflush( &options, A(k, m) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, B(k, n) );
                    }
                }
            }
        }
        /*
         *  MorseLeft / MorseLower / MorseNoTrans
         */
        else {
            if (trans == MorseNoTrans) {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == B->mt-1 ? B->m-k*B->mb : B->mb;
                    lda = BLKLDD(A, k);
                    ldb = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(k, k), lda,
                                    B(k, n), ldb);
                    }
                    MORSE_TASK_dataflush( &options, A(k, k) );

                    for (m = k+1; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        lda = BLKLDD(A, m);
                        ldb = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                mzone,  A(m, k), lda,
                                        B(k, n), B->mb,
                                lalpha, B(m, n), ldb);
                        }
                        MORSE_TASK_dataflush( &options, A(m, k) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, B(k, n) );
                    }
                }
            }
            /*
             *  MorseLeft / MorseLower / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == 0 ? B->m-(B->mt-1)*B->mb : B->mb;
                    lda = BLKLDD(A, B->mt-1-k);
                    ldb = BLKLDD(B, B->mt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempkm, tempnn, A->mb,
                            lalpha, A(B->mt-1-k, B->mt-1-k), lda,
                                    B(B->mt-1-k,        n), ldb);
                    }
                    MORSE_TASK_dataflush( &options, A(B->mt-1-k, B->mt-1-k) );

                    for (m = k+1; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                trans, MorseNoTrans,
                                B->mb, tempnn, tempkm, A->mb,
                                mzone,  A(B->mt-1-k, B->mt-1-m), lda,
                                        B(B->mt-1-k, n       ), ldb,
                                lalpha, B(B->mt-1-m, n       ), B->mb);
                        }
                        MORSE_TASK_dataflush( &options, A(B->mt-1-k, B->mt-1-m) );
                    }
                    for (n = 0; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, B(B->mt-1-k, n) );
                    }
                }
            }
        }
    }
    /*
     *  MorseRight / MorseUpper / MorseNoTrans
     */
    else {
        if (uplo == MorseUpper) {
            if (trans == MorseNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;
                    lda = BLKLDD(A, k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldb = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            lalpha, A(k, k), lda,  /* lda * tempkn */
                                    B(m, k), ldb); /* ldb * tempkn */
                    }
                    MORSE_TASK_dataflush( &options, A(k, k) );

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldb = BLKLDD(B, m);
                        for (n = k+1; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, tempnn, B->mb, A->mb,
                                mzone,  B(m, k), ldb,  /* ldb * B->mb   */
                                        A(k, n), lda,  /* lda * tempnn */
                                lalpha, B(m, n), ldb); /* ldb * tempnn */
                        }
                        MORSE_TASK_dataflush( &options, B(m, k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, A(k, n) );
                    }
                }
            }
            /*
             *  MorseRight / MorseUpper / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n-(B->nt-1)*B->nb : B->nb;
                    lda = BLKLDD(A, B->nt-1-k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldb = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            alpha, A(B->nt-1-k, B->nt-1-k), lda,  /* lda * tempkn */
                                   B(       m, B->nt-1-k), ldb); /* ldb * tempkn */
                        MORSE_TASK_dataflush( &options, A(B->nt-1-k, B->nt-1-k) );

                        for (n = k+1; n < B->nt; n++) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, B->nb, tempkn, A->mb,
                                minvalpha, B(m,        B->nt-1-k), ldb,  /* ldb  * tempkn */
                                           A(B->nt-1-n, B->nt-1-k), A->mb, /* A->mb * tempkn (Never last row) */
                                zone,      B(m,        B->nt-1-n), ldb); /* ldb  * B->nb   */
                        }
                        MORSE_TASK_dataflush( &options, B(m,        B->nt-1-k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, A(B->nt-1-n, B->nt-1-k) );
                    }
                }
            }
        }
        /*
         *  MorseRight / MorseLower / MorseNoTrans
         */
        else {
            if (trans == MorseNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n-(B->nt-1)*B->nb : B->nb;
                    lda = BLKLDD(A, B->nt-1-k);
                    lalpha = k == 0 ? alpha : zone;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldb = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            lalpha, A(B->nt-1-k, B->nt-1-k), lda,  /* lda * tempkn */
                                    B(       m, B->nt-1-k), ldb); /* ldb * tempkn */
                        MORSE_TASK_dataflush( &options, A(B->nt-1-k, B->nt-1-k) );

                        for (n = k+1; n < B->nt; n++) {
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, MorseNoTrans,
                                tempmm, B->nb, tempkn, A->mb,
                                mzone,  B(m,        B->nt-1-k), ldb,  /* ldb * tempkn */
                                        A(B->nt-1-k, B->nt-1-n), lda,  /* lda * B->nb   */
                                lalpha, B(m,        B->nt-1-n), ldb); /* ldb * B->nb   */
                        }
                        MORSE_TASK_dataflush( &options, B(m,        B->nt-1-k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, A(B->nt-1-k, B->nt-1-n) );
                    }
                }
            }
            /*
             *  MorseRight / MorseLower / Morse[Conj]Trans
             */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt-1 ? B->n-k*B->nb : B->nb;
                    lda = BLKLDD(A, k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
                        ldb = BLKLDD(B, m);
                        MORSE_TASK_ztrsm(
                            &options,
                            side, uplo, trans, diag,
                            tempmm, tempkn, A->mb,
                            alpha, A(k, k), lda,  /* lda * tempkn */
                                   B(m, k), ldb); /* ldb * tempkn */
                        MORSE_TASK_dataflush( &options, A(k, k) );

                        for (n = k+1; n < B->nt; n++) {
                            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                            ldan = BLKLDD(A, n);
                            MORSE_TASK_zgemm(
                                &options,
                                MorseNoTrans, trans,
                                tempmm, tempnn, B->mb, A->mb,
                                minvalpha, B(m, k), ldb,  /* ldb  * tempkn */
                                           A(n, k), ldan, /* ldan * tempkn */
                                zone,      B(m, n), ldb); /* ldb  * tempnn */
                        }
                        MORSE_TASK_dataflush( &options, B(m, k) );
                    }
                    for (n = k+1; n < B->nt; n++) {
                        MORSE_TASK_dataflush( &options, A(n, k) );
                    }
                }
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
}
