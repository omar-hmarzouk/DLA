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
 * @file pzgebrd_ge2gb.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "control/common.h"

void morse_pzgebrd_ge2gb(MORSE_desc_t A, MORSE_desc_t T, MORSE_desc_t D,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    int k;
    int tempkm, tempkn;
    if (A.m >= A.n){
       for (k = 0; k < A.nt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

           morse_pzgeqrf(
               morse_desc_submatrix(&A, k*A.mb, k*A.nb, A.m-k*A.mb, tempkn),
               morse_desc_submatrix(&T, k*T.mb, k*T.nb, T.m-k*T.mb, tempkn),
               morse_desc_submatrix(&D, k*T.mb, k*T.nb, T.m-k*T.mb, tempkn),
               sequence, request);

           morse_pzunmqr(
               MorseLeft,
               MorseConjTrans,
               morse_desc_submatrix(&A, k*A.mb,     k*A.nb, A.m-k*A.mb, tempkn),
               morse_desc_submatrix(&A, k*A.mb, (k+1)*A.nb, A.m-k*A.mb, A.n-(k+1)*A.nb),
               morse_desc_submatrix(&T, k*T.mb,     k*T.nb, T.m-k*T.mb, tempkn),
               morse_desc_submatrix(&D, k*T.mb,     k*T.nb, T.m-k*T.mb, tempkn),
               sequence, request);

           if (k+1 < A.nt){
              tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;

              morse_pzgelqf(
                  morse_desc_submatrix(&A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  morse_desc_submatrix(&T, k*T.mb, (k+1)*T.nb, T.mb,   T.n-(k+1)*T.nb),
                  morse_desc_submatrix(&D, k*T.mb, (k+1)*T.nb, T.mb,   T.n-(k+1)*T.nb),
                  sequence, request);

              morse_pzunmlq(
                  MorseRight, MorseConjTrans,
                  morse_desc_submatrix(&A,     k*A.mb, (k+1)*A.nb, tempkm,         A.n-(k+1)*A.nb),
                  morse_desc_submatrix(&A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  morse_desc_submatrix(&T,     k*T.mb, (k+1)*T.nb, T.mb,           T.n-(k+1)*T.nb),
                  morse_desc_submatrix(&D,     k*T.mb, (k+1)*T.nb, T.mb,           T.n-(k+1)*T.nb),
                  sequence, request);
           }
       }
    }
    else{
       for (k = 0; k < A.mt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

           morse_pzgelqf(
               morse_desc_submatrix(&A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               morse_desc_submatrix(&T, k*T.mb, k*T.nb, T.mb,   T.n-k*T.nb),
               morse_desc_submatrix(&D, k*T.mb, k*T.nb, T.mb,   T.n-k*T.nb),
               sequence, request);

           morse_pzunmlq(
               MorseRight, MorseConjTrans,
               morse_desc_submatrix(&A,     k*A.mb, k*A.nb, tempkm,         A.n-k*A.nb),
               morse_desc_submatrix(&A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, A.n-k*A.nb),
               morse_desc_submatrix(&T,     k*T.mb, k*T.nb, T.mb,           T.n-k*T.nb),
               morse_desc_submatrix(&D,     k*T.mb, k*T.nb, T.mb,           T.n-k*T.nb),
               sequence, request);

           if (k+1 < A.mt){
              tempkm = k+1 == A.mt-1 ? A.m-(k+1)*A.mb : A.mb;

              morse_pzgeqrf(
                   morse_desc_submatrix(&A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, tempkn),
                   morse_desc_submatrix(&T, (k+1)*T.mb, k*T.nb, T.m-(k+1)*T.mb, tempkn),
                   morse_desc_submatrix(&D, (k+1)*T.mb, k*T.nb, T.m-(k+1)*T.mb, tempkn),
                   sequence, request);

              morse_pzunmqr(
                  MorseLeft, MorseConjTrans,
                  morse_desc_submatrix(&A, (k+1)*A.mb,     k*A.nb, A.m-(k+1)*A.mb, tempkn),
                  morse_desc_submatrix(&A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  morse_desc_submatrix(&T, (k+1)*T.mb,     k*T.nb, T.m-(k+1)*T.mb, tempkn),
                  morse_desc_submatrix(&D, (k+1)*T.mb,     k*T.nb, T.m-(k+1)*T.mb, tempkn),
                  sequence, request);
           }
       }
    }
}
