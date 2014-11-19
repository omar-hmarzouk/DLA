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
 * @file global.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @date 2010-11-15
 *
 **/

/*******************************************************************************
 *  LAPACK Constants
 **/
char *morse_lapack_constants[] =
{
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 100
    "Row",                  // 101: MorseRowMajor
    "Column",               // 102: MorseColMajor
    "", "", "", "", "", "", "", "",
    "No transpose",         // 111: MorseNoTrans
    "Transpose",            // 112: MorseTrans
    "Conjugate transpose",  // 113: MorseConjTrans
    "", "", "", "", "", "", "",
    "Upper",                // 121: MorseUpper
    "Lower",                // 122: MorseLower
    "All",                  // 123: MorseUpperLower
    "", "", "", "", "", "", "",
    "Non-unit",             // 131: MorseNonUnit
    "Unit",                 // 132: MorseUnit
    "", "", "", "", "", "", "", "",
    "Left",                 // 141: MorseLeft
    "Right",                // 142: MorseRight
    "", "", "", "", "", "", "", "",
    "",                     // 151:
    "",                     // 152:
    "",                     // 153:
    "",                     // 154:
    "",                     // 155:
    "",                     // 156:
    "Epsilon",              // 157: MorseEps
    "",                     // 158:
    "",                     // 159:
    "",                     // 160:
    "", "", "", "", "", "", "", "", "", "",
    "One norm",             // 171: MorseOneNorm
    "",                     // 172: MorseRealOneNorm
    "",                     // 173: MorseTwoNorm
    "Frobenius norm",       // 174: MorseFrobeniusNorm
    "Infinity norm",        // 175: MorseInfNorm
    "",                     // 176: MorseRealInfNorm
    "Maximum norm",         // 177: MorseMaxNorm
    "",                     // 178: MorseRealMaxNorm
    "",                     // 179
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 200
    "Uniform",              // 201: MorseDistUniform
    "Symmetric",            // 202: MorseDistSymmetric
    "Normal",               // 203: MorseDistNormal
    "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 240
    "Hermitian",            // 241 MorseHermGeev
    "Positive ev Hermitian",// 242 MorseHermPoev
    "NonSymmetric pos sv",  // 243 MorseNonsymPosv
    "Symmetric pos sv",     // 244 MorseSymPosv
    "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 290
    "No Packing",           // 291 MorseNoPacking
    "U zero out subdiag",   // 292 MorsePackSubdiag
    "L zero out superdiag", // 293 MorsePackSupdiag
    "C",                    // 294 MorsePackColumn
    "R",                    // 295 MorsePackRow
    "B",                    // 296 MorsePackLowerBand
    "Q",                    // 297 MorsePackUpeprBand
    "Z",                    // 298 MorsePackAll
    "",                     // 299

    "",                     // 300
    "No vectors",           // 301 MorseNoVec
    "Vectors needed",       // 302 MorseVec
    "I",                    // 303 MorseIvec
    "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "",                     // 390
    "Forward",              // 391
    "Backward",             // 392
    "", "", "", "", "", "", "", "",
    "Columnwise",           // 401
    "Rowwise",              // 402
    "", "", "", "", "", "", "", ""  // Remember to add a coma!
};
