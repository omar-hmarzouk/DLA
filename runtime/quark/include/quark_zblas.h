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
 * @file quark_zblas.h
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
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#ifndef _QUARK_ZBLAS_H_
#define _QUARK_ZBLAS_H_

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif



/** ****************************************************************************
 *  Declarations of QUARK wrappers (called by QUARK) - alphabetical order
 **/
void CORE_dzasum_quark(Quark *quark);
void CORE_dzasum_f1_quark(Quark *quark);
void CORE_zaxpy_quark(Quark *quark);
void CORE_zgeadd_quark(Quark *quark);
void CORE_zbrdalg_quark(Quark *quark);
void CORE_zgelqt_quark(Quark *quark);
void CORE_zgemm_quark(Quark *quark);
void CORE_zgeqrt_quark(Quark *quark);
void CORE_zgessm_quark(Quark *quark);
void CORE_zgessq_quark(Quark *quark);
void CORE_zgetrf_quark(Quark *quark);
void CORE_zgetrf_incpiv_quark(Quark *quark);
void CORE_zgetrf_nopiv_quark(Quark *quark);
void CORE_zgetrf_reclap_quark(Quark *quark);
void CORE_zgetrf_rectil_quark(Quark* quark);
void CORE_zgetrip_quark(Quark *quark);
void CORE_zgetrip_f1_quark(Quark *quark);
void CORE_zgetrip_f2_quark(Quark *quark);
#ifdef COMPLEX
void CORE_zhemm_quark(Quark *quark);
void CORE_zherk_quark(Quark *quark);
void CORE_zher2k_quark(Quark *quark);
#endif
void CORE_zhegst_quark(Quark *quark);
void CORE_zherfb_quark(Quark *quark);
void CORE_zhessq_quark(Quark *quark);
void CORE_zlacpy_quark(Quark *quark);
void CORE_zlatro_quark(Quark *quark);
void CORE_zlange_quark(Quark *quark);
void CORE_zlange_max_quark(Quark *quark);
#ifdef COMPLEX
void CORE_zlanhe_quark(Quark *quark);
#endif
void CORE_zlansy_quark(Quark *quark);
void CORE_zlantr_quark(Quark *quark);
void CORE_zlaset_quark(Quark *quark);
void CORE_zlaset2_quark(Quark *quark);
void CORE_zlatro_quark(Quark *quark);
void CORE_zlauum_quark(Quark *quark);
void CORE_zpamm_quark(Quark *quark);
void CORE_zplghe_quark(Quark *quark);
void CORE_zplgsy_quark(Quark *quark);
void CORE_zplrnt_quark(Quark *quark);
void CORE_zplssq_quark(Quark *quark);
void CORE_zplssq2_quark(Quark *quark);
void CORE_zpotrf_quark(Quark *quark);
void CORE_zshift_quark(Quark *quark);
void CORE_zshiftw_quark(Quark *quark);
void CORE_zssssm_quark(Quark *quark);
void CORE_zsymm_quark(Quark *quark);
void CORE_zsyrk_quark(Quark *quark);
void CORE_zsyr2k_quark(Quark *quark);
void CORE_zsyssq_quark(Quark *quark);
void CORE_zsytrf_nopiv_quark(Quark *quark);
void CORE_zswpab_quark(Quark *quark);
void CORE_zswptr_ontile_quark(Quark *quark);
void CORE_ztrasm_quark(Quark *quark);
void CORE_ztrdalg_quark(Quark *quark);
void CORE_ztrmm_quark(Quark *quark);
void CORE_ztrsm_quark(Quark *quark);
void CORE_ztrssq_quark(Quark *quark);
void CORE_ztrtri_quark(Quark *quark);
void CORE_ztslqt_quark(Quark *quark);
void CORE_ztsmlq_quark(Quark *quark);
void CORE_ztsmlq_hetra1_quark(Quark *quark);
void CORE_ztsmlq_corner_quark(Quark *quark);
void CORE_ztsmqr_quark(Quark *quark);
void CORE_ztsmqr_hetra1_quark(Quark *quark);
void CORE_ztsmqr_corner_quark(Quark *quark);
void CORE_ztsqrt_quark(Quark *quark);
void CORE_ztstrf_quark(Quark *quark);
void CORE_zttmqr_quark(Quark *quark);
void CORE_zttqrt_quark(Quark *quark);
void CORE_zttmlq_quark(Quark *quark);
void CORE_zttlqt_quark(Quark *quark);
void CORE_zunmlq_quark(Quark *quark);
void CORE_zunmqr_quark(Quark *quark);

void CORE_zlaswp_quark(Quark* quark);
void CORE_zlaswp_f2_quark(Quark* quark);
void CORE_zlaswp_ontile_quark(Quark *quark);
void CORE_zlaswp_ontile_f2_quark(Quark *quark);
void CORE_zlaswpc_ontile_quark(Quark *quark);
void CORE_ztrmm_p2_quark(Quark* quark);
void CORE_zgemm_f2_quark(Quark* quark);
void CORE_zgemm_p2_quark(Quark* quark);
void CORE_zgemm_p2f1_quark(Quark* quark);
void CORE_zgemm_p3_quark(Quark* quark);

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
