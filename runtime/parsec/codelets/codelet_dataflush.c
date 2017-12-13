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
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks.h"

void MORSE_TASK_flush_data( const MORSE_option_t *options,
                            const MORSE_desc_t *A, int Am, int An )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_data_flush( PARSEC_dtd_taskpool, RTBLKADDR( A, MORSE_Complex64_t, Am, An ) );
}

void MORSE_TASK_flush_desc( const MORSE_option_t *options,
                            MORSE_enum uplo, const MORSE_desc_t *A )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_data_flush_all( PARSEC_dtd_taskpool, (parsec_data_collection_t*)(A->schedopt) );

    (void)uplo;
}

void MORSE_TASK_flush_all()
{
}
