/**
 *
 * @copyright (c) 2009-2015 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2015 Inria. All rights reserved.
 * @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
#include "chameleon_parsec.h"
#include "chameleon_timer.h"

double RUNTIME_get_time(){
    return CHAMELEON_timer();
}

void RUNTIME_start_profiling()
{
    morse_warning("RUNTIME_start_profiling()", "FxT profiling is not available with PaRSEC\n");
}

void RUNTIME_stop_profiling()
{
    morse_warning("RUNTIME_stop_profiling()", "FxT profiling is not available with PaRSEC\n");
}

void RUNTIME_start_stats()
{
    morse_warning("RUNTIME_start_stats()", "pruning stats are not available with PaRSEC\n");
}

void RUNTIME_stop_stats()
{
    morse_warning("RUNTIME_stop_stats()", "pruning stats are not available with PaRSEC\n");
}

void RUNTIME_schedprofile_display(void)
{
    morse_warning("RUNTIME_schedprofile_display(parsec)", "Scheduler profiling is not available with PaRSEC\n");
}

void RUNTIME_kernelprofile_display(void)
{
    morse_warning("RUNTIME_kernelprofile_display(parsec)", "Kernel profiling is not available with PaRSEC\n");
}
