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
 * @file runtime_codelets.h
 *
 *  MORSE codelets kernel
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver,
 *  and INRIA Bordeaux Sud-Ouest
 *
 * @version 0.9.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 **/

#ifndef _CODELETS_H_
#define _CODELETS_H_

#include "runtime_codelet_profile.h"

//#undef STARPU_CUDA_ASYNC
#ifdef STARPU_CUDA_ASYNC
#define CODELET_CUDA_FLAGS(flags) .cuda_flags = {(flags)},
#else
#define CODELET_CUDA_FLAGS(flags)
#endif

#define CODELETS_ALL(cl_name, _nbuffers, cpu_func_name, cuda_func_name, _original_location_, cuda_flags)	\
    struct starpu_perfmodel cl_##cl_name##_fake = {                           \
        .type   = STARPU_HISTORY_BASED,                                       \
        .symbol = "fake_"#cl_name                                             \
    };                                                                        \
                                                                              \
    struct starpu_perfmodel cl_##cl_name##_model = {                          \
        .type   = STARPU_HISTORY_BASED,                                       \
        .symbol = ""#cl_name                                                  \
    };                                                                        \
                                                                              \
    struct starpu_codelet cl_##cl_name = {                                    \
        .where     = (_original_location_),                                   \
        .cpu_func  = ((cpu_func_name)),                                       \
        CODELET_CUDA_FLAGS(cuda_flags)                                        \
        .cuda_func = ((cuda_func_name)),                                      \
        .nbuffers  = ((_nbuffers)),                                           \
        .model     = &cl_##cl_name##_model,                                   \
        .name      = #cl_name                                                 \
    };                                                                        \
                                                                              \
    void cl_##cl_name##_restrict_where(uint32_t where)                        \
    {                                                                         \
      if ( cl_##cl_name.where & where )                                       \
        cl_##cl_name.where = (cl_##cl_name.where & where);                    \
    }                                                                         \
                                                                              \
    void cl_##cl_name##_restore_where(void)                                   \
    {                                                                         \
        cl_##cl_name.where = (_original_location_);                           \
    }                                                                         \
                                                                              \
    void cl_##cl_name##_restore_model(void)                                   \
    {                                                                         \
        cl_##cl_name.model = &cl_##cl_name##_model;                           \
    }

#define CODELETS_CPU(name, _nbuffers, cpu_func_name)                          \
  CODELETS_ALL( name, _nbuffers, cpu_func_name, NULL, STARPU_CPU, 0 )

#define CODELETS_GPU(name, _nbuffers, cpu_func_name, cuda_func_name, cuda_flags) \
  CODELETS_ALL( name, _nbuffers, cpu_func_name, cuda_func_name, STARPU_CPU  | STARPU_CUDA, cuda_flags )


#define CODELETS_ALL_HEADER(name)                                             \
     CHAMELEON_CL_CB_HEADER(name)                                             \
     void cl_##name##_load_fake_model(void);                                  \
     void cl_##name##_restore_model(void);                                    \
     extern struct starpu_codelet cl_##name;                                  \
     void cl_##name##_restrict_where(uint32_t where);                         \
     void cl_##name##_restore_where(void);

#if defined(CHAMELEON_USE_CUDA)
#define CODELETS(name, _nbuffers, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_GPU(name, _nbuffers, cpu_func_name, cuda_func_name, cuda_flags)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#elif defined(CHAMELEON_SIMULATION)
#define CODELETS(name, _nbuffers, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_GPU(name, _nbuffers, cpu_func_name, (starpu_cuda_func_t) 1, cuda_flags)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#else
#define CODELETS(name, _nbuffers, cpu_func_name, cuda_func_name, cuda_flags) \
    CODELETS_CPU(name, _nbuffers, cpu_func_name)

#define CODELETS_HEADER(name)  CODELETS_ALL_HEADER(name)
#endif

#define SCODELETS_HEADER(name)                CODELETS_HEADER(s##name)
#define DCODELETS_HEADER(name)                CODELETS_HEADER(d##name)
#define CCODELETS_HEADER(name)                CODELETS_HEADER(c##name)
#define ZCODELETS_HEADER(name)                CODELETS_HEADER(z##name)

#define SCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(s##name)
#define DCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(d##name)
#define CCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(c##name)
#define ZCODELETS_CPU_HEADER(name)        CODELETS_CPU_HEADER(z##name)

#endif /* _CODELETS_H_ */
