#
# Check BLAS/Lapack subroutines
#

set(TEST_CMD_shm    )
set(TEST_CMD_shmgpu --gpus=1)
set(TEST_CMD_mpi    --p=2)
set(TEST_CMD_mpigpu --p=2 --gpus=1)

set(MPI_CMD_shm )
set(MPI_CMD_shmgpu )
set(MPI_CMD_mpi mpirun -np 4)
set(MPI_CMD_mpigpu mpirun -np 4)

set( TEST_CATEGORIES shm )
if (CHAMELEON_USE_CUDA AND HAVE_CUDA)
   set( TEST_CATEGORIES ${TEST_CATEGORIES} shmgpu )
endif()

# if (CHAMELEON_USE_MPI AND HAVE_MPI)
#    set( TEST_CATEGORIES ${TEST_CATEGORIES} mpi )
#    if (CHAMELEON_USE_CUDA AND HAVE_CUDA)
#       set( TEST_CATEGORIES ${TEST_CATEGORIES} mpigpu )
#    endif()
# endif()

set(TESTLIST 
    gels
    gemm
    getrf_incpiv
    getrf_nopiv
    geqrf
    posv
    potrf
    potri
    )

set(CHAMELEON_PRECISIONS_ZC "c;z")
set(TESTLIST_ZC
    sytrf
    )

foreach(cat ${TEST_CATEGORIES})
    foreach(prec ${RP_CHAMELEON_PRECISIONS})
        foreach(test ${TESTLIST})
            add_test(time_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --check --nowarmup)
        endforeach()
    endforeach()
    foreach(prec ${CHAMELEON_PRECISIONS_ZC})
        foreach(test ${TESTLIST_ZC})
            add_test(time_${cat}_${prec}${test} ${MPI_CMD_${cat}} ./time_${prec}${test}_tile ${TEST_CMD_${cat}} --check --nowarmup)
        endforeach()
    endforeach()
endforeach()

