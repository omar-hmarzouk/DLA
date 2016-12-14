#!/bin/sh

cd spack
export SPACK_ROOT=$PWD
. $SPACK_ROOT/share/spack/setup-env.sh
cd -
spack install -v openblas
spack install -v starpu
spack load openblas
spack load starpu
