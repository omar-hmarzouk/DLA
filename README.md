Chameleon: A dense linear algebra software for heterogeneous architectures
====================

Chameleon is a C library providing parallel algorithms to perform
BLAS/LAPACK operations exploiting fully modern architectures.

Chameleon dense linear algebra software relies on sequential
task-based algorithms where sub-tasks of the overall algorithms are
submitted to a Runtime system. Such a system is a layer between the
application and the hardware which handles the scheduling and the
effective execution of tasks on the processing units. A Runtime system
such as [StarPU](http://starpu.gforge.inria.fr/) is able to manage
automatically data transfers between not shared memory area
(CPUs-GPUs, distributed nodes).

This kind of implementation paradigm allows to design high performing
linear algebra algorithms on very different type of architecture:
laptop, many-core nodes, CPUs-GPUs, multiple nodes. For example,
Chameleon is able to perform a Cholesky factorization
(double-precision) at 80 TFlop/s on a dense matrix of order 400 000
(i.e. 4 min). Chameleon is a sub-project of
[MORSE](http://icl.cs.utk.edu/morse/) specifically dedicated to dense
linear algebra.

Get Chameleon
---------------------

To use last development states of Chameleon, please clone the master
branch hosted here, *i.e.*:

    git clone git@gitlab.inria.fr:solverstack/chameleon.git

Last releases of Chameleon are hosted on the
[gforge.inria.fr](https://gforge.inria.fr/frs/?group_id=2884) for now.
Future releases will be available on this gitlab project.

Documentation
---------------------

There is no up-to-date documentation of Chameleon.  We would like to
provide a doxygen documentation hosted on
[gitlab](https://about.gitlab.com/2016/04/07/gitlab-pages-setup/) in
the future.

The documentation of Chameleon's last release is available here:
[chameleon-0.9.1
documentation](http://morse.gforge.inria.fr/chameleon/0.9.1/chameleon_users_guide-0.9.1.html)

Installation
---------------------

### Build and install with CMake

Chameleon can be built using [CMake](https://cmake.org/). This
installation requires to have some library dependencies already
installed on the system.

Please refer to the
[chameleon-0.9.1](http://morse.gforge.inria.fr/chameleon/0.9.1/chameleon_users_guide-0.9.1.html#Installing-CHAMELEON)
to get configuration information.

### Distribution of Chameleon
To get support to install a full distribution (Chameleon +
dependencies) we encourage users to use the morse branch of
**Spack**.

Please read these documentations:

* [Spack Morse](http://morse.gforge.inria.fr/spack/spack.html).
* [Section Chameleon](http://morse.gforge.inria.fr/spack/spack.html#sec-2-1).

Get involved!
---------------------

### Mailing list

TODO

### Contributions

https://gitlab.inria.fr/solverstack/chameleon/blob/master/CONTRIBUTING.md

### Authors

First, since the Chameleon library started as an extension of the PLASMA library 
to support multiple runtime systems, all developpers of the PLASMA library are 
developpers of the Chameleon library.

The following people contributed to the development of Chameleon:
  * Emmanuel Agullo, PI
  * Olivier Aumage
  * Cedric castagnede
  * Terry Cojean
  * Mathieu Faverge, PI
  * Nathalie Furmento
  * Reazul Hoque, PaRSEC support
  * Gregoire Pichon, Two-sided algorithms
  * Florent Pruvost, PI, Cmake and spack 
  * Marc Sergent
  * Guillaume Sylvand
  * Samuel Thibault, StarPU support
  * Omar Zenati

If we forgot your name, please let us know that we can fix that mistake.

### Citing Chameleon

TODO

### Licence

https://gitlab.inria.fr/solverstack/chameleon/blob/master/LICENCE.txt
