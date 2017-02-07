README for developers
====================

This page is dedicated to rules and conventions that Chameleon's
developers must follow and that should be read by contributors.

Gitlab flow: how to contribute to Chameleon
---------------------

Please read and follow guidelines given here:
https://gitlab.inria.fr/solverstack/chameleon/blob/master/CONTRIBUTING.md

### Update submodules

Chameleon git project depends on a "submodule" git, located in
cmake_modules/morse_cmake and hosted here
https://gitlab.inria.fr/solverstack/morse_cmake.

To update this submodule to the last development state, follow these
steps:

```
#!shell
git submodule update --remote cmake_modules/morse_cmake
git commit cmake_modules/morse_cmake -m "update morse_cmake submodule"
git push --recurse-submodules=check
```

Documentation
---------------------

### Generate the documentation

#### Prerequisites

To generate the documentation you need to have
[Doxygen](http://www.stack.nl/~dimitri/doxygen/) and
[Texinfo](https://www.gnu.org/software/texinfo/) installed on your
system.

For example, on Debian systems:
```
#!shell
sudo apt install doxygen texinfo
```

#### configure + make documentation

Enter into the Chameleon's source root directory and configure with
*CHAMELEON_ENABLE_DOCS=ON*, then generate the documentation with `make
docs`

```
#!shell
cd chameleon
mkdir build && cd build
cmake .. -DCHAMELEON_ENABLE_DOCS=ON
make docs
```

### Rules about source code documentation

TODO - Write here the rules to document source code inside chameleon.

Packaging
---------------------

TODO - Write here the steps to create a new release:
 - what needs to be checked (ctest, docummentation, version, ...)
 - how to generate the package .tar.gz
 - what and where to publish externally

Naming and writting conventions
---------------------

TODO - Write here the rules and conventions used in the source code.