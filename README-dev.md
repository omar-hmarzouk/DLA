README for developers
====================

This page is dedicated to rules and conventions that Chameleon's
developers must follow and that should be read by contributors.

Gitlab flow
---------------------

TODO - Write here how to use git and gitlab to push modifications into
Chameleon.

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

TODO - Write here how to generate the documentation.

Packaging
---------------------

TODO - Write here the steps to create a new release:
 - what needs to be checked (ctest, docummentation, version, ...)
 - how to generate the package .tar.gz
 - what and where to publish externally

Naming and writting conventions
---------------------

TODO - Write here the rules and conventions used in the source code.