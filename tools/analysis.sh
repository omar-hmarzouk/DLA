#!/bin/bash

# Performs an analysis of Chameleon source code
# We consider to be in Chameleon's source code root

# build with proper options
mkdir -p build
cd build
rm * -rf
export CFLAGS="-O0 -g -fPIC -fdiagnostics-show-option --coverage -fno-inline -Wall"
export LDFLAGS="--coverage"
cmake .. -DCHAMELEON_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=$PWD/install -DCMAKE_C_FLAGS="$CFLAGS" -DCMAKE_EXE_LINKER_FLAGS="$LDFLAGS" -DCMAKE_VERBOSE_MAKEFILE=ON
make -j5 | tee ../chameleon-build.log

# run tests
STARPU_SILENT=1 ctest --no-compress-output || /usr/bin/true

# capture coverage
lcov --directory . --capture --output-file ../chameleon.lcov
cd ..
genhtml -o coverage chameleon.lcov
lcov_cobertura.py chameleon.lcov --output chameleon-coverage.xml

# filter sources:
# - consider generated files in build
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long
export SOURCES_TO_ANALYZE=`bash -c 'find ./build -type d -name CMakeFiles -prune -o -type d -name Testing -prune -o -type f -regex ".*\.c\|.*\.h" -print && \
                                    find . -path ./build -prune -o -type f -regex "^[^z]*\.c" -print && \
                                    find . -path ./build -prune -o -type f -regex "^[^z]*\.h" ! -name 'lapacke*.h' ! -name 'cblas*.h' -print | xargs'`
# actually we need to remove cblas/lapacke headers (cppcheck analysis too long)
# we will get them back after cppcheck
rm coreblas/include/cblas.h coreblas/include/lapacke.h coreblas/include/lapacke_config.h coreblas/include/lapacke_mangling.h

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_CUBLAS_V2 -UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"
# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem ${UNDEFINITIONS} ${SOURCES_TO_ANALYZE} 2> chameleon-cppcheck.xml
# run rats analysis
rats -w 3 --xml ${SOURCES_TO_ANALYZE} > chameleon-rats.xml

# get back cblas/lapacke headers after static analysis
git checkout coreblas/include/cblas.h coreblas/include/lapacke.h coreblas/include/lapacke_config.h coreblas/include/lapacke_mangling.h

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://hpclib-sed.bordeaux.inria.fr/sonarqube-dev
sonar.links.homepage=https://gitlab.inria.fr/solverstack/chameleon
sonar.links.ci=https://gitlab.inria.fr/solverstack/chameleon/pipelines
sonar.links.scm=https://gitlab.inria.fr/solverstack/chameleon/
sonar.projectKey=chameleon
sonar.projectName=Chameleon
sonar.projectDescription=Dense linear algebra subroutines for heterogeneous and distributed architectures
sonar.projectVersion=master
sonar.language=c++
sonar.sources=`bash -c 'echo $SOURCES_TO_ANALYZE | sed -e "s/ /, /g"'`
sonar.sourceEncoding=UTF-8
sonar.cxx.compiler.charset=UTF-8
sonar.cxx.compiler.regex=^(.*):([0-9]+):[0-9]+: warning: (.*)\[(.*)\]$
sonar.cxx.compiler.reportPath=chameleon-build.log
sonar.cxx.coverage.reportPath=chameleon-coverage.xml
sonar.cxx.cppcheck.reportPath=chameleon-cppcheck.xml
sonar.cxx.rats.reportPath=chameleon-rats.xml
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X -Dsonar.login=077ef775a8b4fe9ba722497ef0511ca6dcfb3fcd > sonar.log
