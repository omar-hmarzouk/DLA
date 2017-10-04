#!/bin/bash

# Performs an analysis of Chameleon source code
# We consider to be in Chameleon's source code root

# build with proper options
#mkdir -p build
#cd build
#rm * -rf
#cmake .. -DCHAMELEON_USE_MPI=ON -DCMAKE_INSTALL_PREFIX=$PWD/install -DCMAKE_VERBOSE_MAKEFILE=ON -DMORSE_ENABLE_WARNING=ON -DMORSE_ENABLE_COVERAGE=ON
#make -j5 | tee ../chameleon-build.log

# run tests
#STARPU_SILENT=1 ctest --no-compress-output || /usr/bin/true

# capture coverage
#lcov --directory . --capture --output-file ../chameleon.lcov
lcov --directory build --capture --output-file chameleon.lcov
#cd ..
#genhtml -o coverage chameleon.lcov
lcov_cobertura.py chameleon.lcov --output chameleon-coverage.xml

# filter sources:
# - consider generated files in build
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long
./tools/find_sources.sh

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"
# run cppcheck analysis
#cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem ${UNDEFINITIONS} --file-list=./filelist.txt 2> chameleon-cppcheck.xml
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem ${UNDEFINITIONS} `cat filelist.txt` 2> chameleon-cppcheck.xml
# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > chameleon-rats.xml

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
sonar.sources=build, compute, control, coreblas, example, include, runtime, testing, timing
sonar.inclusions=`cat filelist.txt | xargs echo | sed 's/ /, /g'`
sonar.sourceEncoding=UTF-8
sonar.cxx.compiler.charset=UTF-8
sonar.cxx.compiler.parser=GCC
sonar.cxx.compiler.regex=^(.*):(\\d+):\\d+: warning: (.*)\\[(.*)\\]$
sonar.cxx.compiler.reportPath=chameleon-build.log
sonar.cxx.coverage.reportPath=chameleon-coverage.xml
sonar.cxx.cppcheck.reportPath=chameleon-cppcheck.xml
sonar.cxx.rats.reportPath=chameleon-rats.xml
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X -Dsonar.login=077ef775a8b4fe9ba722497ef0511ca6dcfb3fcd > sonar.log
