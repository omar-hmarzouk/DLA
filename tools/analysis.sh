#!/bin/bash

# Performs an analysis of Chameleon source code:
# - we consider to be in Chameleon's source code root
# - we consider having the coverage file chameleon_coverage.xml in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in build
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of chameleon and make cppcheck analysis too long
./tools/find_sources.sh

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UCHAMELEON_USE_OPENCL -UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"
# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingIncludeSystem ${UNDEFINITIONS} --file-list=./filelist.txt 2> chameleon_cppcheck.xml
# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > chameleon_rats.xml

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
sonar.cxx.compiler.reportPath=chameleon_starpu.log, chameleon_starpu_simgrid.log, chameleon_quark.log
sonar.cxx.coverage.reportPath=chameleon_coverage.xml
sonar.cxx.cppcheck.reportPath=chameleon_cppcheck.xml
sonar.cxx.rats.reportPath=chameleon_rats.xml
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X -Dsonar.login=077ef775a8b4fe9ba722497ef0511ca6dcfb3fcd > sonar.log
