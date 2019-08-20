#!/bin/bash

{DIRECTIVES}

HOMEDIR={HOME_FOLDER}
OUTDIR=$HOMEDIR
FOLDER={FOLDER}

stagein()
{{
    echo "Running in: $PWD  @ $HOSTNAME"
    echo "Time : " $(date -u)

    uname -a
    ulimit -S -s 20000
    ulimit -Sc 0
    ulimit -Hc 0
    ulimit -c 0
    ulimit -d unlimited
    ulimit -f unlimited
    ulimit -l unlimited
    ulimit -n unlimited
    ulimit -s unlimited
    ulimit -t unlimited
    cd ${{HOMEDIR}} 2> /dev/null || {{ echo "The directory does not exist."; exit -1; }}

    echo Current folder is
    pwd
    ls -l

    echo "Setting up environment..."
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
    source ${{ATLAS_LOCAL_ROOT_BASE}}/user/atlasLocalSetup.sh

    setupATLAS
    source setup.sh

    cd $HOMEDIR;
}}

runcode()
{{
    {COMMAND}
}}

stageout()
{{
    cd ${{OUTDIR}}; ls -l
}}

stagein
runcode
stageout

exit
