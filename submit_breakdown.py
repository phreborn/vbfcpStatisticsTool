
from __future__ import print_function

import argparse
import yaml
import ctypes
import os
import subprocess
import sys

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
    if os.path.isfile("lib/libStatisticsTools.so"):
        path = "lib/libStatisticsTools.so"
    elif os.path.isfile("lib/libStatisticsTools.dylib"):
        path = "lib/libStatisticsTools.dylib"

    ROOT.gSystem.Load(path)
    tools = ctypes.cdll.LoadLibrary(path)
    prototype = ctypes.CFUNCTYPE(ctypes.c_void_p)
    loadCustom = prototype(('loadCustom', tools))
    loadCustom()
except Exception:
    print("Could not load shared library. Make sure that it was compiled.")


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=str, required=True, help="File to run over.")
    p.add_argument("--poi", type=str, required=True, help="POIs to measure.")
    p.add_argument("--snapshot", type=str, default="nominalNuis", help="Initial snapshot.")
    p.add_argument("--folder", type=str, required=True, help="Output folder.")
    p.add_argument("--profile", type=str, help="Parameters to profile.")
    p.add_argument("--fix", type=str, help="Parameters to fix.")
    p.add_argument("--workspace", type=str, default="combined", help="WS to grab.")
    p.add_argument("--modelconfig", type=str, default="ModelConfig", help="MC to load.")
    p.add_argument("--data", type=str, help="Data to use.")
    p.add_argument("--minimizerType", type=str, help="Minimizer type.")
    p.add_argument("--minimizerAlgo", type=str, help="Minimizer algorithm.")
    p.add_argument("--strategy", type=int, help="Default strategy.")
    p.add_argument("--numCPU", type=int, help="Number of CPUs.")
    p.add_argument("--binned", type=int, help="Binned likelihood.")
    p.add_argument("--starfix", type=int, help="Fix StarMomentMorph cache.")
    p.add_argument("--multifix", type=int, help="Fix MultiPdf level 2.")
    p.add_argument("--precision", type=float, help="Precision for scan.")
    p.add_argument("--eps", type=float, help="Convergence criterium.")
    p.add_argument("--offset", type=int, help="Offset likelihood.")
    p.add_argument("--optimize", type=int, help="Optimize constant terms.")
    p.add_argument("--loglevel", type=str, help="Control verbosity.")
    p.add_argument("--classification", type=str, help="Definition of uncertainty categories.")
    p.add_argument("--subtractFromTotal", action='store_true', help="Subtract uncertainties from total.")
    p.add_argument("--doIndividual", action='store_true', help="Compute uncertainty for individual sources.")

    # submission arguments
    p.add_argument("--queue", type=str, default='short', help="queue to submit to.")
    p.add_argument("--submit", action='store_true', help="queue to submit to.")

    args = p.parse_args(argv)

    with open(args.classification, 'r') as inf:
        yaml_config = yaml.load(inf)
    submission_args = dict()
    submission_args["queue"] = args.queue
    submission_args["submit"] = args.submit
    del args.queue
    del args.submit

    return args, submission_args, yaml_config


def getJobDef(replacedir):
    text = """
#!/bin/bash

#PBS -j oe
#PBS -N pulls_{folder}_{category}
#PBS -o qsub/{folder}/stdout_{category}.out
#PBS -q {queue}
#PBS -u $USER
#PBS -l pvmem=16000MB

WORKDIR=$TMPDIR/PBS_$PBS_JOBID
HOMEDIR={home_folder}
OUTDIR=$HOMEDIR
FOLDER={folder}

stagein()
{{
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
    mkdir -vp ${{WORKDIR}}
    cd ${{HOMEDIR}} 2> /dev/null || {{ echo "The directory does not exist."; exit -1; }}

    echo Current folder is
    pwd
    ls -l

    # ATLAS environment
    source /project/atlas/nikhef/cvmfs/setup.sh
    source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase;
    echo ${{ATLAS_LOCAL_ROOT_BASE}}

    # LCG environment
    lsetup "lcgenv -p LCG_93 x86_64-slc6-gcc62-opt ROOT";

    cd {home_folder};
}}

runcode()
{{
    {command}
}}
stageout()
{{
    cd ${{OUTDIR}}; ls -l
}}

stagein
runcode
stageout

exit
"""
    return text.format(**replacedir)


def buildCommand(base, these_args):
    base += " "
    for key, value in vars(these_args).iteritems():
        if value is not None:
            if isinstance(value, bool):
                if value:
                    base += "--" + key + " "
            elif isinstance(value, str) or isinstance(value, int) or isinstance(value, float):
                base += "--" + key + " " + str(value) + " "
    return base


def main(argv):
    args, submission_args, yaml_config = parse_args(argv)

    submission_path = os.path.join("qsub", args.folder)
    if not os.path.exists(submission_path):
        os.makedirs(submission_path)

    cmd = buildCommand("./bin/run_breakdown.exe", args)

    replacedir = dict()
    replacedir['folder'] = args.folder
    replacedir['home_folder'] = os.getcwd()
    replacedir['queue'] = submission_args["queue"]

    for category in yaml_config.keys():
        this_command = cmd + " --category " + category
        print(this_command)
        replacedir['command'] = this_command
        replacedir['category'] = category
        job_content = getJobDef(replacedir)
        path_to_job = os.path.join(submission_path, category + ".sh")
        with open(path_to_job, "w") as f:
            f.write(job_content)
        if submission_args["submit"]:
            submission_cmd = "qsub " + path_to_job
            os.system(submission_cmd)


if __name__ == "__main__":
    exit(main(sys.argv[1:]))

    cmd = "./bin/run_breakdown.exe --input ../Workspaces_couplings/v4/combination/WS-Comb-mu.root --workspace combWS --data combData --poi 'mu' --eps 0.1 --loglevel DEBUG --strategy 0 --classification ../config/classification_couplings.yaml --subtractFromTotal > breakdown_couplings_fine_180609.txt"
