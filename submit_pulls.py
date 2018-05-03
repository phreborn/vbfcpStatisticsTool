from __future__ import print_function

import argparse
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


__author__ = "Stefan Gadatsch"
__credits__ = ["Stefan Gadatsch"]
__version__ = "0.0.1"
__maintainer__ = "Stefan Gadatsch"
__email__ = "stefan.gadatsch@cern.ch"


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("workspace", type=str, help="Path to workspace to run on.")
    p.add_argument("--folder", type=str, default="test", help="Identifier for the workspace.")
    p.add_argument("--queue", type=str, default="generic", help="Queue to submit to.")
    p.add_argument("--workspaceName", type=str, default="combined", help="Name of the workspace.")
    p.add_argument("--ModelConfigName", type=str, default="ModelConfig", help="Name of the ModelConfig.")
    p.add_argument("--dataName", type=str, default="combData", help="Name of the dataset.")
    p.add_argument("--poi", type=str, default="mu", help="Name of the POI.")
    p.add_argument("--snapshotName", type=str, default="ucmles", help="Name of the snapshot from which all fits start.")
    p.add_argument("--loglevel", type=str, default="INFO", help="Control the printout.")
    p.add_argument("--profile", type=str, default="", help="Parameters that should be profiled.")
    p.add_argument("--fix", type=str, default="", help="Parameters that should be fixed.")
    p.add_argument("--strategy", type=str, default="0", help="Defaul strategy.")
    p.add_argument("--eps", type=str, default="1.0", help="Eps.")
    p.add_argument("--precision", type=str, default="0.01", help="Precision of uncertainty evaluation.")

    args = p.parse_args()

    config = {}

    config["workspace"] = args.workspace
    config["folder"] = args.folder
    config["queue"] = args.queue
    config["workspaceName"] = args.workspaceName
    config["ModelConfigName"] = args.ModelConfigName
    config["dataName"] = args.dataName
    config["poi"] = args.poi
    config["snapshotName"] = args.snapshotName
    config["loglevel"] = args.loglevel
    config["profile"] = args.profile
    config["fix"] = args.fix
    config["strategy"] = args.strategy
    config["eps"] = args.eps
    config["precision"] = args.precision

    if ("<" in config["profile"] or ">" in config["profile"]):
        config["profile"] = "'" + config["profile"] + "'"

    return config


def main(argv):
    config = parse_args(argv)

    workspace = config["workspace"]
    folder = config["folder"]
    queue = config["queue"]
    workspaceName = config["workspaceName"]
    ModelConfigName = config["ModelConfigName"]
    dataName = config["dataName"]
    poi = config["poi"]
    snapshotName = config["snapshotName"]
    loglevel = config["loglevel"]
    profile = config["profile"]
    fix = config["fix"]
    strategy = config["strategy"]
    eps = config["eps"]
    precision = config["precision"]

    os.system("mkdir -vp bsub/%s" % folder)
    os.system("mkdir -vp root-files/%s/pulls" % folder)

    home_folder = os.getcwd()

    f = ROOT.TFile.Open(workspace)
    w = f.Get(workspaceName)
    mc = w.obj(ModelConfigName)
    nuis = mc.GetNuisanceParameters()
    niter = nuis.createIterator()

    var = niter.Next()
    while var:
        varName = str(var.GetName())

        print(varName)

        if "gamma_stat" in varName:
            var = niter.Next()
            continue

        if varName.startswith('scale_'):
            var = niter.Next()
            continue

        if varName.startswith('unconst_'):
            var = niter.Next()
            continue

        if varName.startswith('u_'):
            var = niter.Next()
            continue

        if varName.endswith('_COMB'):
            var = niter.Next()
            continue

        if varName in fix:
            var = niter.Next()
            continue

        submitJob(parameter=varName, folder=folder, workspace=workspace, poi=poi, workspaceName=workspaceName, ModelConfigName=ModelConfigName, dataName=dataName, snapshotName=snapshotName, loglevel=loglevel, profile=profile, strategy=strategy, fix=fix, eps=eps, precision=precision, home_folder=home_folder, queue=queue)

        var = niter.Next()


def submitJob(parameter, folder, workspace, poi, workspaceName, ModelConfigName, dataName, snapshotName, loglevel, profile, strategy, fix, eps, precision, home_folder, queue):
    bsubFileName = "bsub/" + folder + "/" + parameter + ".sh"
    bsubFile = open(bsubFileName, "w")
    text = getJobDef(parameter=parameter, folder=folder, workspace=workspace, poi=poi, workspaceName=workspaceName, ModelConfigName=ModelConfigName, dataName=dataName, snapshotName=snapshotName, loglevel=loglevel, profile=profile, strategy=strategy, fix=fix, eps=eps, precision=precision, home_folder=home_folder, queue=queue)
    bsubFile.write(text)
    bsubFile.close()
    os.system("chmod -R 775 bsub/" + folder)
    command = "bsub < " + bsubFileName
    print(command)
    os.system(command)


def getJobDef(parameter, folder, workspace, poi, workspaceName, ModelConfigName, dataName, snapshotName, loglevel, profile, strategy, fix, eps, precision, home_folder, queue):
    command = ""
    command += "./bin/pulls.exe --input %s" % (workspace)
    if (poi != ""):
        command += " --poi %s" % (poi)
    if (parameter != ""):
        command += " --parameter %s" % (parameter)
    if (workspaceName != ""):
        command += " --workspace %s" % (workspaceName)
    if (ModelConfigName != ""):
        command += " --modelconfig %s" % (ModelConfigName)
    if (dataName != ""):
        command += " --data %s" % (dataName)
    if (snapshotName != ""):
        command += " --snapshot %s" % (snapshotName)
    if (folder != ""):
        command += " --folder %s" % (folder)
    if (loglevel != ""):
        command += " --loglevel %s" % (loglevel)
    if (profile != ""):
        command += " --profile %s" % (profile)
    if (strategy != ""):
        command += " --strategy %s" % (strategy)
    if (fix != ""):
        command += " --fix \"%s\"" % (fix)
    if (eps != ""):
        command += " --eps %s" % (eps)
    if (precision != ""):
        command += " --precision %s" % (precision)
    command += ";\n"
    command = command.replace("--poi ',", "--poi '")

    print(command)

    text = """
#!/bin/bash

#BSUB -J pulls_%s_%s
#BSUB -o bsub/%s/stdout_%s.out
#BSUB -q %s
#BSUB -u $USER@cern.ch
#BSUB -R "select[pool>1000 && mem>2000 && swap>2000]"
#BSUB -R "rusage[mem=2000]"

WORKDIR=$TMPDIR/LSF_$LSB_JOBID
HOMEDIR=%s
OUTDIR=$HOMEDIR
FOLDER=%s

stagein()
{
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
    mkdir -vp ${WORKDIR}
    cd ${HOMEDIR} 2> /dev/null || { echo "The directory does not exist."; exit -1; }

    echo Current folder is
    pwd
    ls -l

    # ATLAS environment
    export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase;
    echo ${ATLAS_LOCAL_ROOT_BASE}
    alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh';
    setupATLAS

    # LCG environment
    lsetup "lcgenv -p LCG_88 x86_64-slc6-gcc62-opt ROOT";

    cd %s;
}

runcode()
{
""" % (folder, parameter, folder, parameter, queue, home_folder, folder, home_folder)

    text += command
    text += """
}

stageout()
{
    cd ${OUTDIR}; ls -l
}

stagein
runcode
stageout

exit
"""
    return text


if __name__ == '__main__':
    exit(main(sys.argv[1:]))
