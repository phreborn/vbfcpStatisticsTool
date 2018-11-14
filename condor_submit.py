#!/usr/bin/env python
from __future__ import print_function

import sys, os, glob
import fileinput
import argparse
import ctypes
import itertools
from decimal import Decimal

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

__author__ = "Stefan Gadatsch"
__credits__ = ["Stefan Gadatsch"]
__version__ = "0.0.1"
__maintainer__ = "Stefan Gadatsch"
__email__ = "stefan.gadatsch@cern.ch"
# migration to Condor: Tim Michael Heinz Wolf

# set general parameters
memory = 8000
MaxRuntime = 36*60*60
mail = "tim.michael.heinz.wolf@cern.ch"

# template scripts
sub_temp = "template_submit.sub"
template_pulls = "template_pulls.sh"
template_scan = "template_scan.sh"

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
    p.add_argument("--scanRange", type=str, default="0.0:2.0", help="Range which should be scanned.")
    p.add_argument("--bins", type=str, default="10", help="Number of bins for the scanned range.")
    p.add_argument("--pointsPerJob", type=str, default="1", help="Points to scan per job.")

    p.add_argument("--template", type=str, required=True)
    p.add_argument("--debug", dest='debug', action='store_true')
    p.add_argument("--missing", dest='missing', action='store_true')

    args = p.parse_args()

    config = vars(p.parse_args())

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
    scanRange = config["scanRange"]
    bins = config["bins"]
    pointsPerJob = config["pointsPerJob"]
    template = config["template"]
    debug = config["debug"]
    missing = config["missing"]

    submission_dir = "bsub/{0}".format(folder)
    resultdir = "root-files/{0}/{1}".format(folder, template)
    os.system("mkdir -vp " + submission_dir)
    os.system("mkdir -vp " + resultdir)

    home_folder = os.getcwd()

    f = ROOT.TFile.Open(workspace)
    w = f.Get(workspaceName)
    mc = w.obj(ModelConfigName)
    nuis = mc.GetNuisanceParameters()
    niter = nuis.createIterator()

    if template == "pulls":
      var = niter.Next()
      jobNames = str()
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

          if missing and os.path.exists("root-files/" + folder + "/pulls/" + varName + ".root"):
              var = niter.Next()
              continue

          parameters_to_pass = ["workspace", "workspaceName", "poi", "ModelConfigName", "dataName", "snapshotName", "folder", "loglevel", "profile", "strategy", "fix", "eps", "precision"]
          parameter_dict = dict()
          for parameter in parameters_to_pass:
              parameter_dict[parameter] = config[parameter]

          jobName = getJobsPulls(debug, varName, folder, home_folder, submission_dir, parameter_dict)
          jobNames += jobName + "\n"

          var = niter.Next()
      submitScript(template, sub_temp, submission_dir, jobNames, memory, mail, MaxRuntime)

    elif template == "scan":
        ranges = scanRange.split(',')
        allbins = bins.split(',')

        allpoiVals = []
        for thisRange, thisbins in itertools.izip(ranges, allbins):
            doLog = False
            if "L" in thisbins:
                doLog = True
                thisbins = thisbins.replace("L", "")
            poiVals = set()
            low, high = thisRange.split(':')
            counter = Decimal(low)
            step = (Decimal(high) - Decimal(low)) / Decimal(thisbins)
            while counter <= Decimal(high):
                thisVal = counter
                if doLog:
                    thisVal = pow(10, counter)
                thiscounter = "%.13f" % thisVal
                poiVals.add(thiscounter)
                counter += step
            if doLog:
                poiVals.add("%.13f" % pow(10, Decimal(high)))
            else:
                poiVals.add(high)
            allpoiVals.append(poiVals)

        a = list(itertools.product(*allpoiVals))

        totaljobs = 0
        for thisA in a:
            totaljobs += 1

        submitArray = []
        pointsInArray = 0
        submittedalready = 0
        jobNames = str()
        for thisA in a:
            myString = ",".join(thisA)
            submitArray.append(myString)
            pointsInArray += 1

            if pointsInArray == int(pointsPerJob):
                jobName = getJobsScan(debug, submitArray, folder, workspace, poi, workspaceName, ModelConfigName, dataName, snapshotName, loglevel, profile, strategy, fix, eps, precision, home_folder, submission_dir)
                jobNames += jobName + "\n"

                submittedalready += pointsInArray
                pointsInArray = 0
                submitArray = []

        if pointsInArray == totaljobs-submittedalready and totaljobs-submittedalready > 0:
            jobName = getJobsScan(debug, submitArray, folder, workspace, poi, workspaceName, ModelConfigName, dataName, snapshotName, loglevel, profile, strategy, fix, eps, precision, home_folder, submission_dir)
            jobNames += jobName + "\n"

            submittedalready += pointsInArray
            pointsInArray = 0
            submitArray = []
        submitScript(template, sub_temp, submission_dir, jobNames, memory, mail, MaxRuntime)


def getJobsPulls(debug, parameter, folder, home_folder, submitdir, parameter_dict):
    jobName = "pulls_{0}".format(parameter)
    parameter_dict["parameter"] = parameter
    command = "./bin/pulls.exe --input {workspace} --poi {poi} --workspace {workspaceName} --modelconfig {ModelConfigName} --data {dataName} --snapshot {snapshotName} --folder {folder} --loglevel {loglevel} --strategy {strategy} --eps {eps} --precision {precision}"
    command = command.format(**parameter_dict)
    if profile != "":
        command += " --profile {}".format(parameter_dict["profile"])
    if fix != "":
        command += " --fix {}".format(parameter_dict["fix"])
    command += " 2>&1;\n"

    with open(template_pulls) as f:
        text = f.read()
    replacedict = dict()
    replacedict["FOLDER"] = folder
    replacedict["HOME_FOLDER"] = home_folder
    replacedict["COMMAND"] = command
    text = text.format(**replacedict)

    scriptName = jobName + ".sh"
    scriptToCreate = os.path.join(submitdir, scriptName)
    with open(scriptToCreate, "w") as f:
        f.write(text)
    os.system("chmod +x " + scriptToCreate)
    print(command)
    if debug:
        exit(0)
    return jobName


def getJobsScan(debug, pointsInArray, folder, workspace, poi, workspaceName, ModelConfigName, dataName, snapshotName, loglevel, profile, strategy, fix, eps, precision, home_folder, submission_dir):
    print(pointsInArray)
    pois = poi.split(',')
    poiVal = pointsInArray[0]
    poiVals = poiVal.split(',')
    jobName = "scan_{0}_{1}".format(poi, poiVal)

    command = ""
    for thisPoiVals in pointsInArray:
        thisPoiVal = thisPoiVals.split(',')
        command += "./bin/scan.exe --input %s" % (workspace)

        if (poi != ""):
            command += " --poi '"
            for poi_tmp, val_tmp in itertools.izip(pois, thisPoiVal):
                command += ",%s[%s]" % (poi_tmp, val_tmp)
            command += "'"
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

    with open(template_scan, "r") as f:
        text = f.read()
    replacedict = dict()
    replacedict["FOLDER"] = folder
    replacedict["HOME_FOLDER"] = home_folder
    replacedict["COMMAND"] = command
    text = text.format(**replacedict)

    jobName = str()
    for thisPoi, thisPoiVal in itertools.izip(pois, poiVals):
        jobName += thisPoi + thisPoiVal
    subFileName = os.path.join(submission_dir, jobName + ".sh")
    with open(subFileName, "w") as f:
        f.write(text)
    if debug:
        exit(0)

    return jobName


def submitScript(template, sub_temp, submission_dir, jobNames, memory, mail, MaxRuntime):
    script_sub = os.path.join(submission_dir, "run" + template + ".sub")
    with open(sub_temp) as f:
        text = f.read()
    replacedict = dict()
    replacedict["RUN"] = os.path.join(submission_dir, "$(job)")
    replacedict["MEM"] = str(memory)
    replacedict["MAXRUNTIME"] = str(MaxRuntime)
    replacedict["MAIL"] = mail
    replacedict["JOBS"] = jobNames[:-1]
    text = text.format(**replacedict)

    with open(script_sub, "w") as f:
        f.write(text)
    CondorSubmitCommand = "condor_submit " + script_sub
    print(CondorSubmitCommand)
    os.system(CondorSubmitCommand)


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
