from __future__ import print_function

import argparse
import ctypes
import os
import subprocess
import sys
from decimal import Decimal
import itertools

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
__credits__ = ["Stefan Gadatsch", "Tim Wolf"]
__version__ = "0.0.1"
__maintainer__ = "Stefan Gadatsch"
__email__ = "stefan.gadatsch@cern.ch"


def parse_args(argv):
    p = argparse.ArgumentParser()
    p.add_argument("workspace", type=str, help="Path to workspace to run on.")
    p.add_argument("--poi", type=str, default="mu", help="Name of the POI.")
    p.add_argument("--snapshotName", type=str, default="ucmles", help="Name of the snapshot from which all fits start.")
    p.add_argument("--folder", type=str, default="test", help="Identifier for the workspace.")
    p.add_argument("--profile", type=str, default="", help="Parameters that should be profiled.")
    p.add_argument("--fix", type=str, default="", help="Parameters that should be fixed.")
    p.add_argument("--workspaceName", type=str, default="combined", help="Name of the workspace.")
    p.add_argument("--ModelConfigName", type=str, default="ModelConfig", help="Name of the ModelConfig.")
    p.add_argument("--dataName", type=str, default="combData", help="Name of the dataset.")
    p.add_argument("--strategy", type=str, default="0", help="Defaul strategy.")
    p.add_argument("--precision", type=str, default="0.01", help="Precision of uncertainty evaluation.")
    p.add_argument("--eps", type=str, default="1.0", help="Eps.")
    p.add_argument("--loglevel", type=str, default="INFO", help="Control the printout.")
    p.add_argument("--queue", type=str, default="espresso", help="Queue to submit to.")
    p.add_argument("--scanRange", type=str, default="0.0:2.0", help="Range which should be scanned.")
    p.add_argument("--bins", type=str, default="10", help="Number of bins for the scanned range.")
    p.add_argument("--pointsPerJob", type=str, default="1", help="Points to scan per job.")

    p.add_argument("--mail", type=str, default="$USER@cern.ch", help="Notifications.")
    p.add_argument("--mem", type=str, default="8000", help="Memory.")
    p.add_argument("--pool", type=str, default="1000", help="Pool.")
    p.add_argument("--rmem", type=str, default="2000", help="RMemory.")
    p.add_argument("--swap", type=str, default="8000", help="Swap.")
    p.add_argument("--time", type=str, default="129600", help="Run time.")

    batchsystem = p.add_mutually_exclusive_group()
    batchsystem.add_argument('--condor', help='CERN condor', action='store_true')
    batchsystem.add_argument('--lsf', help='CERN LSF', action='store_true')

    args = p.parse_args()

    if ("<" in args.profile or ">" in args.profile):
        args.profile = "'" + args.profile + "'"

    return args


def main(argv):
    # Parse commandline arguments
    args = parse_args(argv)

    if args.condor:
        args.cluster = "condor"
    elif args.lsf:
        args.cluster = "lsf"
        print("CERN lsf is not active anymore, consider using condor")
        sys.exit(-1)
    else:
        print("Must specify cluster type: condor or lsf")
        sys.exit(-1)

    args.home_folder = os.getcwd()
    args.submission_folder = "{0}/{1}".format(args.cluster, args.folder)
    args.result_folder = "root-files/{0}/scan".format(args.folder)

    # Create directories
    os.system("mkdir -vp " + args.submission_folder)
    os.system("mkdir -vp " + args.result_folder)

    # Get POI values
    poi_values = get_poi_values(args)

    # Split in requested number of jobs and write corresponding job scripts
    jobs = collect_jobs(args, poi_values)

    # Submit jobs
    if args.cluster == "lsf":
        submit_lsf(args, jobs)
    elif args.cluster == "condor":
        submit_condor(args, jobs)


def get_poi_values(args):
    ranges = args.scanRange.split(',')
    allbins = args.bins.split(',')

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

    poi_values = list(itertools.product(*allpoiVals))

    return poi_values


def collect_jobs(args, poi_values):
    submit_array = []
    points_in_array = 0
    submitted_already = 0
    total_jobs = 0
    jobs = str()

    for this_val in poi_values:
        total_jobs += 1

    for this_val in poi_values:
        myString = ",".join(this_val)
        submit_array.append(myString)
        points_in_array += 1

        if points_in_array == int(args.pointsPerJob):
            job = write_job(submit_array, args)
            jobs += job + "\n"
            submitted_already += points_in_array
            points_in_array = 0
            submit_array = []

    if points_in_array == total_jobs - submitted_already and total_jobs - submitted_already > 0:
        job = write_job(submit_array, args)
        jobs += job + "\n"
        submitted_already += points_in_array

    return jobs


def write_job(pointsInArray, args):
    pois = args.poi.split(',')
    poiVal = pointsInArray[0]
    poiVals = poiVal.split(',')

    command = get_command(pointsInArray, args)

    if args.cluster == "lsf":
        with open("share/lsf_submit.sub") as f:
            directives = f.read()

        replacedict = dict()
        replacedict["POI"] = args.poi
        replacedict["POIVAL"] = poiVal
        replacedict["FOLDER"] = args.folder
        replacedict["QUEUE"] = args.queue
        replacedict["MAIL"] = args.mail
        replacedict["POOL"] = args.pool
        replacedict["MEM"] = args.mem
        replacedict["SWAP"] = args.swap
        replacedict["RMEM"] = args.rmem

        directives = directives.format(**replacedict)
    elif args.cluster == "condor":
        directives = ""

    with open("share/batchjob.sh") as f:
        batchjob = f.read()

    replacedict = dict()
    replacedict["DIRECTIVES"] = directives
    replacedict["FOLDER"] = args.folder
    replacedict["HOME_FOLDER"] = args.home_folder
    replacedict["COMMAND"] = command

    batchjob = batchjob.format(**replacedict)

    job_name = "scan_{0}_{1}".format(args.poi, poiVal)
    batchjob_file = os.path.join(args.submission_folder, job_name + ".sh")

    with open(batchjob_file, "w") as f:
        f.write(batchjob)

    if args.cluster == "lsf":
        return batchjob_file
    elif args.cluster == "condor":
        return job_name


def get_command(vals, args):
    command = ""

    for thisPoiVals in vals:
        thisPoiVal = thisPoiVals.split(',')

        command += "./bin/scan.exe --input %s" % (args.workspace)

        if (args.poi != ""):
            command += " --poi '"
            for poi_tmp, val_tmp in itertools.izip(args.poi.split(','), thisPoiVal):
                command += ",%s[%s]" % (poi_tmp, val_tmp)
            command += "'"
        if (args.workspaceName != ""):
            command += " --workspace %s" % (args.workspaceName)
        if (args.ModelConfigName != ""):
            command += " --modelconfig %s" % (args.ModelConfigName)
        if (args.dataName != ""):
            command += " --data %s" % (args.dataName)
        if (args.snapshotName != ""):
            command += " --snapshot %s" % (args.snapshotName)
        if (args.folder != ""):
            command += " --folder %s" % (args.folder)
        if (args.loglevel != ""):
            command += " --loglevel %s" % (args.loglevel)
        if (args.profile != ""):
            command += " --profile %s" % (args.profile)
        if (args.strategy != ""):
            command += " --strategy %s" % (args.strategy)
        if (args.fix != ""):
            command += " --fix \"%s\"" % (args.fix)
        if (args.eps != ""):
            command += " --eps %s" % (args.eps)
        if (args.precision != ""):
            command += " --precision %s" % (args.precision)
        command += ";\n    "
        command = command.replace("--poi ',", "--poi '")

    return command[:-5]


def submit_lsf(args, jobs):
    for submission_script in jobs.split("\n"):
        command = "bsub < " + submission_script

        print(command)

        os.system(command)


def submit_condor(args, jobs):
    submission_script = os.path.join(args.submission_folder, "scan.sub")

    with open("share/condor_submit.sub") as f:
        text = f.read()

    replacedict = dict()
    replacedict["RUN"] = os.path.join(args.submission_folder, "$(job)")
    replacedict["MEM"] = args.mem
    replacedict["MAXRUNTIME"] = args.time
    replacedict["MAIL"] = args.mail
    replacedict["JOBS"] = jobs[:-1]

    text = text.format(**replacedict)

    with open(submission_script, "w") as f:
        f.write(text)

    command = "condor_submit " + submission_script

    print(command)

    os.system(command)


if __name__ == '__main__':
    exit(main(sys.argv[1:]))
