#!/usr/bin/env python3
import sys, os, glob, datetime, json
from pprint import pprint
import numpy as np
import pandas as pd
from shutil import move


def main(argv):
    """
    Automatic processing of COHERENT NaI(Tl) crystal data.
    Uses "runDB.json" as a database.  When a crystal is scanned,
    we input its serial number and run numbers, and run this code
    to process the data and do file handling (moving)
    Usage:
    $ python auto_process.py [options]
        default (no arguments): process all runs on this computer
        -crys [serial no]: process data from a particular crystal
        -sync : upload data to cenpa-rocks and delete the raw ORCA files
    """

    # -- load our run database and make it global --
    global runDB
    with open("runDB.json") as f:
        runDB = json.load(f)

    # -- parse command line args --
    for i, opt in enumerate(argv):

        if opt=="-crys":
            sn = argv[i+1]
            process_crystal(sn)
            exit()

    sns = [id for id in runDB.keys() if "UW" in id]
    print(sns)




def process_crystal(sn):
    print(sn)


def sandbox():
    """ deleteme! """

    with open("runDB.json") as f:
        runDB = json.load(f)

    runs_to_process = runDB["UW3"]["voltage"]

    raw_files = list(set(glob.glob("{}/**/**/Data/*Run*".format(raw_path), recursive=True)))

    for f in raw_files:
        run = int(f.split('Run')[-1])
        if run in runs_to_process:

            mod_path = "/Users/ccenpa/Desktop/coherent/ORCA\ Experiments"

            tmp = f.split("/")
            fname = tmp[-1]
            mon = tmp[-3]
            yr = tmp[-4]

            newfile = "{}/{}/{}/Data/{}".format(mod_path,yr,mon,fname)

            # print(newfile)

            cmd = "./getSpectrum {}".format(newfile)
            sh(cmd)

            # out_file = "NaI_ET_run{}.root".format(run)
            # os.mkdir("{}/UW3".format(built_path))
            # os.mkdir("{}/UW3/voltage/".format(built_path))
            # os.mkdir("{}/UW3/voltage/600V".format(built_path))
            # move(out_file, "{}/UW3/voltage/600V/{}".format(built_path, out_file))

            # nope, do this on rocks
            # sh("gzip {}".format(f))


            exit()

    exit()

    # change these to the working coherent ones
    # data_path = "/Users/mjcenpa/Data/C1"
    # workdir = "/Users/mjcenpa/Data/C1/Processed/built_runs"
    # cwd = os.getcwd()
    # os.chdir(workdir)


    # note: if this list ever gets too huge, you can switch to using os.walk
    raw_files = list(set(glob.glob("{}/**/**/Data/*Run*".format(data_path), recursive=True)))
    built_files = glob.glob("{}/Processed/built_runs/*.root".format(data_path))
    built_runs = list(set([int(f.split("OR_run")[-1].split(".")[0]) for f in built_files]))

    for i, file in enumerate(raw_files):

        # get the run number
        run = int(file.split("Run")[-1])

        if run not in built_runs:
            print("\nProcessing run", run, " ...")

            # use a quieter build command to limit log file size
            cmd1 = """majorcaroot -v error --sis --donotbuild {}""".format(raw_files[i])
            print("  ",cmd1)
            sh(cmd1)

    now = datetime.datetime.now()
    print("Processing is up to date!", now.strftime("%Y-%m-%d %H:%M"))
    os.chdir(cwd)


def sh(cmd, sh=False):
    """ Wraps a shell command."""
    import shlex
    import subprocess as sp
    if not sh: sp.call(shlex.split(cmd))  # "safe"
    else: sp.call(cmd, shell=sh)  # "less safe"


if __name__=="__main__":
    main(sys.argv[1:])
