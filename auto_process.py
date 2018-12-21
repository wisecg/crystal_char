#!/usr/bin/env python3
import sys, os, glob, json, time
import pathlib
import datetime
from pprint import pprint
from shutil import move

def main(argv):
    """
    Automatic processing of COHERENT NaI(Tl) crystal data.
    Uses "runDB.json" as a database.  When a crystal is scanned,
    we input its serial number and run numbers to the database,
    and run this code to process the data and handle the output files.
    Usage:
    $ python auto_process.py [options]
        default (no arguments): process all runs on this computer
        -crys [serial no]: process data from a particular crystal
        -sync : upload data to cenpa-rocks and delete the raw ORCA files
        -ovr : overwrite ROOT output files
    """
    # -- load our run database and make it global --
    global runDB
    with open("runDB.json") as f:
        runDB = json.load(f)

    overwrite = False

    # -- parse command line args --
    for i, opt in enumerate(argv):

        if opt=="-ovr":
            overwrite = True

        if opt=="-crys":
            sn = argv[i+1]
            process_crystal(sn, overwrite)
            exit()

    # default: run all crystals we haven't already processed
    sns = [id for id in runDB.keys() if "UW" in id]
    for sn in sns:
        process_crystal(sn, overwrite)


def process_crystal(sn, overwrite=False):
    """
    given a serial number, create the directories for Calibration,
    run getSpectrum (aka orcaroot), and move the files
    """
    print("Processing crystal:", sn)

    if sn not in runDB.keys():
        print("Crystal {} not found! Exiting ...".format(sn))
        exit()

    # get a list of all raw files on this computer
    fstr = runDB["raw_path"]
    raw_files = list(set(glob.glob("{}/**/**/Data/*Run*".format(fstr), recursive=True)))

    # get a list of all processed (built) run numbers on this computer
    built_files = list(set(glob.glob("{}/{}/*".format(runDB["built_path"], sn), recursive=True)))
    built_runs = []
    for f in built_files:
        if "OR_run" not in f:
            continue
        run = int(f.split("OR_run")[-1].split(".")[0])
        built_runs.append(run)
        print("found processed run:", run) # debug, delete this

    # get a list of run numbers for this crystal from our JSON file
    crys_runs = []
    for run_type in runDB[sn]:
        crys_runs.extend(runDB[sn][run_type])

    test_vals = {"Position":runDB["pos_vals"],
                 "Voltage":runDB["HV_vals"]}

    print("Run list:", crys_runs, "\nTest vals:")
    pprint(test_vals)

    # -- create the directory of folders needed by the Calibration code. --
    # if the directory already exists, preserve any files already there.
    for pos in range(1,6):
        path = "{}/{}/Position/position_{}".format(runDB["built_path"], sn, pos)
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    for vol in [600, 700, 800, 900, 1000]:
        path = "{}/{}/Voltage/{}_V".format(runDB["built_path"], sn, vol)
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    # -- loop over the raw files --
    # process only the runs for this crystal, and sort into built directories
    for i, f in enumerate(sorted(raw_files)):

        run = int(f.split("Run")[-1])
        if run not in crys_runs:
            continue

        print("Processing run {}, started at: {}".format(run, datetime.datetime.now()))

        out_file = "NaI_ET_run{}.root".format(run)

        if os.path.isfile(out_file) and overwrite is False:
            print("built file already exists, use -ovr to overwrite")
            continue

        # -- actually process the ORCA file and create a ROOT one --
        # this step can take 8-10 minutes
        t_start = time.time()
        ftmp = f.replace(' ', '\ ')
        cmd = "./getSpectrum --verbosity error {}".format(ftmp)
        print(cmd)
        sh(cmd)
        print("Done processing: {:.2f} min".format((time.time()-t_start)/60))

        # check output file
        if not os.path.isfile(out_file):
            print("Error, I expected to find an output file:", out_file)
            exit()

        # # this check doesn't work, but could be fixed
        # from ROOT import TFile
        # tf = TFile(out_file)
        # if tf.IsZombie():
        #     print("CAAARRL, it's a goddamn Walker! ... file."
        #           "    Found corrupted file: {}".format(out_file))
        #     print("    Processing on this file was probably interrupted.")
        #     exit()

        # now figure out which folder to move it to
        for run_type in runDB[sn]:

            dest_folder = runDB["built_path"]

            # get the index value
            idxs = [i for i, x in enumerate(runDB[sn][run_type]) if x==run]
            if len(idxs) < 1:
                continue
            idx = idxs[0]

            test_val = test_vals[run_type][idx]

            if run_type == "Position":
                folder_name = "position_{}".format(test_val)
            elif run_type == "Voltage":
                folder_name = "{}_V".format(test_val)

            cmd = "mv {} {}/{}/{}/{}/{}".format(out_file,
                  runDB["built_path"],sn, run_type,folder_name, out_file)

            print(cmd)
            sh(cmd)

    # add a last check that we have all files we expect
    # if so, let's do the sync option and upload to cenpa-rocks

    now = datetime.datetime.now()
    print("Processing is up to date!", now.strftime("%Y-%m-%d %H:%M"))
    # os.chdir(cwd)


def sh(cmd, sh=False):
    """ Wraps a shell command."""
    import shlex
    import subprocess as sp
    if not sh: sp.call(shlex.split(cmd))  # "safe"
    else: sp.call(cmd, shell=sh)  # "less safe"


if __name__=="__main__":
    main(sys.argv[1:])
