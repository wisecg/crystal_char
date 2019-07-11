#!/usr/bin/env python3
import sys, os, glob, json, time
import pathlib
import datetime
import argparse
import shlex
import subprocess as sp
import sys
from pprint import pprint
from shutil import move
import time
import datetime as datetime
import numpy as np

def main(argv):
    """
    Automatic processing of COHERENT NaI(Tl) crystal data.
    Uses "crysDB.json" as a database.  When a crystal is scanned,
    we input its serial number and run numbers to the database,
    and run this code to process the data and handle the output files.
    For usage help, run: $ python auto_process.py -h

    TODO:
    - run number checker (verify all expected runs exist & files are OK)
    """
    # -- load our run database and make it global --
    global crysDB
    with open("crysDB.json") as f:
        crysDB = json.load(f)

    # -- parse args --
    par = argparse.ArgumentParser(description="coherent crystal characterization suite")
    arg = par.add_argument
    arg("-c", "--crys", type=str, help="set crystal S/N")
    arg("-p", "--proc", type=str, help="process a crystal")
    arg("-t", "--temp", type=str, help='start temperature data taking')
    arg("-pt", "--printtemp", type=str, help='print current temperature')
    arg("-a", "--all", action="store_true", help="process all crystals in the DB")
    arg("-o", "--over", action="store_true", help="overwrite existing files")
    arg("-t", "--temp", type=str, help='start temperature data taking')
    arg("-z", "--zip", action="store_true", help='run gzip on raw files (on cenpa-rocks)')
    arg("-s", "--sync", action="store_true", help='sync DAQ with cenpa-rocks')
    args = vars(par.parse_args())

    # -- set parameters --
    crys_sn, overwrite = None, False

    if args["crys"]:
        crys_sn = args["crys"]

    if args["over"]:
        overwrite = args["over"]

    # -- run analysis --
    if args["proc"]:
        sn = args["proc"]
        process_crystal(sn, overwrite)

    if args["all"]:
        all_sns = [k for k in crysDB if "SN" in k]
        for sn in all_sns:
            process_crystal(sn, overwrite)

    if args["sync"]:
        sync_data()

    if args["zip"]:
        # clean_gzip()
        zip_data(overwrite)

    if args["temp"]:
        run_num = args["temp"]
        measure_temp(run_num)

    if args["printtemp"]:
        print_temp()

def process_crystal(sn, overwrite=False):
    """
    given a serial number, create the directories for Calibration,
    run getSpectrum (aka orcaroot), and move the files
    """
    print("Processing crystal:", sn)

    if sn not in crysDB.keys():
        print("Crystal {} not found! Exiting ...".format(sn))
        exit()

    # get a list of all raw files on this computer
    fstr = crysDB["raw_path"]
    raw_files = list(set(glob.glob("{}/**/**/Data/*Run*".format(fstr), recursive=True)))

    # get a list of all processed (built) run numbers on this computer
    built_files = list(set(glob.glob("{}/{}/*".format(crysDB["built_path"], sn), recursive=True)))
    built_runs = []
    for f in built_files:
        if "OR_run" not in f:
            continue
        run = int(f.split("OR_run")[-1].split(".")[0])
        built_runs.append(run)
        print("found processed run:", run) # debug, delete this

    # get a list of run numbers for this crystal from our JSON file
    crys_runs = []
    for run_type in crysDB[sn]:
        crys_runs.extend(crysDB[sn][run_type])

    test_vals = {"Position":crysDB["pos_vals"],
                 "Voltage":crysDB["HV_vals"]}

    print("Run list:", crys_runs, "\nTest vals:")
    pprint(test_vals)

    # -- create the directory of folders needed by the Calibration code. --
    # if the directory already exists, preserve any files already there.
    for pos in range(1,6):
        path = "{}/{}/Position/Position_{}".format(crysDB["built_path"], sn, pos)
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)

    for vol in [600, 700, 800, 900, 1000]:
        path = "{}/{}/Voltage/{}_V".format(crysDB["built_path"], sn, vol)
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
        for run_type in crysDB[sn]:

            dest_folder = crysDB["built_path"]

            # get the index value
            idxs = [i for i, x in enumerate(crysDB[sn][run_type]) if x==run]
            if len(idxs) < 1:
                continue
            idx = idxs[0]

            test_val = test_vals[run_type][idx]

            if run_type == "Position":
                folder_name = "Position_{}".format(test_val)
            elif run_type == "Voltage":
                folder_name = "{}_V".format(test_val)

            cmd = "mv {} {}/{}/{}/{}/{}".format(out_file,
                  crysDB["built_path"], sn, run_type,folder_name, out_file)

            print(cmd)
            sh(cmd)

    # add a last check that we have all files we expect
    print("Listing output files:")
    sh("find {}/{}".format(crysDB["built_path"], sn))

    # TODO: if the sync option is set, upload to cenpa-rocks and delete the raw files

    now = datetime.datetime.now()
    print("Processing is up to date!", now.strftime("%Y-%m-%d %H:%M"))


def sync_data():
    """
    to run: `python auto_process.py -s`
    rsync between DAQ machine and cenpa-rocks
    """
    if "DAQ" not in os.environ:
        print("Error, we're not on the CENPA DAQ machine.  Exiting...")
        exit()

    print("Syncing data ...")

    # raw data
    raw_path = crysDB["raw_path"].replace(" ", "\ ")
    raw_loc = "{}/".format(raw_path)
    raw_rocks = "{}:{}/".format(crysDB["rocks_login"], crysDB["rocks_data2"])

    # run rsync for raw data (can take a while ...)
    cmd = "rsync -av {} {}".format(raw_loc, raw_rocks)
    print("Syncing raw data directories ...\n",cmd)
    sh(cmd)

    # -- built data
    built_loc = "{}/".format(crysDB["built_path"])
    built_rocks = "{}:{}/".format(crysDB["rocks_login"], crysDB["rocks_built"])

    # run rsync for built data
    cmd = "rsync -av {} {}".format(built_loc, built_rocks)
    print("Syncing built data directories ...\n",cmd)
    sh(cmd)

    # make sure every (raw and built) local file is found on cenpa-rocks

    # local list
    r_list = glob.glob(crysDB["raw_path"] + "/**", recursive=True)
    b_list = glob.glob(crysDB["built_path"] + "/**", recursive=True)
    f_list = r_list + b_list

    # remote list
    ls = sp.Popen(['ssh', crysDB["rocks_login"], 'ls -R {}'.format(crysDB["rocks_data2"])],
                  stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = ls.communicate()
    out = out.decode('utf-8')
    remote_list = out.split("\n")
    remote_list = [r for r in remote_list if ":" not in r and len(r)!=0]

    # make sure all files have successfully transferred
    for f in f_list:
        fname = f.split("/")[-1]
        if len(fname) == 0:
            continue
        if fname not in remote_list:
            print("whoa, ", fname, "not found in remote list!")
            exit()

    print("All files in:\n    {}\n    {}\nhave been backed up to cenpa-rocks."
          .format(crysDB["raw_path"], crysDB["built_path"]))
    print("It should be OK to delete local files.")

    # don't delete these files, orca needs them
    ignore_list = [".Orca", "RunNumber"]

    # now delete old files, ask for Y/N confirmation
    print("OK to delete local files? [Y/N]")
    if input() == "Y":
        for f in f_list:
            f.replace(" ", "\ ")
            if os.path.isfile(f):
                if any(ig in f for ig in ignore_list):
                    continue
                os.remove(f)

    now = datetime.datetime.now()
    print("Processing is up to date!", now.strftime("%Y-%m-%d %H:%M"))


def zip_data(overwrite=False):
    """
    to run: `python auto_process.py -z`
    NOTE: this is the part of auto_process that is executed via cron:
    $ [login as coherent to cenpa-rocks]
    $ crontab -e
    * * * * * ~/analysis/crystal_char/task.sh >> ~/analysis/crystal_char/logs/cron.log 2>&1
    (then change the *'s to be an appropriate time interval, say 4 hours)
    """
    raw_rocks = "{}/".format(crysDB["rocks_data2"])

    print("Getting file list ...")
    all_files = glob.glob(raw_rocks + "/**", recursive=True)

    print("Scanning files ...")
    for f_name in all_files:

        # if it has "Run", no file extension, and ends in a number,
        # it's an ORCA raw file.  brilliant
        f = f_name.split("/")[-1]
        if not("Run" in f and "." not in f and "RunNumber" not in f):
            continue
        raw_file = f_name

        # grab the run number
        run_num = int(f.split("Run")[-1])

        # check for already existing .tar.gz files for this run
        t_idx = [i for i,x in enumerate(all_files) if "Run{}.tar.gz".format(run_num) in x]

        # if a .tar.gz file exists for this run, check its integrity
        # before deleting the raw file
        for idx in t_idx:
            tar_file = all_files[idx]
            print("Found tar file and raw file for run", run_num, "running gunzip")
            p = sp.Popen(["gunzip","-t",tar_file], stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = p.communicate()
            if err is None:
                print("Zipped file passes checks. Deleting raw file:\n    ", raw_file)
                # os.remove(raw_file)
            else:
                print("Zipped file failed checks.  Deleting zipped file:\n    ", tar_file)
                # os.remove(tar_file)

        # compress the file
        this_dir = os.getcwd()
        dest_dir = "/".join(f_name.split("/")[:-1])
        os.chdir(dest_dir)
        cmd = "tar cvzf {}.tar.gz {} --remove-files".format(f, f)
        print("Zipping file:", cmd)
        sh(cmd)
        os.chdir(this_dir)

        # exit()


def clean_gzip():
    """
    oops, i gzipped when i should have tar'd. clean it up
    """
    this_dir = os.getcwd()
    os.chdir("/data/COHERENT2/data/CrystalChar/raw")
    all_files = glob.glob("./**", recursive=True)
    for f in all_files:
        if ".gz" in f and "tar" not in f:
            print(f)
            sh("gunzip " + f)
    os.chdir(this_dir)


def sh(cmd, sh=False):
    """ Wraps a shell command."""
    import shlex
    import subprocess as sp
    if not sh: sp.call(shlex.split(cmd))  # "safe"
    else: sp.call(cmd, shell=sh)  # "less safe"


def measure_temp(run_num):
    runinfo = open("runinfo.txt", "w")
    runinfo.write(run_num)
    runinfo.close()
    sp.Popen("scp runinfo.txt pi@10.66.193.80:~/tempstudy", shell=True)
    sp.Popen("ssh pi@10.66.193.80 'cd ~/tempstudy && source env/bin/activate && python3 characterization_run.py'", shell=True)


def print_temp():
    sp.Popen("ssh pi@10.66.193.80 'cd ~/tempstudy && source env/bin/activate && python3 print_temp.py'", shell=True)


if __name__=="__main__":
    main(sys.argv[1:])
