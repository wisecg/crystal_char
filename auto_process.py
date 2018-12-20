#!/usr/bin/env python3
import sys, os, glob, datetime

def main():

    # change these to the working coherent ones
    data_path = "/Users/mjcenpa/Data/C1"
    workdir = "/Users/mjcenpa/Data/C1/Processed/built_runs"
    cwd = os.getcwd()
    # os.chdir(workdir)

    print("hi clint")
    print("data_path:",data_path)
    print("workdir:",workdir)

    exit()

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
    main()
