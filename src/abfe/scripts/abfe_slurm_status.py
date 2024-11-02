#!/usr/bin/env python
# example script by snakemake.
import subprocess
import sys


def main():
    jobid = sys.argv[1]

    output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

    running_status = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
    if "COMPLETED" in output:
        print("success")
    elif any(r in output for r in running_status):
        print("running")
    else:
        print("failed")


if __name__ == "__main__":
    main()
