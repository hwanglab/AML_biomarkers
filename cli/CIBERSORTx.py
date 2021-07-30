#!/usr/bin/env python3

import subprocess
import argparse
import os
import shutil

parser = argparse.ArgumentParser(description = "Deconvolute Samples")

parser.add_argument("--dir", "-d", help="path to run directory", default=None)
parser.add_argument("--ref", "-r", help="path to CIBERSORTx GEP", required=True)
parser.add_argument("--id", "-i", help="ID to use for outputs", required=True)
parser.add_argument("--cores", "-c", help="number of cores to use", default=1)
parser.add_argument("--verbose", "-v", help="verbose level to use", default=2)
parser.add_argument("--container-software", "-C", help="what container software to use", default="singularity")

argv = parser.parse_args(["-r", "CIBERSORTx_cibersort_ref_input_inferred_refsample.txt", "-i", "flt3_cebpa"])

docker_singularity_cmd = vars(argv).get("container_software")

if vars(argv).get("container_software") == "singularity":
    use_singularity = True
elif vars(argv).get("container_software") == "docker":
    use_singularity = False
else:
    raise SyntaxError("Only Singularity and Docker are supported")

# set up outputs
if vars(argv).get("dir") == None:
    output_path = "outs/" + vars(argv).get("id")
    plots_path = "plots/" + vars(argv).get("id")
else:
    output_path = vars(argv).get("dir") + "/outs/" + vars(argv).get("id")
    plots_path = vars(argv).get("dir") + "/plots/" + vars(argv).get("id")

name_of_output_directory = "/cibersort_results"

# set up bind mounts
input_bind = output_path + ":/src/data"
output_bind = output_path + name_of_output_directory + ":/src/outdir"

# Code to make reference here

user_email = os.getenv("EMAIL")
user_token = os.getenv("TOKEN")

if use_singularity:
    container_cmd = ["singularity", "run", "--pwd", "/src", "-B", input_bind, "-B", output_bind, "docker://cibersortx/fractions"]
else:
    container_cmd = "docker" + " " + "run" + "-v " + input_bind + "-v " + output_bind + " cibersortx/fractions"

cibersort_cmd = ["--username", user_email, "--verbose", "FALSE", "--token", user_token, " --single_cell", "TRUE", "--outdir", output_path + name_of_output_directory]

ref_filename = vars(argv).get("ref")

def run_cmd(mixture):
    specific_args = ["--sigmatrix", ref_filename, "--mixture", mixture]
    res = [container_cmd, cibersort_cmd, specific_args]
    flattened = [val for sublist in res for val in sublist]
    return flattened

data = ["tcga_data.txt", "target_data.txt", "beat_aml.txt"]

for dat in data:
    source_path = "{}/{}".format("cibersort_in", dat)
    destination_path = "{}/{}".format(output_path, dat)
    shutil.copyfile(source_path, destination_path)
    subprocess.run(run_cmd(dat))

def write_invocation(argv, output_path):
    import time
    import sys
    args = vars(argv).items()
    args = list(args)
    sys_cmd = sys.argv[0]

    sys_time = time.asctime(time.localtime(time.time()))

    header = "{}: {}------------------".format(sys_cmd, sys_time)

    fname = output_path + "/invocation"
    file1 = open(fname, "a")  # append mode
    file1.write(header)
    file1.write("\n")
    for arg in args:
        formatted_args = "{}: {}".format(arg[0], arg[1])
        file1.write(formatted_args)
        file1.write("\n")
    file1.write("\n")
    file1.close()
 
write_invocation(argv, output_path)
