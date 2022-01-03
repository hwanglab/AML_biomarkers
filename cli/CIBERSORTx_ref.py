#!/usr/bin/env python3

import subprocess
import argparse
import os
import shutil
import logging
import sys
import binascii
import lib

parser = argparse.ArgumentParser(description = "Deconvolute Samples")

parser.add_argument("--dir", "-d", help="path to run directory", default=None)
parser.add_argument("--id", "-i", help="ID to use for outputs", required=True)
parser.add_argument("--cores", "-c", help="number of cores to use", default=1)
parser.add_argument("--verbose", "-v", help="should debug messages be printed", action="count", default=0)
parser.add_argument("--debug", "-D", help="Should the max amount of information from CIBERSORT be printed?", action="store_true")

argv = parser.parse_args()

verbose_level = 2 - argv.verbose
logging.basicConfig(level=verbose_level, format="%(levelname)s [%(asctime)s] %(name)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

logger = logging.getLogger("CIBERSORTx: Reference")

docker_singularity_cmd = "docker"
try: 
    subprocess.run("docker", stderr=subprocess.DEVNULL)
except FileNotFoundError as e:
    logger.info("Docker not Found. Falling back to Singularity")
    docker_singularity_cmd = "singularity"

if docker_singularity_cmd == "singularity":
    try: 
        subprocess.run("singularity", stderr=subprocess.DEVNULL)
    except FileNotFoundError as e:
        logger.critical("Singularity Not Found")
        raise RuntimeError("Docker or Singularity cannot be located")
    use_singularity = True
    logger.info("Using Singularity")
elif docker_singularity_cmd == "docker":
    use_singularity = False
    logger.info("Using Docker")
else:
    raise SyntaxError("Only Singularity and Docker are supported")

# set up outputs
if argv.dir == None:
    output_path = "outs/" + argv.id
    plots_path = "plots/" + argv.id
else:
    output_path = argv.dir + "/outs/" + argv.dir
    plots_path = argv.dir + "/plots/" + argv.dir

random_string = binascii.b2a_hex(os.urandom(15))
bc_ext = "no_bc"


name_of_output_directory = "{}_cibersort_results".format(bc_ext)

logger.debug("Output Directory: {}/{}".format(output_path, name_of_output_directory))

try:
    os.mkdir("{}/{}".format(output_path, name_of_output_directory))
except FileExistsError as e:
    logger.debug("Output directory already created")

# set up bind mounts
if docker_singularity_cmd == "docker":
    path_to_run_folder = os.path.abspath(os.getcwd())
    output_path = "{}/{}".format(path_to_run_folder, output_path)
input_bind = output_path + ":/src/data"
output_bind = "{}/{}:/src/outdir".format(output_path, name_of_output_directory)
   
user_email = os.getenv("EMAIL")
user_token = os.getenv("TOKEN")

logger.info("Preparing Run Commands")
if use_singularity:
    container_cmd = ["singularity", "run", "--pwd", "/src", "-B", input_bind, "-B", output_bind, "docker://cibersortx/fractions"]
else:
    container_cmd = ["docker", "run", "-v", input_bind, "-v", output_bind, "cibersortx/fractions"]

logger.debug("The container command is: {}".format(" ".join(container_cmd)))

logger.debug("Testing if Reference can be found")

ref_filename = "CIBERSORTx_cell_type_sourceGEP.txt"
ref_filename_end = "CIBERSORTx_cell_type_sourceGEP.txt"
cibersort_cmd = ["--username", user_email, "--verbose", "FALSE", "--token", user_token, " --single_cell", "TRUE", "--outdir", "{}/{}".format(output_path, name_of_output_directory)]
logger.debug("The CIBERSORTx command is: {}".format(" ".join(cibersort_cmd)))

if os.path.isfile("{}/{}".format(output_path, ref_filename_end)):
    logger.info("Existing Reference Found")
elif os.path.isfile("{}/cibersort_results/{}".format(output_path, ref_filename)):
    logger.info("Existing Reference Found")
    logger.info("Copying Reference to CIBERSORTx input directory")
    source_path = "{}/cibersort_results/{}".format(output_path, ref_filename)
    destination_path = "{}/{}".format(output_path, ref_filename)
    shutil.copyfile(source_path, destination_path)
else:
    logger.info("Reference Not Found. Creating Reference")
    if not os.path.isfile("{}/cibersort_ref_input.txt".format(output_path)):
        raise FileNotFoundError("CIBERSORTx Reference Inputs not Found -- please run make_CIBERSORT_ref.R first")
    ref_create_cmd = ["--refsample", "cibersort_ref_input.txt", "--fraction", "0.50"]
    logger.debug("Joining commands")
    ref_run_cmd = [container_cmd, cibersort_cmd, ref_create_cmd]
    logger.debug("Unnesting commands")
    ref_run_cmd = [val for sublist in ref_run_cmd for val in sublist]
    logger.debug("The run command is: {}".format(" ".join(ref_run_cmd)))
    if argv.debug:
        subprocess_fun = subprocess.STDOUT
    else:
        subprocess_fun = subprocess.DEVNULL
    subprocess.run(" ".join(ref_run_cmd), shell=True, stdout=subprocess_fun)
    logger.info("Reference Creation Finished")
    logger.info("Copying Reference to CIBERSORTx input directory")
    source_path = "{}/{}".format(output_path, ref_filename)
    destination_path = "{}/{}".format(output_path, ref_filename_end)
    shutil.copyfile(source_path, destination_path)
 
lib.write_invocation(argv, output_path)