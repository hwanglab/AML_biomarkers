#!/usr/bin/env python3

import subprocess
import argparse
import os
import shutil
import logging
import sys
import uuid
from datetime import datetime
import lib

parser = argparse.ArgumentParser(description = "Deconvolute Samples", add_help=False)

runtime = parser.add_argument_group("Runtime Arguments")
runtime.add_argument("-d", "--dir", help="path to run directory", default=None)
runtime.add_argument("-i", "--id", help="ID to use for outputs", required=True)
runtime.add_argument("-m", "--mixture", help="what mixture file to use", required=True, nargs="+")

flags = parser.add_argument_group("optional argumuments")
flags.add_argument("-h", "--help", help = "show this help message and exit", action = "help")
flags.add_argument("-v", "--verbose", help="should debug messages be printed", action="count", default=0)
flags.add_argument("-D", "--debug", help="Should stdout be printed from CIBERSORTx", action="store_true")
flags.add_argument("--batch-correct", help="Should the integrated data from seurat be used?", action="store_true")
flags.add_argument("--batch-correct-method", "-M", help="which method from Seurat should be used?", choices=["CCA", "RPCA"])

group = flags.add_mutually_exclusive_group()
group.add_argument("-X", "--B-mode", help="Should B mode batch correction be applied?", action="store_true")
group.add_argument("-S", "--S-mode", help="Should S mode batch correction be applied?", action="store_true")

argv = parser.parse_args()
verbose_level = 2 - argv.verbose
logging.basicConfig(level=verbose_level, format="%(levelname)s [%(asctime)s] %(name)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

if (len(argv.mixture) == 1):
    name_to_use = argv.mixture[0]
else:
    name_to_use = lib.concatenate_list_data(argv.mixture)

logger = logging.getLogger("CIBERSORTx: [{}]".format(name_to_use))

if (argv.batch_correct): 
    if (argv.batch_correct_method is None):
        logger.critical("Batch correction method not set")
        sys.exit(1)

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
        sys.exit(1)

if docker_singularity_cmd == "singularity":
    use_singularity = True
    logger.info("Using Singularity")
elif docker_singularity_cmd == "docker":
    use_singularity = False
    logger.info("Using Docker")
else:
    logger.critical("Only Singularity and Docker are supported")
    sys.exit(1)

# set up outputs
if argv.dir == None:
    output_path = "outs/" + argv.id
    plots_path = "plots/" + argv.id
else:
    output_path = "{}/outs/{}".format(argv.dir, argv.id)
    plots_path = "{}/plots/{}".format(argv.dir, argv.id)

random_string = uuid.uuid1()
cibersort_output = lib.setup_results_dir(argv)
name_of_output_directory = "{}/{}".format(cibersort_output, random_string)

logger.debug("Output Directory: {}".format(name_of_output_directory))

try:
    os.mkdir("{}/{}".format(output_path, cibersort_output))
except FileExistsError as e:
    logger.debug("Output directory already created")

try: 
    os.mkdir("{}/{}".format(output_path, name_of_output_directory))
except FileExistsError as e:
    files = os.listdir("{}/{}".format(output_path, name_of_output_directory))
    for f in files:
        if f.startswith("temp"):
            logger.warning("Temporary files already in temporary directory. Two CIBERSORTx instances may share an run directory.")
    logger.debug("Temporary directory already exists")

# set up bind mounts
if docker_singularity_cmd == "docker":
    path_to_run_folder = os.path.abspath(os.getcwd())
    output_path = "{}/{}".format(path_to_run_folder, output_path)
input_bind = output_path + ":/src/data"
output_bind = "{}/{}:/src/outdir".format(output_path, name_of_output_directory)
   
user_email = os.getenv("EMAIL")
user_token = os.getenv("TOKEN")

if user_email is None or user_token is None:
    logger.critical("Enviorment Variables not found!")
    logger.debug("User Email: {}".format(user_email))
    logger.debug("User Token: {}".format(user_token))
    raise RuntimeError("Env vars not found!")

logger.info("Preparing Run Commands")
if use_singularity:
    container_cmd = ["singularity", "run", "--pwd", "/src", "-B", input_bind, "-B", output_bind, "docker://cibersortx/fractions"]
else:
    container_cmd = ["docker", "run", "-v", input_bind, "-v", output_bind, "cibersortx/fractions"]

logger.debug("Testing if Reference can be found")

bc_ext = "no_bc"
if argv.batch_correct:
    bc_ext = argv.batch_correct_method

ref_filename = "CIBERSORTx_cell_type_sourceGEP.txt"
cibersort_cmd = ["--username", user_email, "--verbose", str(argv.debug).upper(), "--token", user_token, " --single_cell", "TRUE", "--outdir", "{}/{}".format(output_path, name_of_output_directory)]

if os.path.isfile("{}/{}".format(output_path, ref_filename)):
    logger.info("Existing Reference Found")
elif os.path.isfile("{}/{}_cibersort_results/{}".format(output_path, bc_ext, ref_filename)):
    logger.info("Existing Reference Found")
    logger.info("Copying Reference to CIBERSORTx input directory")
    source_path = "{}/{}_cibersort_results/{}".format(output_path, bc_ext, ref_filename)
    destination_path = "{}/{}".format(output_path, ref_filename)
    shutil.copyfile(source_path, destination_path)
else:
    logger.critical("Reference Not Found. Please run CIBERSORTx_ref.py")
    sys.exit(1)
    
def run_cmd(mixture):
    specific_args = ["--sigmatrix", ref_filename, "--mixture", mixture, "--label", mixture]
    bc_args = ["--rmbatchBmode", str(argv.B_mode).upper(), "--rmbatchSmode", str(argv.S_mode).upper(), "--refsample", "cibersort_ref_input.txt"]
    res = [container_cmd, cibersort_cmd, specific_args, bc_args]
    flattened = [val for sublist in res for val in sublist]
    return flattened

for dat in argv.mixture:
    logger.info("Running CIBERSORTx on {}".format(dat))
    source_path = "{}/{}".format("cibersort_in", dat)
    destination_path = "{}/{}".format(output_path, dat)
    shutil.copyfile(source_path, destination_path)
    logger.debug("The run command is: {}".format(" ".join(run_cmd(dat))))
    if argv.debug:
        subprocess_fun = None
    else:
        subprocess_fun = subprocess.DEVNULL
    start_time = datetime.now()
    subprocess.run(" ".join(run_cmd(dat)), shell=True, stdout=subprocess_fun)
    end_time = datetime.now()
    elapsed_time = end_time - start_time
    logger.info("CIBERSORTx completed in {}".format(elapsed_time))

files = os.listdir("{}/{}/{}".format(output_path, name_of_output_directory, random_string))

remove = True
for f in files:
    source_path = "{}/{}/{}".format(output_path, name_of_output_directory, f)
    destination_path = "{}/{}/{}".format(output_path, cibersort_output, f)
    try: 
        shutil.copyfile(source_path, destination_path)
    except IsADirectoryError as e:
        logger.error("Error moving files. Temporary Directory will not be deleted!")
        remove = False

if remove:
    shutil.rmtree(source_path)
 
lib.write_invocation(argv, output_path)
