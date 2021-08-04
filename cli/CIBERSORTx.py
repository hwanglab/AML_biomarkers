#!/usr/bin/env python3

import subprocess
import argparse
import os
import shutil
import logging
import sys
import uuid
from datetime import datetime

parser = argparse.ArgumentParser(description = "Deconvolute Samples")

parser.add_argument("--dir", "-d", help="path to run directory", default=None)
parser.add_argument("--id", "-i", help="ID to use for outputs", required=True)
parser.add_argument("--cores", "-c", help="number of cores to use", default=1)
parser.add_argument("--verbose", "-v", help="verbose level to use [1 (DEBUG) - 5 (CRITICAL)]", default=2, type=int)
parser.add_argument("--mixture", "-m", help="what mixture file to use", required=True, nargs="+")
parser.add_argument("--batch-correct", help="Should B mode batch correction be applied?", action="store_true")
parser.add_argument("--debug-cibersort", "-D", help="Should the max amount of information from CIBERSORT be printed?", action="store_true")

argv = parser.parse_args()

logging.basicConfig(level=vars(argv).get("verbose"), format="%(levelname)s [%(asctime)s] %(name)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

def concatenate_list_data(list):
    result= ''
    for element in list:
        result += str(element)
        result += " "
    return result
if (len(vars(argv).get("mixture")) == 1):
    name_to_use = vars(argv).get("mixture")[0]
else:
    name_to_use = concatenate_list_data(vars(argv).get("mixture"))

logger = logging.getLogger("CIBERSORTx: [{}]".format(name_to_use))

docker_singularity_cmd = "docker"
try: 
    subprocess.run("docker", stderr=subprocess.DEVNULL)
except FileNotFoundError as e:
    logger.info("Docker not Found. Falling back to Singularity")
    docker_singularity_cmd = "singularity"

try: 
    subprocess.run("singularity", stderr=subprocess.DEVNULL)
except FileNotFoundError as e:
    logger.critical("Singularity Not Found")
    raise RuntimeError("Docker or Singularity cannot be located")

if docker_singularity_cmd == "singularity":
    use_singularity = True
    logger.info("Using Singularity")
elif docker_singularity_cmd == "docker":
    use_singularity = False
    logger.info("Using Docker")
else:
    raise SyntaxError("Only Singularity and Docker are supported")

# set up outputs
if vars(argv).get("dir") == None:
    output_path = "outs/" + vars(argv).get("id")
    plots_path = "plots/" + vars(argv).get("id")
else:
    output_path = vars(argv).get("dir") + "/outs/" + vars(argv).get("id")
    plots_path = vars(argv).get("dir") + "/plots/" + vars(argv).get("id")

random_string = uuid.uuid1()
name_of_output_directory = "/cibersort_results/{}".format(random_string)

logger.debug("Output Directory:{}".format(name_of_output_directory))

try:
    os.mkdir(output_path + "/cibersort_results")
except FileExistsError as e:
    logger.debug("Output directory already created")

try: 
    os.mkdir(output_path + name_of_output_directory)
except FileExistsError as e:
    files = os.listdir(output_path + name_of_output_directory)
    for f in files:
        if f.startswith("temp"):
            logger.warning("Temporary files already in temporary directory. Two CIBERSORTx instances may share an run directory.")
    logger.debug("Temporary directory already exists")

# set up bind mounts
if docker_singularity_cmd == "docker":
    path_to_run_folder = os.path.abspath(os.getcwd())
    output_path = "{}/{}".format(path_to_run_folder, output_path)
input_bind = output_path + ":/src/data"
output_bind = output_path + name_of_output_directory + ":/src/outdir"
   
user_email = os.getenv("EMAIL")
user_token = os.getenv("TOKEN")

if user_email is None or user_token is None:
    logger.critical("Enviorment Variables not found!")
    raise RuntimeError("Env vars not found!")

logger.info("Preparing Run Commands")
if use_singularity:
    container_cmd = ["singularity", "run", "--pwd", "/src", "-B", input_bind, "-B", output_bind, "docker://cibersortx/fractions"]
else:
    container_cmd = ["docker", "run", "-v", input_bind, "-v", output_bind, "cibersortx/fractions"]

logger.debug("Testing if Reference can be found")

ref_filename = "CIBERSORTx_cell_type_sourceGEP.txt"
cibersort_cmd = ["--username", user_email, "--verbose", str(vars(argv).get("debug_cibersort")).upper(), "--token", user_token, " --single_cell", "TRUE", "--outdir", output_path + name_of_output_directory]

if os.path.isfile("{}/{}".format(output_path, ref_filename)):
    logger.info("Existing Reference Found")
elif os.path.isfile("{}/cibersort_results/{}".format(output_path, ref_filename)):
    logger.info("Existing Reference Found")
    logger.info("Copying Reference to CIBERSORTx input directory")
    source_path = "{}/cibersort_results/{}".format(output_path, ref_filename)
    destination_path = "{}/{}".format(output_path, ref_filename)
    shutil.copyfile(source_path, destination_path)
else:
    logger.info("Reference Not Found. Please run CIBERSORTx_ref.py")
    raise RuntimeError("Reference Not Found")
    
def run_cmd(mixture):
    specific_args = ["--sigmatrix", ref_filename, "--mixture", mixture, "--label", mixture, "--rmbatchBmode", str(vars(argv).get("batch_correct")).upper()]
    res = [container_cmd, cibersort_cmd, specific_args]
    flattened = [val for sublist in res for val in sublist]
    return flattened

data = vars(argv).get("mixture")

for dat in data:
    logger.info("Running CIBERSORTx on {}".format(dat))
    source_path = "{}/{}".format("cibersort_in", dat)
    destination_path = "{}/{}".format(output_path, dat)
    shutil.copyfile(source_path, destination_path)
    logger.debug("The run command is: {}".format(" ".join(run_cmd(dat))))
    if vars(argv).get("debug_cibersort"):
        subprocess_fun = None
    else:
        subprocess_fun = subprocess.DEVNULL
    start_time = datetime.now()
    subprocess.run(" ".join(run_cmd(dat)), shell=True, stdout=subprocess_fun)
    end_time = datetime.now()
    elapsed_time = end_time - start_time
    logger.info("CIBERSORTx completed in {}".format(elapsed_time))

files = os.listdir("{}/cibersort_results/{}".format(output_path, random_string))

for f in files:
    source_path = "{}/cibersort_results/{}/{}".format(output_path, random_string, f)
    destination_path = "{}/cibersort_results".format(output_path)
    shutil.copyfile(source_path, destination_path)

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
