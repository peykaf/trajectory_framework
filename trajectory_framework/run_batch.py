import os
import argparse
import subprocess

def create_slurm_script(filename, queue=None, mem=None, nodes=None):
    script_filename = "tempscript.sh"
    with open(script_filename, "w") as myscript:
        myscript.write("#!/bin/bash\n")
        myscript.write("#SBATCH -n 1\n")
        if mem and mem > 0:
            myscript.write("#SBATCH --mem=%i\n" % mem)
        if queue:
            myscript.write("#SBATCH -p %s\n" % queue)
        if nodes:
            myscript.write("#SBATCH -w, --nodelist=%s\n" % nodes)


        myscript.write("./trajectory_framework %s" % filename)

        return script_filename

def submit_job(filename, batch, queue=None, mem=None, nodes=None):
    if batch == "slurm":
        script_file = create_slurm_script(filename, queue, mem, nodes)
    else:
        raise Exception("This script does not support the %s batch system\n" % batch)

    subprocess.call("sbatch -J %s %s" % (filename, script_file), shell=True)
    os.remove(script_file)



default_batch = 'slurm'
default_queue = None
default_mem = None
default_nodes = None

parser = argparse.ArgumentParser(description='Submit optimisation jobs to batch queuing system')
parser.add_argument('filename', help='The input filename for optimisation', type=str)
parser.add_argument('-b', '--batch', default=default_batch, type=str, help='the batch queuing system to use on the server')
parser.add_argument('-q', '--queue', default=default_queue, type=str, help='the queue in which to submit the job')
parser.add_argument('-m', '--memory', default=default_mem, type=int, help='memory requested from the queuing system')
parser.add_argument('-n', '--nodes', default=default_nodes, type=str, help='nodes requested from the queuing system')

args = parser.parse_args()

submit_job(args.filename, args.batch, args.queue, args.memory, args.nodes)

