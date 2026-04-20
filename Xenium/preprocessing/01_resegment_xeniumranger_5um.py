

import glob
import os, sys
import subprocess

DATA_DIR    = os.path.abspath("../../../data")
RESULTS_DIR = os.path.abspath("../../../results")
SCRIPTS_DIR = os.path.abspath("../../../scripts")


dataDir = f'{DATA_DIR}/raw/akdemirlab'
resDir = f'{RESULTS_DIR}/resegment/expan5um'
srcDir = f'{SCRIPTS_DIR}/resegment/src'
outDir = f'{RESULTS_DIR}/resegment/out'
errorDir = f'{RESULTS_DIR}/resegment/error'
exp_dis = 5
# bsub parameters

def generate_commands(inn, out, sName, src, exp):
    command = f"""
    module load xeniumranger/1.7.0.2
    mkdir -p {out} && cd {out}

    xeniumranger resegment \
        --xenium-bundle {inn} \
        --id {sName} \
        --expansion-distance {exp} \
        --resegment-nuclei false
    """

    with open(os.path.join(src, f'{sName}.sh'),'w') as f:
        f.write(command)


def main():

    for xenium_output in glob.glob(f'{dataDir}/202*'):
        # os.chdir(resDir)
        new_dir = os.path.basename(xenium_output)
        # if not os.path.exists(new_dir):
        #     os.mkdir(new_dir)
        # os.chdir(new_dir)

        outputs = glob.glob(f'{xenium_output}/output*')
        for regionDir in outputs:
            fileName = os.path.basename(regionDir)
            generate_commands(regionDir,os.path.join(resDir, new_dir), fileName, srcDir,exp_dis)

            # submit jobs to lsf
            queue     = 'medium'
            job_name  = fileName
            wall_time = '24:00'
            memory    = 64
            n_node    = 12

            job_sub = f"bsub \
                            -J {job_name} \
                            -o {outDir}/{job_name}.o.txt \
                            -e {errorDir}/{job_name}.e.txt \
                            -cwd {srcDir} \
                            -q {queue}\
                            -W {wall_time} \
                            -n {n_node} \
                            -M {memory} \
                            -R rusage[mem={memory}] \
                            /bin/bash {os.path.join(srcDir, f'{fileName}.sh')}"

            subprocess.call(job_sub, shell=True)

if __name__ =='__main__':
    main()

