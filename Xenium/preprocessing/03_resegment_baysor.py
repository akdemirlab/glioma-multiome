

import glob
import os, sys
import subprocess
import pandas as pd

DATA_DIR    = os.path.abspath("../../../data")
RESULTS_DIR = os.path.abspath("../../../results")
SCRIPTS_DIR = os.path.abspath("../../../scripts")

dataDir          = f'{DATA_DIR}/raw/akdemirlab'
resDir           = f'{RESULTS_DIR}/resegment/baysor'
srcDir           = f'{SCRIPTS_DIR}/resegment/baysor'
filtered_csv_dir = f'{RESULTS_DIR}/resegment/baysor/filtered_csv'
outDir           = f'{RESULTS_DIR}/resegment/out'
errorDir         = f'{RESULTS_DIR}/resegment/error'

BAYSOR_BIN = '/rsrch4/home/genomic_med/bzhao2/pipelines/baysor/bin'


# bsub parameters

def save_filtered_transcript(transcript_csv, sname):
    transcript = pd.read_csv(transcript_csv)
    transcript = transcript[
        (~transcript['feature_name'].str.startswith('BLANK')) &
        (~transcript['feature_name'].str.startswith('DeprecatedCodeword')) &
        (~transcript['feature_name'].str.startswith('NegControlCodeword')) &
        (~transcript['feature_name'].str.startswith('NegControlProbe')) &
        (transcript['qv'] >= 20)]
    transcript['cell_id_1'] = transcript['cell_id'].map(dict(zip(list(transcript['cell_id'].unique()), list(range(len(transcript['cell_id'].unique()))))))

    transcript.rename(columns={'x_location':'x', 'y_location':'y', 'z_location':'z'}).\
        to_csv(f'{filtered_csv_dir}/{sname}_filtered_qv20_transcripts.csv',index=False)


def generate_commands(inn, out, sName, src):
    command = f"""

mkdir -p {out} && cd {out}

# run baysor
export PATH=$PATH:{BAYSOR_BIN}
baysor run -x x -y y -z z -g feature_name -m 25 --save-polygons GeoJSON -p --prior-segmentation-confidence 0.5 -o {out}/{sName}_baysor_seg.csv {filtered_csv_dir}/{sName}_filtered_qv20_transcripts.csv :cell_id_1

#run xeniumranger
module load xeniumranger/1.7.0.2

xeniumranger import-segmentation --id={sName} --xenium-bundle={inn} --transcript-assignment={out}/{sName}_baysor_seg.csv --viz-polygons={out}/{sName}_baysor_seg_polygons.json --units=microns --localcores=32 --localmem=138

"""
    with open(os.path.join(src, f'{sName}.sh'),'w') as f:
        f.write(command)


def main():

    for xenium_output in glob.glob(f'{dataDir}/202*'):

        new_dir = os.path.basename(xenium_output)
        outputs = glob.glob(f'{xenium_output}/output*')

        for regionDir in outputs:
            fileName = os.path.basename(regionDir)
            save_filtered_transcript(f'{regionDir}/transcripts.csv.gz', fileName)
            generate_commands(regionDir, os.path.join(resDir, new_dir), fileName, srcDir)

            # submit jobs to lsf
            queue     = 'medium'
            job_name  = fileName
            wall_time = '24:00'
            memory    = 150
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

if __name__ == '__main__':
    main()
