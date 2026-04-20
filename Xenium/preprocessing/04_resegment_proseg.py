import spatialdata
import spatialdata_io
import pandas as pd
import dill
import glob
import subprocess
import os

import argparse

DATA_DIR    = os.path.abspath("../../../data")
RESULTS_DIR = os.path.abspath("../../../results")
SCRIPTS_DIR = os.path.abspath("../../../scripts")

pkl_dir  = f'{RESULTS_DIR}/resegment/split_pkl/expan5um'
csv_dir  = f'{RESULTS_DIR}/resegment/proseg/prior_expan5/filtered_csv'
res_dir  = f'{RESULTS_DIR}/resegment/proseg/prior_expan5/proseg'
srcDir   = f'{SCRIPTS_DIR}/resegment/proseg/prior_expan5'
outDir   = f'{RESULTS_DIR}/resegment/out'
errorDir = f'{RESULTS_DIR}/resegment/error'


# ---------------------------------------------------------------------------
# Step 1: filter transcripts from dill files and write per-sample CSVs
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Filter transcripts from output of xenium')
parser.add_argument('inputdir',  type=str, help='Input directory containing .dill files')
parser.add_argument('outputdir', type=str, help='Output directory for filtered transcript CSVs')
args = parser.parse_args()

dills = glob.glob(f"{args.inputdir}/*.dill")

for dil in dills:
    with open(dil,'rb') as f:
        spdata = dill.load(f)
    sname = os.path.basename(dil).split('.')[0]
    transcripts = spdata.points['transcripts'].compute()
    transcripts.rename(columns={'x':'x_location','y':'y_location','z':'z_location'}, inplace=True)
    filter_transcripts = transcripts[
        (transcripts['qv'] >= 20) &
        (~transcripts["feature_name"].str.startswith("BLANK_")) &
        (~transcripts["feature_name"].str.startswith("NegControlProbe_")) &
        (~transcripts["feature_name"].str.startswith("DeprecatedCodeword_")) &
        (~transcripts["feature_name"].str.startswith("NegControlCodeword_"))
    ]
    filter_transcripts.to_csv(f'{args.outputdir}/{sname}.csv', index=False)


# ---------------------------------------------------------------------------
# Step 2: submit proseg jobs to LSF
# ---------------------------------------------------------------------------

## run proseg
## didn't debug!!

def bsub(
        job_name,
        src_file,
        queue    = 'e40medium',
        wall_time= '24:00',
        n_node   = 12,
        memory   = 300,
        ):

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
        /bin/bash -c 'mkdir -p {res_dir}/{job_name} && cd {res_dir}/{job_name}; proseg --xenium {src_file}'"
    subprocess.Popen(job_sub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


for csv in glob.glob(f'{csv_dir}/*.csv'):
    sname = os.path.basename(csv).split('.')[0]
    bsub(sname, csv)
