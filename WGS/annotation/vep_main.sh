#BSUB -W 240:00
#BSUB -q long
#BSUB -n 1 
#BSUB -M 8
#BSUB -R rusage[mem=8]
#BSUB -J vep
#BSUB -cwd /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/scripts/summary_v2/annotation
#BSUB -o /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/scripts/summary_v2/annotation/vep_main.out

module load nextflow
cd /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/scripts/summary_v2/annotation
nextflow run vep_main.nf