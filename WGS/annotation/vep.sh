#BSUB -W 200:00
#BSUB -q long
#BSUB -n 2
#BSUB -M 16
#BSUB -R rusage[mem=16]
#BSUB -J $<tumor>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<tumor>.out

#use
TUMOR=$<tumor>
module load perl/5.28.1

cd /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results_new01/summary/summary_snv
mkdir -p vep
vep --cache --offline --everything -i snv_overlap/$TUMOR.vcf -o vep/$TUMOR.tsv