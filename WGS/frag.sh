#BSUB -W 240:00
#BSUB -q long
#BSUB -n 1
#BSUB -M 32
#BSUB -R rusage[mem=32]
#BSUB -J $<Sequencing_ID>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<Sequencing_ID>.out

###############
#parameters
PRJ=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs
RUN=$<Sequencing_ID>
BAM=$PRJ/data/${RUN}.bwa_sort.bam
FRG_RESULT=$PRJ/results/frag/$RUN
cd $PRJ

###############
#frag
GCMAP_PATH=/rsrch4/home/neurosurgery/akdemirlab/references/jabba/gcmap_hg38
GCMAP_WINDOW=200

module load jabba/1.0
export PATH=${PATH}:$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))')
module load samtools/1.15
mkdir -p $FRG_RESULT

frag -b $BAM -d $GCMAP_PATH -w $GCMAP_WINDOW -o $FRG_RESULT

###############
#generate merged coverage
COV_FILE=$FRG_RESULT/cov.rds
COV_50K_FILE=$FRG_RESULT/cov_50k.rds
SEQINFO=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/seqinfo_hg38.rds
RS=/rsrch4/home/neurosurgery/akdemirlab/scripts/jabba_scripts/frag_bin_merge.R

module load R/4.2.1
Rscript --vanilla $RS $SEQINFO $COV_FILE 50000 $COV_50K_FILE
