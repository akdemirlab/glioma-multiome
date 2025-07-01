#BSUB -W 200:00
#BSUB -q long
#BSUB -n 4
#BSUB -M 32
#BSUB -R rusage[mem=32]
#BSUB -J $<tumor>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<tumor>.out

###############
#modules
module load smoove
module load singularity/3.7.0

###############
#parameters
TUMOR=$<tumor>
NORMAL=$<normal>
MOUNT1=/rsrch4
MOUNT2=/rsrch6
DATA_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/data/bwa_bam_new01
TUM_BAM=$DATA_DIR/$TUMOR.bwa_sort.bam
NORM_BAM=$DATA_DIR/$NORMAL.bwa_sort.bam
threads=4
RESULT_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results_new01/smoove/$TUMOR
reference_fasta=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa
BED=/rsrch4/home/neurosurgery/akdemirlab/references/smoove/exclude.cnvnator_100bp.GRCh38.20170403.bed

mkdir -p $RESULT_DIR
cd $RESULT_DIR
singularity run -B $MOUNT1:$MOUNT1 -B $MOUNT2:$MOUNT2 /risapps/singularity/repo/smoove/0.2.5/smoove_latest.sif smoove \
 call -x --name $TUMOR --exclude $BED --fasta $reference_fasta -p $threads \
 --genotype $TUM_BAM $NORM_BAM
