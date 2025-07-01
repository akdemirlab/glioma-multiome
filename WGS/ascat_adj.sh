#BSUB -W 200:00
#BSUB -q long
#BSUB -n 10
#BSUB -M 64
#BSUB -R rusage[mem=64]
#BSUB -J $<tumor>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<tumor>.out

###############
#parameters
PRJ=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs
PRJ_MAIN=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main
TUMOR=$<tumor>
NORMAL=$<normal>
GENDER=$<gender>
PURITY=$<purity>
PLOIDY=$<ploidy>
TUM_BAM=$PRJ_MAIN/*/data/bwa_bam*/$TUMOR.bwa_sort.bam
NORM_BAM=$PRJ_MAIN/*/data/bwa_bam*/$NORMAL.bwa_sort.bam
CORES=10
REF=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa

cd $PRJ


###############
#ASCAT
MOUNT1=/rsrch4
MOUNT2=/rsrch6
SNPGC=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv
ASCAT_OUTPUT=$PRJ/results_adj/ascat/$TUMOR

cd $PRJ_MAIN
mkdir -p $ASCAT_OUTPUT
module load ascatNgs/4.5.0-sing
singularity exec  -B $MOUNT1:$MOUNT1 -B $MOUNT2:$MOUNT2 $ascatngs_sif ascat.pl \
   -r $REF \
   -t $TUM_BAM \
   -n $NORM_BAM \
   -sg $SNPGC \
   -pr WGS \
   -g $GENDER \
   -gc chrY \
   -pl ILLUMINA \
   -rs Human \
   -ra GRCh38 \
   -c $CORES \
   -purity $PURITY \
   -ploidy $PLOIDY \
   -o $ASCAT_OUTPUT