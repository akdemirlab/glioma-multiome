#BSUB -W 200:00
#BSUB -q long
#BSUB -n 2
#BSUB -M 32
#BSUB -R rusage[mem=32]
#BSUB -J $<tumor>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<tumor>.out

###############
#parameters
PRJ=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs
BAM_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/data/bwa_bam_new01
TUMOR=$<tumor>
NORMAL=$<normal>

TUM_BAM=$BAM_DIR/$TUMOR.bwa_sort.bam
NORM_BAM=$BAM_DIR/$NORMAL.bwa_sort.bam

VCF=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/human_9606_b151_GRCh38p7_vcf/00-common_all.vcf.gz
REF=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa

PILEUP_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results_new01/facets/snp-pileup
PILEUP_FILE=$PILEUP_DIR/$TUMOR.csv.gz
mkdir -p $PILEUP_DIR

FACET_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results_new01/facets/facets/$TUMOR
mkdir -p $FACET_DIR

cd $PRJ
###############
#snp-pileup
module load htslib/1.9
snp-pileup -g -q15 -Q20 -P100 -r25,0 $VCF $PILEUP_FILE $NORM_BAM $TUM_BAM

###############
#facets
cd $FACET_DIR
module load R/4.2.1
RSCRIPT=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/scripts/sv/facets.R
Rscript --vanilla $RSCRIPT $PILEUP_FILE $FACET_DIR $TUMOR
