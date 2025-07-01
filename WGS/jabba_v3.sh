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
TUMOR=$<tumor>
NORMAL=$<normal>
GENDER=$<gender>
TUM_BAM=$PRJ/data/bwa_bam_new01/$TUMOR.bwa_sort.bam
NORM_BAM=$PRJ/data/bwa_bam_new01/$NORMAL.bwa_sort.bam
CORES=10
REF=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa

cd $PRJ

###############
#svaba
SVABA_PATH=/rsrch4/home/genomic_med/lye1/bin/svaba/bin
SVABA_RESULT=$PRJ/results_new01/svaba/$TUMOR

module load jabba/1.0
export PATH=${PATH}:$SVABA_PATH
mkdir -p $SVABA_RESULT
cd $SVABA_RESULT
svaba run -t $TUM_BAM -n $NORM_BAM -p $CORES -a somatic_run -G $REF --override-reference-check
cd -

TUMOR_NF_VCF=$SVABA_RESULT/somatic_run.svaba.unfiltered.somatic.sv.vcf
TUMOR_NO_CHR_NF_VCF=$SVABA_RESULT/somatic_run.svaba.unfiltered.somatic.no_chr.sv.vcf
if test ! -f $TUMOR_NO_CHR_NF_VCF; then
awk '{gsub(/^chr/,""); print}' $TUMOR_NF_VCF | awk '{gsub(/ID=chr/,"ID="); print}' > $TUMOR_NO_CHR_NF_VCF
fi

###############
#frag
GCMAP_PATH=/rsrch4/home/neurosurgery/akdemirlab/references/jabba/gcmap_hg38
GCMAP_WINDOW=200
FRG_NORMAL=$PRJ/results_new01/frag/$NORMAL
FRG_TUMOR=$PRJ/results_new01/frag/$TUMOR

module load R/4.2.1
RS=/rsrch4/home/neurosurgery/akdemirlab/scripts/jabba_scripts/frag_tnratio.R

COV_NORMAL=$FRG_NORMAL/cov.rds
COV_TUMOR=$FRG_TUMOR/cov.rds
COV_RATIO=$FRG_TUMOR/cov_ratio.rds
Rscript --vanilla $RS $COV_TUMOR $COV_NORMAL $COV_RATIO

###############
#ASCAT
MOUNT=/rsrch4
SNPGC=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv
ASCAT_OUTPUT=$PRJ/results_new01/ascat/$TUMOR

cd $PRJ
mkdir -p $ASCAT_OUTPUT
module load ascatNgs/4.5.0-sing
singularity exec -B $MOUNT:$MOUNT $ascatngs_sif ascat.pl \
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
   -o $ASCAT_OUTPUT


ASCAT=$ASCAT_OUTPUT/${TUMOR}.samplestatistics.txt
PLOIDY=`awk '{if($1=="Ploidy"){print $2}}' $ASCAT`
PURITY=`awk '{if($1=="rho"){print $2}}' $ASCAT`
echo ascat ploidy: $PLOIDY
echo ascat purity: $PURITY

################
#jba
module load jabba/1.0
JABBA_PATH=$(Rscript -e 'cat(paste0(installed.packages()["JaBbA", "LibPath"], "/JaBbA/extdata/"))')
export PATH=${JABBA_PATH}:${PATH}

JABBA_RESULT=$PRJ/results_new01/jba_v3/$TUMOR
mkdir -p $JABBA_RESULT
jba $TUMOR_NO_CHR_NF_VCF $COV_RATIO -o $JABBA_RESULT\
 --cores $CORES --ploidy $PLOIDY --purity $PURITY --tilim 100000