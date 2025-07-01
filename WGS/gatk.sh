#BSUB -W 240:00
#BSUB -q long
#BSUB -n 1 
#BSUB -M 8
#BSUB -R rusage[mem=8]
#BSUB -J gatk_b1
#BSUB -cwd /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/scripts
#BSUB -o /rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/scripts/gatk_b1.out

## ==================================
## Space needed in GATK preprocessing
##
## Input
## bam: 93G+81G+100G+66G=340G
##
## Output
## markduplicates: 231G
## recalibrated: 265G
## work: 2.9T
## -----

module load nextflow/22.10.1
module load singularity/3.7.0
export NXF_SINGULARITY_CACHEDIR=/rsrch4/home/neurosurgery/akdemirlab/pipelines/nf-core-images
NXF_OPTS='-Xms1g -Xmx4g'
IGENOMES=/rsrch4/home/neurosurgery/akdemirlab/references/igenomes
WORKFLOW=/rsrch4/home/neurosurgery/akdemirlab/pipelines/nf-core-sarek-3.1.2/workflow
CONFIG=/rsrch4/home/neurosurgery/akdemirlab/pipelines/nf-core_sarek_seadragon.config
SAMPLE_SHEET=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/scripts/gatk_samplesheet_b1.csv
GATK_RESULT=/rsrch4/home/moonshots/WorkingFolder/catalyst_wgs/results/gatk_b1

mkdir -p $GATK_RESULT
cd $GATK_RESULT
 
nextflow run $WORKFLOW \
 --igenomes_base $IGENOMES \
 --genome GATK.GRCh38 \
 --split_fastq 0 \
 -profile singularity \
 --save_output_as_bam \
 -c $CONFIG \
 --step mapping \
 --skip_tools fastqc,multiqc \
 --tools mutect2,strelka,ascat,cnvkit \
 --input $SAMPLE_SHEET \
 --outdir $GATK_RESULT