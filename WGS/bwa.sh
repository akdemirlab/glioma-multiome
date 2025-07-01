#BSUB -W 200:00
#BSUB -q long
#BSUB -n 16 
#BSUB -M 64
#BSUB -R rusage[mem=64]
#BSUB -J $<Sequencing_ID>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<Sequencing_ID>.out

THREAD=16
DATA_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/data/fastq_new01
OUT_DIR=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/data/bwa_bam_new01

REF=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa
SAMPLE=$<Sequencing_ID>
FAST1=$DATA_DIR/Sample_${SAMPLE}/${SAMPLE}*R1_001.fastq.gz
FAST2=$DATA_DIR/Sample_${SAMPLE}/${SAMPLE}*R2_001.fastq.gz
BAM_38=$OUT_DIR/${SAMPLE}.bwa_sort.bam

if test -f $BAM_38; then
echo bam38 $BAM_38 already exists.
exit
fi

module load bwa/0.7.17
module load samtools/1.15
cd $OUT_DIR
bwa mem -t $THREAD -R "@RG\tID:1\tSM:$<Sequencing_ID>" $REF $FAST1 $FAST2 | samtools sort -@ $THREAD -o $BAM_38 -
samtools index $BAM_38