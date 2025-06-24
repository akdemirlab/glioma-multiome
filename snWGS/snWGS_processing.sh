
#this script was run as a job array for  each snSGS cell fastq.

###################################                              ###################################
###################################    Process Raw scDNA Data    ###################################
###################################                              ###################################


REF=$HOME/ref/GRCh38.primary_assembly.genome.fa
DATA=***{data_path}***
OUT=***{out_path}***
THREAD=16

CELL_LIST=($(ls $DATA/))
CELL=${CELL_LIST[$((LSB_JOBINDEX))-1]}

echo $CELL_LIST
echo $CELL

#generate a new merged fastq
cat $DATA/$(basename $CELL)/*.gz > $DATA/$(basename $CELL)/merged.fastq.gz

#make cell out directory
mkdir $OUT/$(basename $CELL)    


###################################    BWA-MEM2 Align    ###################################

module load bwa-mem2
conda activate bwa-mem2-2.2.1

bwa-mem2 mem $REF $DATA/$(basename $CELL)/merged.fastq.gz > $OUT/$(basename $CELL)/aligned.sam

###################################    Samtools Sort    ###################################


module load samtools
conda activate samtools-1.16.1


samtools view -S -b $OUT/$(basename $CELL)/aligned.sam > $OUT/$(basename $CELL)/aligned.bam
samtools sort -@ $THREAD $OUT/$(basename $CELL)/aligned.bam -o $OUT/$(basename $CELL)/sorted.bam


###################################    Picard MarkDuplicates    ###################################

module load picard
conda activate picard-2.27.4

picard MarkDuplicates -I $OUT/$(basename $CELL)/sorted.bam -O $OUT/$(basename $CELL)/$(basename $CELL).bam -M $OUT/$(basename $CELL)/dup_metrics.txt


###################################    INDEX    ###################################


module load samtools
conda activate samtools-1.16.1


samtools index $OUT/$(basename $CELL)/$(basename $CELL).bam


###################################    Remove unnecessary files    ###################################

rm $OUT/$(basename $CELL)/aligned.sam
rm $OUT/$(basename $CELL)/aligned.bam
rm $OUT/$(basename $CELL)/sorted.bam
