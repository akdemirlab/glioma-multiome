#BSUB -W 200:00
#BSUB -q long
#BSUB -n 12
#BSUB -M 128
#BSUB -R rusage[mem=128]
#BSUB -J $<blood_wgs_seq_id>
#BSUB -cwd $<MSUB_OUTPUT_DIR>
#BSUB -o $<MSUB_OUTPUT_DIR>/$<blood_wgs_seq_id>.out

TUMOR="$<wgs_seq_id>"
NORMAL="$<blood_wgs_seq_id>"

. ~/.bashrc
conda activate hatchet
module load cbc
export HATCHET_COMPUTE_CN_SOLVER=cbc
module load shapeit/v2.904

PRJ=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs
INI_FILE=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/scripts/sv/hatchet.ini
read -a array <<< "$TUMOR"
TUM_BAM=""
for i in "${array[@]}"
do
  TUM_BAM="$TUM_BAM `echo $PRJ/data/bwa_bam*/$i.bwa_sort.bam | awk '{print $1}'`"
done
#TUM_BAM=`echo $PRJ/data/bwa_bam*/$TUMOR.bwa_sort.bam | awk '{print $1}'`
NORM_BAM=`echo $PRJ/data/bwa_bam*/$NORMAL.bwa_sort.bam | awk '{print $1}'`
CORES=10
REF=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa
OUTPUT=/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_wgs/results/hatchet_v3/$NORMAL

mkdir -p $OUTPUT
cd $OUTPUT
cp $INI_FILE $NORMAL.ini
echo reference = \"$REF\" >>$NORMAL.ini
echo normal = \"$NORM_BAM\" >>$NORMAL.ini
echo bams = \"$TUM_BAM\" >>$NORMAL.ini
echo samples = \"$TUMOR\" >>$NORMAL.ini
echo output = \"$OUTPUT\" >>$NORMAL.ini
echo processes = $CORES >>$NORMAL.ini
hatchet run $NORMAL.ini >$NORMAL.log
