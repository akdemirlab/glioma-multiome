
MODEL=/rsrch9/home/genetics/vanloolab/tmitchell/lrs/dorado_0.9/models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0
MODEL_5mcs=/rsrch9/home/genetics/vanloolab/tmitchell/lrs/dorado_0.9/models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mCG_5hmCG@v3/
REF=/tmitchell/REF/hg38/GRCh38.primary_assembly.genome.fa

/rsrch9/home/genetics/vanloolab/tmitchell/lrs/dorado_0.9/dorado-0.9.1-linux-x64/bin/dorado basecaller $MODEL $FILEPATH --reference $REF --modified-bases-models $MODEL_5mcs -b 1500 > ${OUT_PATH}/${FILENAME}/${FILENAME}_aligned.bam

conda activate samtools-1.16.1

samtools sort -@ $THREAD ${OUT_PATH}/${FILENAME}/${FILENAME}_aligned.bam -o ${OUT_PATH}/${FILENAME}/${FILENAME}_sorted.bam
samtools index ${OUT_PATH}/${FILENAME}/${FILENAME}_sorted.bam

samtools calmd -@ $THREAD -b ${OUT_PATH}/${FILENAME}/${FILENAME}_sorted.bam $REF > ${OUT_PATH}/${FILENAME}/${FILENAME}_md.bam
samtools index ${OUT_PATH}/${FILENAME}/${FILENAME}_md.bam

mkdir ${OUT_PATH}/${FILENAME}/QC_${FILENAME}

conda activate nanoplot-1.40.2

NanoPlot -t $THREAD --bam ${OUT_PATH}/${FILENAME}/${FILENAME}_md.bam -o ${OUT_PATH}/${FILENAME}/QC_${FILENAME}
