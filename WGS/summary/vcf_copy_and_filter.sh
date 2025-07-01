module load vcftools/0.1.16
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"
conda activate vcftools-0.1.16
module load htslib

function filter_snv() {
    INPUT_VCF=$1
    OUTPUT_VCF=$2
    if [ -e $OUTPUT_VCF ]
    then
        echo $OUTPUT_VCF exists
        exit
    fi
    vcftools --gzvcf $INPUT_VCF --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 \
    --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chr10 \
    --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 \
    --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 \
    --chr chr21 --chr chr22 --chr chrX --chr chrY \
    --remove-filtered-all --remove-indels --recode --stdout | bgzip -c > "$OUTPUT_VCF"
    tabix -p vcf  $OUTPUT_VCF
}

function filter_indel() {
    INPUT_VCF=$1
    OUTPUT_VCF=$2
    if [ -e $OUTPUT_VCF ]
    then
        echo $OUTPUT_VCF exists
        exit
    fi
    vcftools --gzvcf $INPUT_VCF --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 \
    --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chr10 \
    --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 \
    --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 \
    --chr chr21 --chr chr22 --chr chrX --chr chrY \
    --remove-filtered-all --keep-only-indels --recode --stdout | bgzip -c > "$OUTPUT_VCF"
    tabix -p vcf  $OUTPUT_VCF
}

function cp_sv() {
    INPUT_VCF=$1
    OUTPUT_VCF=$2
    if [ -e $OUTPUT_VCF ]
    then
        echo $OUTPUT_VCF exists
        exit
    fi
    cp $INPUT_VCF $OUTPUT_VCF
}

function filter_strelka() {
    SAMPLE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    filter_snv $INPUT_DIR/${SAMPLE}*strelka.somatic_snvs.vcf.gz $OUTPUT_DIR/$SAMPLE.snv.vcf.gz
    filter_indel $INPUT_DIR/${SAMPLE}*strelka.somatic_indels.vcf.gz $OUTPUT_DIR/$SAMPLE.indel.vcf.gz
}

function filter_mutect2() {
    SAMPLE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    filter_snv $INPUT_DIR/${SAMPLE}*.mutect2.filtered.vcf.gz $OUTPUT_DIR/$SAMPLE.snv.vcf.gz
    filter_indel $INPUT_DIR/${SAMPLE}*.mutect2.filtered.vcf.gz $OUTPUT_DIR/$SAMPLE.indel.vcf.gz
}

function cp_smoove(){
    SAMPLE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    cp_sv $INPUT_DIR/${SAMPLE}-smoove.genotyped.vcf.gz $OUTPUT_DIR/$SAMPLE.sv.vcf.gz
}

function cp_purple(){
    SAMPLE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    cp_sv $INPUT_DIR/$SAMPLE.purple.sv.vcf.gz $OUTPUT_DIR/$SAMPLE.sv.vcf.gz
}

function cp_svaba(){
    SAMPLE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    cp_sv $INPUT_DIR/somatic_run.svaba.somatic.sv.vcf $OUTPUT_DIR/$SAMPLE.sv.vcf
    gzip $OUTPUT_DIR/$SAMPLE.sv.vcf
}

function cp_jabba(){
    SAMPLE=$1
    INPUT_DIR=$2
    OUTPUT_DIR=$3
    cp_sv $INPUT_DIR/jabba.simple.vcf $OUTPUT_DIR/$SAMPLE.sv.vcf
    gzip $OUTPUT_DIR/$SAMPLE.sv.vcf
}