#!/usr/bin/env nextflow

params.meta = "${projectDir}/meta.csv"
params.annotation_path = "/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2/annotation"

workflow {
    chan1 = Channel.fromPath(params.meta)
        .splitCsv(header:true)
        .map { row-> row.wgs_seq_id }

    snv_overlap_result = file(params.annotation_path + "/vep_snv_overlap")
    mutect2_result = file(params.annotation_path + "/vep_mutect2")
    mutect2_unfiltered_result = file(params.annotation_path + "/vep_mutect2_unfiltered")

    snv_overlap_result.mkdirs()
    mutect2_result.mkdirs()
    mutect2_unfiltered_result.mkdirs()

    vep_snv_overlap(chan1, snv_overlap_result) | view()
    vep_mutect2(chan1, mutect2_result) | view()
    //vep_mutect2_unfiltered(chan1, mutect2_unfiltered_result)
}

process vep_snv_overlap {
    tag "${wgs_seq_id}s"
    executor 'lsf'
    cpus 2
    time '200.h'
    queue 'long'
    clusterOptions '-M 64 -R rusage[mem=64]'

    input:
    val wgs_seq_id
    file "result_dir" 

    output:
    stdout
    
    script:
    def snv_path = "/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2/summary_snv/snv_overlap/"+wgs_seq_id+".vcf"
    """
    export PATH=/rsrch4/home/neurosurgery/akdemirlab/pipelines/ensembl-vep/:\$PATH
    vep --cache --offline --force --everything -i ${snv_path} -o result_dir/${wgs_seq_id}.snv.tsv
    """

    stub:
    """
    #echo "running ${wgs_seq_id}"
    #head ${snv_path}
    """
}



process vep_mutect2 {
    tag "${wgs_seq_id}m"
    executor 'lsf'
    cpus 2
    time '200.h'
    queue 'long'
    clusterOptions '-M 64 -R rusage[mem=64]'

    input:
    val wgs_seq_id
    file "result_dir" 

    output:
    stdout
    
    script:
    def snv_path = "/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2/mutect2/"+wgs_seq_id+".snv.vcf.gz"
    """
    export PATH=/rsrch4/home/neurosurgery/akdemirlab/pipelines/ensembl-vep/:\$PATH
    vep --cache --offline --force --everything -i ${snv_path} -o result_dir/${wgs_seq_id}.snv.tsv
    """

    stub:
    """
    #echo "running ${wgs_seq_id}"
    #head ${snv_path}
    """
}
