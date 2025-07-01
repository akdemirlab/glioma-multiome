#!/usr/bin/env nextflow

params.directories = "${projectDir}/sample_directories.csv"
params.summary = "/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2"
workflow {
    chan1 = Channel.fromPath(params.directories)
        .splitCsv(header:true)
        .map { row-> tuple(row.patient_id, row.sample_id, row.wgs_seq_id, file(row.strelka_dir), file(row.mutect2_dir), file(row.smoove_dir), file(row.purple_dir), file(row.svaba_dir), file(row.jabba_dir), file(row.ascat_dir), file(row.facets_dir)) }

    make_summary_dir(params.summary)
    summary_snv_sv_cnv(chan1, make_summary_dir.out) | view()
    copy_purple_driver(chan1, make_summary_dir.out) | view()
}

process make_summary_dir {
    input:
    path "summary"

    output:
    path "summary"
    
    script:
    """
    mkdir -p `readlink -f summary`
    mkdir -p summary/strelka summary/mutect2 summary/smoove summary/purple summary/svaba summary/jabba 
    mkdir -p summary/summary_snv/snv_overlap summary/summary_snv/debug
    mkdir -p summary/summary_sv/sv_overlap
    mkdir -p summary/ascat summary/facets summary/purple_purity
    mkdir -p summary/purple_driver
    """
}

process summary_snv_sv_cnv {
    debug true
    input:
    tuple val(patient_id), val(sample_id), val(wgs_seq_id), \
    file('strelka_dir'), file('mutect2_dir'), file('smoove_dir'), file('purple_dir'), file('svaba_dir'), file('jabba_dir'), \
    file('ascat_dir'), file('facets_dir')
    path "summary"

    output:
    stdout

    script:
    """
    source ${projectDir}/vcf_copy_and_filter.sh
    filter_strelka ${wgs_seq_id} strelka_dir summary/strelka
    filter_mutect2 ${wgs_seq_id} mutect2_dir summary/mutect2
    cp_smoove ${wgs_seq_id} smoove_dir summary/smoove
    cp_purple ${wgs_seq_id} purple_dir summary/purple
    cp_svaba ${wgs_seq_id} svaba_dir summary/svaba
    cp_jabba ${wgs_seq_id} jabba_dir summary/jabba
    cp -rfL ascat_dir summary/ascat/${wgs_seq_id}
    cp -rfL facets_dir summary/facets/${wgs_seq_id}
    cp purple_dir/${wgs_seq_id}.purple.purity.tsv summary/purple_purity/

    source ${projectDir}/conda_reset.sh
    python ${projectDir}/snv_overlap.py summary/mutect2/${wgs_seq_id}.snv.vcf.gz summary/strelka/${wgs_seq_id}.snv.vcf.gz \
    summary/summary_snv/snv_overlap/${wgs_seq_id}.vcf >summary/summary_snv/debug/${wgs_seq_id}.txt
    python ${projectDir}/sv_overlap.py \
    "summary/smoove/${wgs_seq_id}.sv.vcf.gz summary/purple/${wgs_seq_id}.sv.vcf.gz summary/svaba/${wgs_seq_id}.sv.vcf.gz summary/jabba/${wgs_seq_id}.sv.vcf.gz" \
    summary/summary_sv/sv_overlap/${wgs_seq_id}.tsv "smoove purple svaba jabba"
    """
}

process copy_purple_driver{
    debug true
    input:
    tuple val(patient_id), val(sample_id), val(wgs_seq_id), \
    file('strelka_dir'), file('mutect2_dir'), file('smoove_dir'), file('purple_dir'), file('svaba_dir'), file('jabba_dir'), \
    file('ascat_dir'), file('facets_dir')
    path "summary"

    output:
    stdout

    script:
    """
    cp purple_dir/${wgs_seq_id}.driver.catalog.somatic.tsv summary/purple_driver/
    """
}