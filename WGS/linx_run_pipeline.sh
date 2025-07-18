#!/bin/bash

#the script was copied from https://github.com/hartwigmedical/hmftools/blob/master/pipeline/scripts/run_pipeline
# example command:
# ./pipeline/scripts/run_pipeline ./pipeline/scripts /sample_data_dir/ /ref_data_dir/ /tools_dir/ "COLO829T,COLO829R" V37 WGS/PANEL 10 16

# required arguments
scripts_dir=$1 && shift
samples_dir=$1 && shift
resources_dir=$1 && shift
tools_dir=$1 && shift
sample_data=$1 && shift
ref_genome_version=$1 && shift
run_mode=$1 && shift
threads=$1 && shift
max_memory=$1 && shift
output_dir=$1 && shift
hg38_ref=/rsrch4/home/neurosurgery/akdemirlab/references/GRCh38/GRCh38.primary_assembly.genome.fa

# input argument checks
if [[ ! -d "${scripts_dir}" ]]; then
  echo "Missing scripts directory: ${scripts_dir}"
  exit 1
fi

if [[ ! -d "${samples_dir}" ]]; then
  echo "Missing samples directory: ${samples_dir}"
  exit 1
fi

if [[ ! -d "${resources_dir}" ]]; then
  echo "Missing resources directory: ${resources_dir}"
  exit 1
fi

if [[ ! -d "${tools_dir}" ]]; then
  echo "Missing tools directory: ${tools_dir}"
  exit 1
fi

if [ ! "${run_mode}" == "WGS" ] && [ ! "${run_mode}" == "PANEL" ]; then
  echo "Invalid run mode: ${run_mode}, must be 'WGS' or 'PANEL'"
  exit 1
fi


sample_array=($(echo $sample_data | tr "," " "))
tumor_id=${sample_array[0]}
sample_count=${#sample_array[@]}

#sample_dir=${samples_dir}/${tumor_id}
mkdir -p ${output_dir}/$tumor_id
sample_dir=${output_dir}/$tumor_id

# assume BAM names match tumor and reference names
tumor_bam=${samples_dir}/${tumor_id}*.bam

# set up reference if provided
reference_id="none"
reference_bam="none"

if [ "${sample_count}" == 2 ] && [ "${run_mode}" == "WGS" ]; then
  reference_id=${sample_array[1]}
  reference_bam=${samples_dir}/${reference_id}*.bam
fi


if [ "${max_memory}" == "" ]; then
  max_memory=16
fi

echo "Running HMF pipeline in ${run_mode} mode"

echo "Environment:"
echo "  scripts dir: ${scripts_dir}"
echo "  samples dir: ${samples_dir}"
echo "  resources dir: ${resources_dir}"
echo "  tools dir: ${tools_dir}"
echo "  sample count: ${sample_count}"
echo "  tumorId: ${tumor_id}"
echo "  referenceId: ${reference_id}"
echo "  ref genome version: ${ref_genome_version}"
echo "  run mode: ${run_mode}"
echo "  threads: ${threads}"
echo "  memory: ${max_memory}"

echo "Sample output dir: ${sample_dir}"

# set resource files
if [ "${ref_genome_version}" == "V37" ]; then
  echo "Reference genome version GRCh37"

  # Reference genome
  ref_genome=${resources_dir}/ref_genome/Homo_sapiens.GRCh37.GATK.illumina.fasta

  # Common
  ensembl_dir=${resources_dir}/common/ensembl_data/
  driver_gene_panel=${resources_dir}/common/DriverGenePanel.37.tsv

  # Point mutations (for Sage, Pave)
  sage_somatic_hotspots=${resources_dir}/variants/KnownHotspots.somatic.37.vcf.gz
  sage_germline_hotspots=${resources_dir}/variants/KnownHotspots.germline.37.vcf.gz
  sage_panel_bed=${resources_dir}/variants/ActionableCodingPanel.37.bed.gz
  sage_coverage_bed=${resources_dir}/variants/CoverageCodingPanel.37.bed.gz
  sage_germline_panel_bed=${resources_dir}/variants/ActionableCodingPanel.germline.37.bed.gz
  high_confidence_bed=${resources_dir}/variants/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz
  somatic_pon_file=${resources_dir}/variants/SageGermlinePon.1000x.37.tsv.gz
  mappability_bed=${resources_dir}/variants/mappability_150.37.bed.gz
  clinvar_vcf=${resources_dir}/variants/clinvar.37.vcf.gz
  germline_blacklist_bed=${resources_dir}/variants/KnownBlacklist.germline.37.bed
  germline_blacklist_vcf=${resources_dir}/variants/KnownBlacklist.germline.37.vcf.gz
  gnomad_path=${resources_dir}/variants/gnomad_variants_v37.csv.gz
  pon_artefact_file=${resources_dir}/variants/PanelArtefacts.37.tsv.gz

  # SVs (for Gripss, Linx)
  sv_hotspot_file=${resources_dir}/sv/known_fusions.37.bedpe
  sv_blacklist_bed=${resources_dir}/sv/sv_prep_blacklist.37.bed
  gridss_blacklist_bed=${resources_dir}/sv/gridss_blacklist.37.bed.gz
  gridss_config=${resources_dir}/sv/gridss.properties
  sv_pon_file=${resources_dir}/sv/sv_pon.37.bedpe.gz
  sgl_pon_file=${resources_dir}/sv/sgl_pon.37.bed.gz
  repeat_mask_file=${resources_dir}/sv/repeat_mask_data.37.fa.gz
  fragile_site_file=${resources_dir}/sv/fragile_sites_hmf.37.csv
  line_element_file=${resources_dir}/sv/line_elements.37.csv
  known_fusion_file=${resources_dir}/sv/known_fusion_data.37.csv
  
  # Copy number (for Amber, Cobalt, Purple)
  amber_loci_vcf=${resources_dir}/copy_number/GermlineHetPon.37.vcf.gz
  gc_profile=${resources_dir}/copy_number/GC_profile.1000bp.37.cnp
  tumor_only_diploid_bed=${resources_dir}/copy_number/DiploidRegions.37.bed.gz
  germline_del_freq_file=${resources_dir}/copy_number/cohort_germline_del_freq.37.csv

  target_region_normalisation=${resources_dir}/copy_number/target_regions_normalisation.37.tsv
  target_regions_definition=${resources_dir}/variants/CoverageCodingPanel.37.bed.gz
  target_regions_ratios=${resources_dir}/copy_number/target_regions_ratios.37.tsv
  target_regions_msi_indels=${resources_dir}/copy_number/target_regions_msi_indels.37.tsv

  # Immune (for Lilac)
  lilac_resource_dir=${resources_dir}/immune/

else
  echo "Reference genome version GATK GRCh38"
  ref_genome=$hg38_ref
  ensembl_dir=${resources_dir}/common/ensembl_data/
  driver_gene_panel=${resources_dir}/common/DriverGenePanel.38.tsv

  # Point mutations (for Sage, Pave)
  sage_somatic_hotspots=${resources_dir}/variants/KnownHotspots.somatic.38.vcf.gz
  sage_germline_hotspots=${resources_dir}/variants/KnownHotspots.germline.38.vcf.gz
  sage_panel_bed=${resources_dir}/variants/ActionableCodingPanel.38.bed.gz
  sage_coverage_bed=${resources_dir}/variants/CoverageCodingPanel.38.bed.gz
  high_confidence_bed=${resources_dir}/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz
  somatic_pon_file=${resources_dir}/variants/SageGermlinePon.98x.38.tsv.gz
  mappability_bed=${resources_dir}/variants/mappability_150.38.bed.gz
  clinvar_vcf=${resources_dir}/variants/clinvar.38.vcf.gz
  germline_blacklist_bed=${resources_dir}/variants/KnownBlacklist.germline.38.bed
  germline_blacklist_vcf=${resources_dir}/variants/KnownBlacklist.germline.38.vcf.gz
  gnomad_path=${resources_dir}/variants/gnomad/
  pon_artefact_file=${resources_dir}/variants/PanelArtefacts.38.tsv.gz

  # SVs (for Gripss, Linx)
  sv_hotspot_file=${resources_dir}/sv/known_fusions.38.bedpe
  sv_blacklist_bed=${resources_dir}/sv/sv_prep_blacklist.38.bed
  gridss_blacklist_bed=${resources_dir}/sv/gridss_blacklist.38.bed.gz
  gridss_config=${resources_dir}/sv/gridss.properties
  sv_pon_file=${resources_dir}/sv/sv_pon.38.bedpe.gz
  sgl_pon_file=${resources_dir}/sv/sgl_pon.38.bed.gz
  repeat_mask_file=${resources_dir}/sv/repeat_mask_data.38.fa.gz
  fragile_site_file=${resources_dir}/sv/fragile_sites_hmf.38.csv
  line_element_file=${resources_dir}/sv/line_elements.38.csv
  known_fusion_file=${resources_dir}/sv/known_fusion_data.38.csv

  # Copy number (for Amber, Cobalt, Purple)
  amber_loci_vcf=${resources_dir}/copy_number/GermlineHetPon.38.vcf.gz
  gc_profile=${resources_dir}/copy_number/GC_profile.1000bp.38.cnp
  tumor_only_diploid_bed=${resources_dir}/copy_number/DiploidRegions.38.bed.gz
  germline_del_freq_file=${resources_dir}/copy_number/cohort_germline_del_freq.38.csv

  target_region_normalisation=${resources_dir}/copy_number/target_regions_normalisation.38.tsv
  target_regions_definition=${resources_dir}/copy_number/target_regions_definition.38.bed
  target_regions_ratios=${resources_dir}/copy_number/target_regions_ratios.38.tsv
  target_regions_msi_indels=${resources_dir}/copy_number/target_regions_msi_indels.38.tsv

  # Immune (for Lilac)
  lilac_resource_dir=${resources_dir}/immune/
fi

# set tool links
sage_jar=${tools_dir}/sage.jar
pave_jar=${tools_dir}/pave.jar
amber_jar=${tools_dir}/amber.jar
gripss_jar=${tools_dir}/gripss.jar
cobalt_jar=${tools_dir}/cobalt.jar
purple_jar=${tools_dir}/purple.jar
linx_jar=${tools_dir}/linx.jar
lilac_jar=${tools_dir}/lilac.jar

#circos=${tools_dir}/circos
circos=$(which circos)

# Amber
amber_dir=${sample_dir}/amber

${scripts_dir}/run_amber ${amber_jar} \
  ${tumor_id} ${tumor_bam} ${reference_id} ${reference_bam} \
  ${amber_dir} \
  ${run_mode} ${ref_genome_version} ${ref_genome} \
  ${amber_loci_vcf} ${threads} ${max_memory} \

# Cobalt
cobalt_dir=${sample_dir}/cobalt

${scripts_dir}/run_cobalt ${cobalt_jar} \
  ${tumor_id} ${tumor_bam} ${reference_id} ${reference_bam} \
  ${cobalt_dir} \
  ${run_mode} ${ref_genome} \
  ${gc_profile} ${tumor_only_diploid_bed} ${target_region_normalisation} \
  ${threads} ${max_memory} \


# Sage somatic

if [ "${reference_id}" != "none" ]; then
  sage_somatic_dir=${sample_dir}/sage_somatic
  sage_vcf=${sage_somatic_dir}/${tumor_id}.sage.somatic.vcf.gz
else
  sage_somatic_dir=${sample_dir}/sage
  sage_vcf=${sage_somatic_dir}/${tumor_id}.sage.vcf.gz
fi

${scripts_dir}/run_sage_somatic ${sage_jar} \
  ${tumor_id} ${tumor_bam} ${reference_id} ${reference_bam} \
  ${sage_somatic_dir} ${sage_vcf} \
  ${run_mode} ${ref_genome_version} ${ref_genome} \
  ${ensembl_dir} ${sage_somatic_hotspots} ${sage_panel_bed} ${sage_coverage_bed} ${high_confidence_bed} \
  ${threads} ${max_memory} \


# Pave somatic
if [ "${reference_id}" != "none" ]; then
  pave_somatic_dir=${sample_dir}/pave_somatic
  pave_vcf=${pave_somatic_dir}/${tumor_id}.sage.somatic.pave.vcf.gz
else
  pave_somatic_dir=${sample_dir}/pave
  pave_vcf=${pave_somatic_dir}/${tumor_id}.sage.pave.vcf.gz
fi

${scripts_dir}/run_pave_somatic ${pave_jar} \
  ${tumor_id}  \
  ${sage_vcf} ${pave_somatic_dir} ${pave_vcf} \
  ${run_mode} ${ref_genome_version} ${ref_genome} \
  ${ensembl_dir} ${driver_gene_panel} ${somatic_pon_file} ${pon_artefact_file} ${mappability_bed} ${gnomad_path} \


# Sage germline
pave_germline_vcf="none"
if [ "${reference_id}" != "none" ]; then
  
  sage_germline_dir=${sample_dir}/sage_germline

  ${scripts_dir}/run_sage_germline ${sage_jar} \
    ${tumor_id} ${tumor_bam} ${reference_id} ${reference_bam} \
    ${sage_germline_dir} \
    ${ref_genome_version} ${ref_genome} \
    ${ensembl_dir} ${sage_germline_hotspots} ${sage_panel_bed} ${high_confidence_bed} \
    ${threads} ${max_memory} \

  # Pave germline
  pave_germline_dir=${sample_dir}/pave_germline
  pave_germline_vcf=${pave_germline_dir}/${tumor_id}.sage.germline.pave.vcf.gz

  ${scripts_dir}/run_pave_germline ${pave_jar} \
    ${tumor_id}  \
    ${sage_vcf} ${pave_germline_dir} ${pave_germline_vcf} \
    ${ref_genome_version} ${ref_genome} \
    ${ensembl_dir} ${driver_gene_panel} ${mappability_bed} ${clinvar_vcf} ${germline_blacklist_bed} ${germline_blacklist_vcf} \

fi


# Gridss 
gridss_dir=${sample_dir}/gridss
gridss_vcf=${gridss_dir}/${tumor_id}.gridss.unfiltered.vcf.gz

${scripts_dir}/run_sv_calling ${tools_dir} \
  ${tumor_id} ${tumor_bam} ${reference_id} ${reference_bam} \
  ${gridss_dir} ${gridss_vcf} \
  ${ref_genome_version} ${ref_genome} \
  ${sv_blacklist_bed} ${sv_hotspot_file} \
  ${gridss_blacklist_bed} ${gridss_config} \
  ${threads} ${max_memory} \


# Gripss somatic
if [ "${reference_id}" != "none" ]; then
  gripss_somatic_dir=${sample_dir}/gripss_somatic
else
  gripss_somatic_dir=${sample_dir}/gripss
fi

${scripts_dir}/run_gripss_somatic ${gripss_jar} \
  ${tumor_id}  ${reference_id} \
  ${gripss_somatic_dir} ${gridss_vcf} \
  ${run_mode} ${ref_genome_version} ${ref_genome} \
  ${sv_hotspot_file} ${sv_pon_file} ${sgl_pon_file} ${repeat_mask_file} ${target_regions_definition} \


if [ "${reference_id}" != "none" ]; then

  # Gripss germline
  gripss_germline_dir=${sample_dir}/gripss_germline

  ${scripts_dir}/run_gripss_germline ${gripss_jar} \
    ${reference_id} \
    ${gripss_germline_dir} ${gridss_vcf} \
    ${ref_genome_version} ${ref_genome} \
    ${sv_hotspot_file} ${sv_pon_file} ${sgl_pon_file}
fi


# Purple
purple_dir=${sample_dir}/purple

sv_vcf=${gripss_somatic_dir}/${tumor_id}.gripss.filtered.somatic.vcf.gz
sv_unfiltered_vcf=${gripss_somatic_dir}/${tumor_id}.gripss.somatic.vcf.gz

${scripts_dir}/run_purple ${purple_jar} \
  ${tumor_id} ${reference_id} \
  ${sv_vcf} ${sv_unfiltered_vcf} ${pave_vcf} ${pave_germline_vcf} ${amber_dir} ${cobalt_dir} ${purple_dir} \
  ${run_mode} ${ref_genome_version} ${ref_genome} \
  ${gc_profile} ${sage_somatic_hotspots} ${sage_germline_hotspots} \
  ${driver_gene_panel} ${ensembl_dir} ${germline_del_freq_file} \
  ${target_regions_definition} ${target_regions_ratios} ${target_regions_msi_indels} \
  ${threads} ${circos} ${max_memory} \


# Linx somatic
if [ "${reference_id}" != "none" ]; then
  linx_somatic_dir=${sample_dir}/linx_somatic
else
  linx_somatic_dir=${sample_dir}/linx
fi

sv_vcf=${purple_dir}/${tumor_id}.purple.sv.vcf.gz

${scripts_dir}/run_linx_somatic ${linx_jar} \
  ${tumor_id} \
  ${sv_vcf} ${purple_dir} ${linx_somatic_dir} \
  ${ref_genome_version} ${fragile_site_file} \
  ${line_element_file} ${ensembl_dir} ${driver_gene_panel} ${known_fusion_file} \
  ${circos} \


if [ "${reference_id}" != "none" ]; then

  # Linx germline
  linx_germline_dir=${sample_dir}/linx_germline

  sv_germline_vcf=${gripss_germline_dir}/${reference_id}.gripss.filtered.germline.vcf.gz

  ${scripts_dir}/run_linx_germline ${linx_jar} \
    ${tumor_id} \
    ${sv_germline_vcf} ${linx_germline_dir} \
    ${ref_genome_version} ${line_element_file} ${ensembl_dir} ${driver_gene_panel}
fi

# Lilac
lilac_dir=${sample_dir}/lilac

somatic_vcf=${purple_dir}/${tumor_id}.purple.somatic.vcf.gz
gene_copy_number=${purple_dir}/${tumor_id}.purple.cnv.gene.tsv

${scripts_dir}/run_lilac ${lilac_jar} \
  ${tumor_id} ${tumor_bam} ${reference_id} ${reference_bam} \
  ${gene_copy_number} ${somatic_vcf} ${lilac_dir} \
  ${ref_genome_version} ${ref_genome} ${lilac_resource_dir} \
  ${threads} \


echo "Pipeline complete"
