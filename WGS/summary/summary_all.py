'''
use metadata as index, summarise SNV, SV, CNV, ploidy and purity of all samples 

###################################################
locations of result files:

All results are in sub-directories of /rsrch4/scratch/neurosurgery/akdemirlab 

metadta: 
catalyst_wgs/data/meta_catalyst_v5.1_lye1.csv

SNV: 
catalyst_wgs/results/summary/mutect2
catalyst_wgs/results/summary/strelka
catalyst_wgs/results/summary/summary_snv/snv_overlap

SV:
catalyst_wgs/results/summary/smoove
catalyst_wgs/results/summary/purple
catalyst_wgs/results/summary/svaba
catalyst_wgs/results/summary/jabba
catalyst_wgs/results/summary/summary_sv/sv_overlap

ploidy & purity:
catalyst_wgs/results/ascat/*
catalyst_wgs/results/facets/facets/*

mtSNV:
catalyst_mt/results/somatic_mtsnv/mtsnv_atac_bowtie
catalyst_mt/results/somatic_mtsnv/mtsnv_atac_bwa
catalyst_mt/results/somatic_mtsnv/mtsnv_rna
catalyst_mt/results/somatic_mtsnv/mtsnv_wgs

###################################################
locations of scripts:

SNV & SV:
analysis/vcf_summary.sh

mtSNV:
catalyst_mt/analysis/somatic_mtSNV.py

'''

import pandas as pd
from cyvcf2 import VCF, Writer
import re
import os
import gzip
import shutil



def vcf_site_counts(vcf_file):
    vcf_obj = VCF(vcf_file)
    count = 0
    for variant in vcf_obj:
        count+=1
    return count

def vcf_site_counts2(vcf_file, passed = True):
    '''count vcf sites with PASS
    '''
    matches = []
    if vcf_file.endswith('.gz'):
        with gzip.open(vcf_file, 'rb') as f:
            for line in f:
                line = line.decode()
                if not line.startswith("#"):
                    if not passed or re.search("PASS", line):
                        matches.append(line)
    else:
        with open(vcf_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if not line.startswith("#"):
                    if not passed or re.search("PASS", line):
                        matches.append(line)
    return len(matches)

def bedpe_site_counts(bedpe_file):
    table = pd.read_csv(bedpe_file, delimiter='\t')
    return len(table)

def bedpe_site_counts_sv(bedpe_file):
    table = pd.read_csv(bedpe_file, delimiter='\t')
    grouped_table = table.groupby('svtype')
    count_df = grouped_table.size().to_frame('Count').reset_index(drop=False)
    return count_df

def ploidy_purity_ascat(stats_file):
    table = pd.read_table(stats_file, delimiter=' ', names=['name', 'value'])
    ploidy = float(table.loc[table['name']=='Ploidy', 'value'])
    purity = float(table.loc[table['name']=='rho', 'value'])
    return [ploidy, purity]

def ploidy_purity_facets(stats_file):
    table = pd.read_table(stats_file, delimiter='\t')
    ploidy = float(table.loc[0, 'ploidy'])
    purity = float(table.loc[0, 'purity'])
    return [ploidy, purity]

def ploidy_purity_purple(stats_file):
    table = pd.read_table(stats_file, delimiter='\t')
    ploidy = float(table.loc[0, 'ploidy'])
    purity = float(table.loc[0, 'purity'])
    return [ploidy, purity]

def find_files(directory, ends):
    files = []
    names = []
    for filename in os.listdir(directory):
        if filename.endswith(ends):
            sample_name = re.sub(ends, "", filename)
            files.append(os.path.join(directory, filename))
            names.append(sample_name)
    return files, names

def find_files2(directory, ends):
    '''find files in subdirectories too
    '''
    files = []
    names = []
    for root, dirs, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(ends):
                sample_name = re.sub(ends, "", filename)
                files.append(os.path.join(root, filename))
                names.append(sample_name)
    return files, names

def vcf_summary(files, names):
    df = pd.DataFrame()
    for i in range(len(names)):
        data = {
            "name": names[i],
            "counts": vcf_site_counts(files[i])
        }
        data = pd.DataFrame(data,  index=[0])
        df = pd.concat([df, data], ignore_index=True)
    return df

def vcf_passed_summary(files, names):
    df = pd.DataFrame()
    for i in range(len(names)):
        data = {
            "name": names[i],
            "counts": vcf_site_counts2(files[i], passed=True)
        }
        data = pd.DataFrame(data,  index=[0])
        df = pd.concat([df, data], ignore_index=True)
    return df

def bedpe_sv_summary(files, names):
    df = pd.DataFrame()
    for i in range(len(names)):
        data = bedpe_site_counts_sv(files[i])
        data['name'] = names[i]
        df = pd.concat([df, data], ignore_index=True)
    return df

def ascat_ploidy_summary(files, names):
    df = pd.DataFrame()
    for i in range(len(names)):
        ploidy, purity = ploidy_purity_ascat(files[i])
        data = {
            "name": names[i],
            "ploidy": ploidy,
            "purity": purity
        }
        data = pd.DataFrame(data,  index=[0])
        df = pd.concat([df, data], ignore_index=True)
    return df

def facets_ploidy_summary(files, names):
    df = pd.DataFrame()
    for i in range(len(names)):
        ploidy, purity = ploidy_purity_facets(files[i])
        data = {
            "name": names[i],
            "ploidy": ploidy,
            "purity": purity
        }
        data = pd.DataFrame(data,  index=[0])
        df = pd.concat([df, data], ignore_index=True)
    return df

def purple_ploidy_summary(files, names):
    df = pd.DataFrame()
    for i in range(len(names)):
        ploidy, purity = ploidy_purity_purple(files[i])
        data = {
            "name": names[i],
            "ploidy": ploidy,
            "purity": purity
        }
        data = pd.DataFrame(data,  index=[0])
        df = pd.concat([df, data], ignore_index=True)
    return df

if __name__ == "__main__":
    summary_path = '/Volumes/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2'
    summary_all_path='/Volumes/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_all_v2'

    def summary_snv():
        mutect_dir = os.path.join(summary_path, 'mutect2') 
        strelka_dir =  os.path.join(summary_path, 'strelka')
        purple_snv_dir = os.path.join(summary_path, 'purple_snv')
        snv_overlap_dir =  os.path.join(summary_path, 'summary_snv/snv_overlap')

        files, names = find_files(mutect_dir, '.snv.vcf.gz')
        df = vcf_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'snv_mutect_counts.csv'), index=False)

        files, names = find_files(strelka_dir, '.snv.vcf.gz')
        df = vcf_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'snv_strelka_counts.csv'), index=False)

        files, names = find_files(snv_overlap_dir, '.vcf')
        df = vcf_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'snv_overlap_counts.csv'), index=False)

    def summary_sv():
        smoove_dir = os.path.join(summary_path, 'smoove') 
        purple_dir = os.path.join(summary_path, 'purple') 
        svaba_dir = os.path.join(summary_path, 'svaba') 
        jabba_dir = os.path.join(summary_path, 'jabba') 
        sv_overlap_dir = os.path.join(summary_path, 'summary_sv/sv_overlap') 

        files, names = find_files(smoove_dir, '.sv.vcf.gz.filtered.tsv')
        df = bedpe_sv_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'sv_smoove_counts.csv'), index=False)

        files, names = find_files(purple_dir, '.sv.vcf.gz.filtered.tsv')
        df = bedpe_sv_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'sv_purple_counts.csv'), index=False)

        files, names = find_files(svaba_dir, '.sv.vcf.gz.filtered.tsv')
        df = bedpe_sv_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'sv_svaba_counts.csv'), index=False)

        files, names = find_files(jabba_dir, '.sv.vcf.gz.filtered.tsv')
        df = bedpe_sv_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'sv_jabba_counts.csv'), index=False)

        files, names = find_files(sv_overlap_dir, '.tsv')
        df = bedpe_sv_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'sv_overlap_counts.csv'), index=False)

    def summary_ploidy():
        ascat_dir = os.path.join(summary_path, 'ascat/') 
        facets_dir = os.path.join(summary_path, 'facets/') 
        purple_dir = os.path.join(summary_path, 'purple_purity/') 

        files, names = find_files2(ascat_dir, '.samplestatistics.txt')
        df = ascat_ploidy_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'ploidy_ascat.csv'), index=False)

        files, names = find_files2(facets_dir, '.stats.txt')
        df = facets_ploidy_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'ploidy_facets.csv'), index=False)

        files, names = find_files2(purple_dir, '.purple.purity.tsv')
        df = purple_ploidy_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'ploidy_purple.csv'), index=False)

    def summary_mtsnv_v2():
        mt_prj_path='/Volumes/scratch/neurosurgery/akdemirlab/catalyst_mt/'
        mtsnv_atac = os.path.join(mt_prj_path, 'results/somatic_mtsnv_v2/mtsnv_atac') 
        mtsnv_wgs = os.path.join(mt_prj_path, 'results/somatic_mtsnv_v2/mtsnv_wgs_merge')

        files, names = find_files(mtsnv_atac, '.vcf')
        df = vcf_passed_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'mtsnv_atac_counts.csv'), index=False)

        files, names = find_files(mtsnv_wgs, '.vcf')
        df = vcf_passed_summary(files, names)
        df.to_csv(os.path.join(summary_all_path, 'mtsnv_wgs_counts.csv'), index=False)

        ############
        #mt read counts)
        mt_depth_wgs = os.path.join(mt_prj_path, 'results/mt_depth_wgs_v2.csv')
        mt_depth_atac = os.path.join(mt_prj_path, 'results/mt_depth_atac_v2.csv')
        #shutil.copyfile(mt_depth_wgs, summary_all_path+'/mt_depth_wgs.csv')
        #shutil.copyfile(mt_depth_atac, summary_all_path+'/mt_depth_atac.csv')

    #summary_snv()
    #summary_sv()
    summary_ploidy()
    #summary_mtsnv_v2()