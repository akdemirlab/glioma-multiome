import sys
import pandas as pd
import glob

def read_vep(path, gene_list=[], max_distance=500):
    column_names = ["Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position",
                    "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra"]
    df = pd.read_csv(path, sep='\t', comment='#', header=None, names=column_names)
    df['Symbol'] = df['Extra'].str.extract(r'SYMBOL=([^;]*);', expand=False)

    df['Distance'] = df['Extra'].str.extract(r'DISTANCE=([^;]*);', expand=False) 
    df['Distance'] = pd.to_numeric(df['Distance'], errors='coerce')

    if len(gene_list)>0:
        df = df[df['Symbol'].isin(gene_list)]
    
    df['value'] = ''
    df['category'] = df['Symbol']
    df.loc[df['Consequence'].str.contains('missense_variant'), 'value'] = 'missense'
    df.loc[df['Consequence'].str.contains('upstream_gene_variant'), 'value'] = 'upstream'
    df.loc[df['Consequence'].str.contains('stop_lost'), 'value'] = 'nonsense'
    df.loc[df['Consequence'].str.contains('start_lost'), 'value'] = 'nonsense'
    df.loc[df['Consequence'].str.contains('stop_gained'), 'value'] = 'nonsense'
    df.loc[df['Consequence'].str.contains('frameshift_variant'), 'value'] = 'frameshift'
    df.loc[df['Consequence'].str.contains('inframe_insertion'), 'value'] = 'inframe_insertion'
    df.loc[df['Consequence'].str.contains('inframe_deletion'), 'value'] = 'inframe_deletion'
    df.loc[df['Consequence'].str.contains('splice_acceptor_variant'), 'value'] = 'splice'
    df.loc[df['Consequence'].str.contains('splice_donor_variant'), 'value'] = 'splice'

    df = df[df['value']!=""]
    filtered_df = df[(df['Distance'].isna()) | (df['Distance'] <= max_distance)]
    summary = filtered_df
    summary = summary[['category', 'value']].drop_duplicates()
    return summary, filtered_df


def mutation_type(purple_files, names):
    matrix = pd.DataFrame()
    #gene_list = ['TP53', 'IDH1', 'PTEN', 'ATRX', 'EGFR', 
    #            'NF1', 'TERT', 'PIK3CA', 'CIC', 'RB1', 'PIK3R1', 'FUBP1']

    for index, file in enumerate(purple_files):
        #print(file)
        path = glob.glob(file)[0]
        summary, df = read_vep(path)
        #summary, df = read_vep(path, gene_list)
        for df_index, row in summary.iterrows():
            #print(row['gene'])
            data = {
                "sample": names[index],
                "category": row['category'],
                "value": row['value']
            }
            data = pd.DataFrame(data,  index=[0])
            matrix = pd.concat([matrix, data], ignore_index=True)
    return matrix

def all_samples():
    vep_dir = '/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2/annotation/vep_mutect2/'
    data = pd.read_csv('/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/meta/meta_catalyst-main_v1.0_12_06.csv')
    samples = data['wgs_seq_id']
    files_snv = vep_dir + samples + '.snv.tsv'
    files_snv = list(files_snv)
    matrix_snv = mutation_type(files_snv, list(samples))

    #files_indel = vep_dir + samples + '.indel.tsv'
    #files_indel = list(files_indel)
    #matrix_indel = mutation_type(files_indel, list(samples))
    #matrix = pd.concat([matrix_snv, matrix_indel], ignore_index=True)
    #matrix = matrix.drop_duplicates()
    
    matrix_snv.to_csv('/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/figure_comut_v2/comut_driver_genes_mutect2_snv.csv', index=False)
    #matrix_indel.to_csv('/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/figure_comut_v2/comut_driver_genes_mutect2_indel.csv', index=False)
    #matrix.to_csv('/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/figure_comut_v2/comut_driver_genes_mutect2.csv', index=False)
    
def test_one_sample():
    vep_file = "235842-WG01.tsv"
    gene_list = ['TP53', 'IDH1', 'PTEN', 'ATRX', 'EGFR', 
                'NF1', 'TERT', 'PIK3CA', 'CIC', 'RB1', 'PIK3R1', 'FUBP1']
    summary, df = read_vep(vep_file, gene_list)

if __name__ == "__main__":
    all_samples()
    #test_one_sample()