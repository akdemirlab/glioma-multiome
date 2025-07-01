import sys
import pandas as pd
import glob

def mutation_type(purple_files, names, likelyhood=0):
    matrix = pd.DataFrame()

    for index, file in enumerate(purple_files):
        path = glob.glob(file)[0]
        df = pd.read_csv(path, sep='\t')
        #print(df)
        for df_index, row in df.iterrows():
            #print(row['gene'])
            data = {
                "sample": names[index],
                "category": row['gene']
            }
            if row['driverLikelihood'] < likelyhood:
                continue
            if row['missense']:
                data['value'] = 'missense'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['nonsense']:
                data['value'] = 'nonsense'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['frameshift']:
                data['value'] = 'frameshift'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['splice']:
                data['value'] = 'splice'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['inframe']:
                data['value'] = 'inframe'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['frameshift']:
                data['value'] = 'frameshift'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['driver']=='DEL':
                data['value'] = 'deletion'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
            if row['driver']=='AMP':
                data['value'] = 'amplification'
                data = pd.DataFrame(data,  index=[0])
                matrix = pd.concat([matrix, data], ignore_index=True)
    return matrix

if __name__ == "__main__":
    linx_dir = '/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/summary_v2/purple_driver/'
    result_dir = "/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/results/figure_comut_v2/"
    data = pd.read_csv('/rsrch4/scratch/neurosurgery/akdemirlab/catalyst_main/meta/meta_catalyst-main_v1.0_12_06.csv')
    samples = data['wgs_seq_id']
    path = ''
    files = linx_dir + samples + '.driver.catalog.somatic.tsv'
    files = list(files)

    matrix = mutation_type(files, list(samples), likelyhood=0.75)
    matrix.to_csv(result_dir+'comut_driver_genes_l075.csv', index=False)
