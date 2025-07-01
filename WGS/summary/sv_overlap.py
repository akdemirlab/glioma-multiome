#!/usr/bin/env python
import sys
import vatect as vt

if __name__ == "__main__":
    file_paths_str = sys.argv[1]
    output_file = sys.argv[2]
    tools_str = sys.argv[3]
    file_paths = file_paths_str.split()
    tools = tools_str.split()

    df_svs = []
    for i in range(len(file_paths)):
        tool = tools[i]
        #print(file_paths[i], tool)
        df_sv = vt.read_sv_vcf(vcf_file=file_paths[i], tool=tool)
        df_filtered = vt.sv_filter(df_sv)
        df_filtered.to_csv(file_paths[i]+'.filtered.tsv', index=False, sep='\t')
        df_svs += [df_filtered]
    merge_filtered, merge = vt.sv_merge(df_svs, tools)
    merge_filtered.to_csv(output_file, index=False, sep='\t')
