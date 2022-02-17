import pandas as pd


"""Function accepts a RSeQC output directory and compiles all sample information from {sample}.GC.xls
Args:
    snakemake.input (list): list of globbed wildcards from RSeQC {sample}/{sample}.GC.xls
    project_title (str): Project title for compiled %GC distributions

Returns:
    Compiled %GC distributions by sample as tab delimited file.
"""


tables = [pd.read_csv(fh, sep = "\t", names = [fh.split('/')[-2]], skiprows = 1) for fh in snakemake.input]
joined_table = pd.concat(tables, axis=1)
joined_table_sorted = joined_table.reindex(sorted(joined_table.columns), axis = 1)
joined_table_sorted.to_csv(snakemake.output[0], sep='\t')