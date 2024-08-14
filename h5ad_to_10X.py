#!/usr/bin/python3

# load in packages needed
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import sys
import os
import argparse
import pathlib
import base64
import scipy
from scipy import io

# declaring parser to store user inputs
parser = argparse.ArgumentParser(description='Takes adata scRNA-seq object and splits it into counts, barcodes, features, to be read into Seurat as a 10X directory.')
parser.add_argument('-f', '--file', help='The file path to the adata object from Scanpy', required=True, type=pathlib.Path)
parser.add_argument('-o', '--outdir', help='Name of directory for the split data to be stored. If it does not exist it will be created', required=True, type=pathlib.Path)

# set the arguments
args = parser.parse_args()

# define paths
input_data_path = str(args.file)
output_data_path = str(args.outdir) + '/matrix_files/'

# check if the directory already exists
if not os.path.exists(output_data_path):
    # create the directory
    os.makedirs(output_data_path)
    print('Directory created successfully! Writing files...')
else:
    print('Directory already exists. Writing files...')


# reading in `adata` and set to raw
adata = sc.read_h5ad(input_data_path)
adata = adata.raw.to_adata()

# write the barcodes to a TSV to be stored in new directory
print('Writing barcodes to barcodes.tsv...')
with open(output_data_path + 'barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')
print('Finished!')
print(' ')

# writing genes to TSV
print('Writing genes to features.tsv...')
with open(output_data_path + 'features.tsv', 'w') as f:
    for item in ['\t'.join([x, x, 'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')
print('Finished!')
print(' ')

# write the actual matrix file
print('Writing counts matrix to matrix.mtx...')
io.mmwrite(output_data_path + 'matrix.mtx', adata.X.T)
print('Finished!')
print(' ')

# gzip files in `matrix files` as it is now
print('gzipping matrix, barcodes, and features...')
os.system(f"gzip {output_data_path}*")
print('Finished gzipping matrix, barcodes, and features')
print(' ')

# get metadata
print('Writing metadata to metadata.csv...')
adata.obs.to_csv(output_data_path + 'metadata.csv')
print('Finished!')
print(' ')

# get UMAP coords and write to CSV
print('Pulling UMAP coordinates, writing to umap_coords.csv...')
umap_coords = pd.DataFrame(data=adata.obsm['X_umap'])
umap_coords.to_csv(output_data_path + 'umap_coords.csv')
print('Finished!')

print(' ')
print(' ')
print('Done! :)')
print(' ')
print('Please check that all files have been written correctly and are not empty.')
