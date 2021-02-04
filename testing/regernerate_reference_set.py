#!/usr/bin/env python
import os
import argparse, os

parser = argparse.ArgumentParser(description='Generates new itesting data_filteres for testing')
parser.add_argument('--inp', help='directory for filtering out improper files 4 testing ', type=str, required=True)
parser.add_argument('--out', help='directory for filtering results. Can by directory or file with tar.gz extention ', type=str, required=True)
args = parser.parse_args()
#regenerates reference data set based on files in testingdata/calculated folder
calcs=args.inp
compress = False
if args.out.find(".tar.gz") > -1:
    compress = True
    args.out = args.out.replace(".tar.gz","")
refs=args.out
include_patterns=["*withr1primer_16S.names.txt",
                  "*on_target.txt",
                  "*cluster_genus_size.csv",
                  "*all_clusters.csv",
                  "*chosen_clusters.csv",
                  "*_counts_per_genus.tsv",
                  "*_info_on_pseudocontig_and_contigs.csv"]
include_cmd_part=""
for p in include_patterns:
    include_cmd_part+=f' --include "{p}" '
os.system(f'rsync --delete -r --progress {calcs}/ {refs}/  {include_cmd_part} --include "*/" --exclude="*" ')
if compress :
  os.system(f"tar cvzf {args.out+'.tar.gz'} {args.out} ; rm -r {args.out}  " )
print(f'rsync --delete -r --progress {calcs}/ {refs}/  {include_cmd_part} --include "*/" --exclude="*" ')
