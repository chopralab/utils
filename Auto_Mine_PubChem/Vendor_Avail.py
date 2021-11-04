import os
import sys
from PubChem_Miner import PubChem_Miner
from tqdm import tqdm
import json
import argparse
import pandas as pd
import numpy as np
from multiprocessing import Pool
from time import sleep

parser = argparse.ArgumentParser()
parser.add_argument("--infile", type=str, required=True)
parser.add_argument("--outfile", type=str)
parser.add_argument("--nprocs", type=int, default=1)
parser.add_argument("--delay", type=float, default=0.5)
parser.add_argument("--printing", action="store_true")
args = parser.parse_args()

def wrapper(smiles):
    return (smiles,args.delay,args.printing)

def parse_json(data):
        smiles = data[0]
        if data[1] is None:
            df_out.loc[len(df_out)] = [smiles,0,False]
        else:
            flag = False
            for i in range(len(data[1]["SourceCategories"]["Categories"])):
                if data[1]["SourceCategories"]["Categories"][i]["Category"] == "Chemical Vendors":
                    flag = True
                    num_vendors = len(data[1]["SourceCategories"]["Categories"][i]["Sources"])
                    if args.printing:
                        print("Number of Vendors: " + str(num_vendors))
                    if num_vendors == 0:
                        df_out.loc[len(df_out)] = [smiles,num_vendors,False]
                    else:
                        df_out.loc[len(df_out)] = [smiles,num_vendors,True]
            if flag == False:
                df_out.loc[len(df_out)] = [smiles,0,False]

df = pd.read_csv(args.infile, sep=',')
df_out = pd.DataFrame(columns=["smiles","is_com_avail","num_vendors"])
smiles_list = df['smiles'].to_list()

with Pool(args.nprocs) as pool:
    it = pool.imap_unordered(PubChem_Miner.get_json_info, smiles_list)
    it = tqdm(it, total=len(smiles_list), desc="Process JSON")
    for data in it:
        parse_json(data)

if args.outfile is None:
    args.outfile = args.infile.split(".")[0] + "_w_num_vend.csv"
df_out.to_csv(args.outfile, index=False, sep=",")
print("Saved outfile as " + args.outfile)