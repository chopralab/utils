import os
import sys
from auto_mine_pubchem import PubChem_Miner
from tqdm import tqdm
import json
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--infile", type=str, required=True)
parser.add_argument("--outfile", type=str)
args = parser.parse_args()

df = pd.read_csv(args.infile, sep=',')
num_vend = []
is_avail = []

for elem in tqdm(df['smiles'].tolist()):
    json_data = PubChem_Miner.get_json_info(elem)
    if json_data["SourceCategories"]["Categories"][0]["Category"] != "Chemical Vendors":
        num_vend.append(0)
        is_avail.append(False)
    else:
        print("Number of Vendors: " + str(len(json_data["SourceCategories"]["Categories"][0]["Sources"])))
        num_vend.append(len(json_data["SourceCategories"]["Categories"][0]["Sources"]))
        if len(json_data["SourceCategories"]["Categories"][0]["Sources"]) == 0:
            is_avail.append(False)
        else:
            is_avail.append(True)

df["has_vendors"] = is_avail
df["number_of_vendors"] = num_vend
if args.outfile is None:
    args.outfile = args.infile.split(".")[0] + "_w_num_vend.csv"
df.to_csv(args.outfile, index=False, sep=",")
print("Saved outfile as " + args.outfile)