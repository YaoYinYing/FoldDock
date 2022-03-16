import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Match OX identifiers for interacting chains.''')

parser.add_argument('--oxdir', nargs=1, type= str, default=sys.stdin, help = 'Path to data.')
parser.add_argument('--pdbmeta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')


def match_ox(ox1,ox2):
    '''Match the organism identifiers "OX" from the interacting pairs and calculate
     - number of possible combinations
     - median and std in hits per species
    '''
    #Remove the zeros (no OX)
    ox1 = ox1[np.argwhere(ox1!=0)]
    ox2 = ox2[np.argwhere(ox2!=0)]
    matching_ox = np.intersect1d(ox1,ox2)
    nunique_ox1 = np.unique(ox1).shape[0]
    nunique_ox2 = np.unique(ox2).shape[0]
    ncombos = []
    for ox in matching_ox:
        match1 = np.argwhere(ox1==ox).shape[0]
        match2 = np.argwhere(ox2==ox).shape[0]
        ncombos.append(match1*match2)

    return nunique_ox1, nunique_ox2, matching_ox.shape[0], np.sum(ncombos), np.median(ncombos), np.std(ncombos), np.array(ncombos)

#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
oxdir = args.oxdir[0]
pdbmeta = pd.read_csv(args.pdbmeta[0])
outdir = args.outdir[0]

oxfiles = glob.glob(oxdir+'*_species.npy')

#Go through all interactions
#Save
results = {'id1':[], 'id2':[],'n_unique1' : [], 'n_unique2' : [], 'n_matching_ox' : [],
            'ncombos_total' : [], 'median_combos' : [], 'std_combos' : []}
all_ind_combos = []

for i in range(len(pdbmeta)):
    row = pdbmeta.loc[i]
    id1 = row.PDB+'_'+row['Chain 1']
    id2 = row.PDB+'_'+row['Chain 2']
    try:
        ox1 = np.load(oxdir+id1+'_species.npy')
        ox2 = np.load(oxdir+id2+'_species.npy')

    except:
        continue

    #Match
    nunique_ox1, nunique_ox2, matching_ox, ncombos_total, median_combos, std_combos, all_n = match_ox(ox1, ox2)
    all_ind_combos.extend([*all_n]) #For the total median and std

    results['id1'].append(id1)
    results['id2'].append(id2)
    results['n_unique1'].append(nunique_ox1)
    results['n_unique2'].append(nunique_ox2)
    results['n_matching_ox'].append(matching_ox)
    results['ncombos_total'].append(ncombos_total)
    results['median_combos'].append(median_combos)
    results['std_combos'].append(std_combos)

#Save
results_df = pd.DataFrame.from_dict(results)
results_df.to_csv(outdir+'oxstats.csv', index=False)
#Get the totals
total_median = np.median(all_ind_combos)
total_std = np.std(all_ind_combos)
print('Total median', total_median)
print('Total std', total_std)
