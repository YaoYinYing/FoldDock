#!/usr/bin/env python3
import argparse
import sys
import os
import numpy as np
import glob
from collections import defaultdict
from scipy.spatial import distance
import pdb

parser = argparse.ArgumentParser(description = '''Fetch the plDDT and calculate several metrics across each chain, both chains and the interface.''')


parser.add_argument('--model', type= str, default=sys.stdin, help = 'Path to modeled complexe.')
parser.add_argument('--metric', type= str, default=sys.stdin, help = 'Path to features for complex.')
parser.add_argument('--it', type= int, default=8, help = 'Interface threshold in Ångström (how close the atoms have to be).')
parser.add_argument('--fetch_atoms', type= str, default='CB', help = 'Atoms to fetch. E.g. CA,CB')
parser.add_argument('--cbr', type= int, default=200, help = 'Chain break residues (how many residues that were inserted as a chain break).')
parser.add_argument('--l1', type= int, default=sys.stdin, help = 'length of chain 1.')
parser.add_argument('--l2', type= int, default=sys.stdin, help = 'length of chain 2.')
parser.add_argument('--outdir', type= str, default='./', help = 'Path to output directory. Include /in end')


################FUNCTIONS#################

def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record



def parse_interface(pdbfile, l1, it, fetch_atoms):
    '''Parse the interface residues based on interactions between the individual chains
    '''

    atoms = []
    residue_numbers = []
    coords = []

    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['atm_name'] not in fetch_atoms:
                continue
            atoms.append(record['atm_name'])
            residue_numbers.append(record['res_no'])
            coords.append([record['x'],record['y'],record['z']])

    #Convert to arrays
    atoms, residue_numbers, coords = np.array(atoms), np.array(residue_numbers), np.array(coords)

    #Get res nos
    #Get chain cut
    ch1_res_nos = np.argwhere(residue_numbers<=l1)[:,0] #All residue numbers
    ch2_res_nos =  np.argwhere(residue_numbers>l1)[:,0]

    #Calculate all distances between chains and fetch interface residues
    #I can't decrease the search space here, since many residues may interact with many others
    interface_residues = []
    for i in range(len(ch1_res_nos)):
        for j in range(len(ch2_res_nos)):
            #Atom-atom distance
            dist = distance.euclidean(coords[i], coords[len(ch1_res_nos)+j]) #Need to add l1 to get the right coords
            if dist < it:
                #Save residues
                interface_residues.append(residue_numbers[i])
                interface_residues.append(residue_numbers[len(ch1_res_nos)+j])


    return np.array(interface_residues)


def calc_plddt_metrics(plDDT, interface_residues, l1):
    '''Parse the plDDT and calculate
    1. Average interface plDDT (the interface is defined as two atoms from different chains being within 6 Å from eachother)
    2. Average plDDT of each single chain
    3. Average plDDT over the entire complex
    '''
    #Interface
    if len(interface_residues)>0:
        if_plddt_av = np.average(plDDT[interface_residues])
        if_plddt_std = np.std(plDDT[interface_residues])
    else:
        if_plddt_av = 0
        if_plddt_std = 0

    #Single chain
    ch1_plddt_av = np.average(plDDT[:l1])
    ch1_plddt_std = np.std(plDDT[:l1])
    ch2_plddt_av = np.average(plDDT[l1:])
    ch2_plddt_std = np.std(plDDT[l1:])
    #Both chains
    plddt_av = np.average(plDDT)
    plddt_std = np.std(plDDT)

    return [if_plddt_av, if_plddt_std, ch1_plddt_av, ch1_plddt_std, ch2_plddt_av, ch2_plddt_std, plddt_av, plddt_std]

#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
model = args.model
metric = args.metric
#meta = pd.read_csv(args.meta[0])
cbr = args.cbr
l1 = args.l1
l2 = args.l2
it = args.it #interface threshold in Å
fetch_atoms = args.fetch_atoms.split(',')
#df_suffix = args.df_suffix[0]
#mode=args.mode[0]
#outdir = args.outdir[0]
#Save data
ids1 = []
ids2 = []
plDDT_metrics = []
interface_atom_num = [] #Number of atoms involved in the interface
interface_res_num = [] #Number of unique residues involved in the interface
metric_names = []
model_names = []
#Parse
metric_pkl = np.load(metric, allow_pickle=True)
metric_name = metric.split('/')[-1][:-4]
plDDT = metric_pkl['plddt']
model_name = model.split('/')[-1][:-4]
#Get interface residues
interface_residues = parse_interface(model, l1, it, fetch_atoms)
#Update the indices of the second chain due to the chain break
#Also decrease with 1 since the resids are 1 indexed and the plDDT 0 indexed
interface_residues[np.argwhere(interface_residues>l1)]-=cbr
interface_residues -= 1
#Calculate metrics
if len(plDDT)!=l1+l2:
    print('Missing residues in plDDT', len(plDDT), 'vs', l1+l2)
    plddt_metric = [0]*8
else:
    plddt_metric = calc_plddt_metrics(plDDT, interface_residues, l1)
#Save
#ids1.append(row.id1)
#ids2.append(row.id2)
#plDDT_metrics.append(plddt_metric)
#interface_atom_num.append(len(interface_residues))
#interface_res_num.append(np.unique(interface_residues).shape[0])
#metric_names.append(metric_name)
#model_names.append(model_name)

#Create df
print ('{},{},{},{},{},{},{},{},{},{}\n'\
       .format(round(plddt_metric[0], 3), round(plddt_metric[1], 3),
               round(plddt_metric[2], 3), round(plddt_metric[3], 3),
               round(plddt_metric[4], 3), round(plddt_metric[5], 3),
               round(plddt_metric[6], 3), round(plddt_metric[7], 3),
               len(interface_residues), 
               np.unique(interface_residues).shape[0]))

#plDDT_metrics = np.array(plDDT_metrics)
#results_df = pd.DataFrame()
#results_df['if_plddt_av']= plDDT_metrics[:,0]
#results_df['if_plddt_std']= plDDT_metrics[:,1]
#results_df['ch1_plddt_av']= plDDT_metrics[:,2]
#results_df['ch1_plddt_std']= plDDT_metrics[:,3]
#results_df['ch2_plddt_av']= plDDT_metrics[:,4]
#results_df['ch2_plddt_std']= plDDT_metrics[:,5]
#results_df['plddt_av']= plDDT_metrics[:,6]
#results_df['plddt_std']= plDDT_metrics[:,7]
#results_df['num_atoms_in_interface'] = interface_atom_num
#results_df['num_res_in_interface'] = interface_res_num
#results_df['metric_name'] = metric_names
#results_df['model_name'] = model_names
##Save
#results_df.to_csv(outdir+'plddt_metrics_'+df_suffix+'.csv')
