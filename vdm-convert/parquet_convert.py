import numpy as np
import pandas as pd
import tempfile
from pymol import cmd
import prody
import os
from os.path import join
import sys
import argparse

########################################
### HELPER FUNCTIONS
########################################
def create_atomgroup(df):
	atoms = prody.AtomGroup('VDM')

	coords = np.array(df[['c_x', 'c_y', 'c_z']])
	res_names = np.array(df['resname'])
	res_nums = np.array(df['resnum'])
	atom_names = np.array(df['name'])
	chain_ids = np.array(df['chain'])

	atoms.setCoords(coords)
	atoms.setNames(atom_names)
	atoms.setResnums(res_nums)
	atoms.setResnames(res_names)
	atoms.setChids(chain_ids)

	return atoms

def generate_filename(dir,name,ext):
	#Check if a file with this name exists. If so, append a number to the end
	name_satisfied = False
	full_name = join(dir, name)+ext
	i = 0
	while not name_satisfied:
		if not os.path.exists(full_name):
			name_satisfied = True
		else: 
			full_name = join(dir, name)+'_'+str(i)+ext
			i+=1
	return full_name
	
########################################
### ARGUMENTS
########################################

run_directory = os.getcwd()
argp = argparse.ArgumentParser()

# # Are we running in pymol? If so, shift the input args
# # Is this a hack? Feels like a hack
# input_command = sys.argv[0]
# if "parquet2pymol" not in input_command:
# 	for i in range(len(sys.argv)-1):
# 		sys.argv[i] = sys.argv[i+1]
# 	sys.argv.pop()

argp.add_argument('--input-dir', '-i', default = ".", help="Parquet input directory, relative to working directory.")
argp.add_argument('--output-dir', '-o', default=None, help="Output PDBs in a persistent directory (default: deletes files after running)")

args = argp.parse_args()

residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
			'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
			'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
			'SER', 'THR', 'TRP', 'TYR', 'VAL']

def main():
	### Make Dataframes from parquest files
	input_dir = args.input_dir
	df_list = []
	for resid in residues:
		parquet_name = resid+'.parquet.gzip'
		if os.path.exists(join(input_dir, parquet_name)):
			df_list.append((resid, pd.read_parquet(join(input_dir, parquet_name))))
		else:
			print('No file for '+resid)

	#print(list(df_list[0][1].columns))

	#open_outdir =  open(args.output_dir) if args.output_dir else tempfile.TemporaryDirectory()
	with tempfile.TemporaryDirectory() as pdbdir:
		#If user specified a persistent output directory, use that:
		if args.output_dir:
			pdbdir = args.output_dir
			if not os.path.exists(pdbdir):
				os.mkdir(pdbdir)

		### Create atom groups, so we can order and output as pdbs later
		print('Creating atom groups...')
		atomgroup_list = []
		for resid, data in df_list:
			last_chain = None
			last_idx = 0
			for i in range(len(data)):
				if (data.loc[i, 'chain'] == 'Y' and last_chain == 'X') or (i == len(data)-1):
					#grab the last idx and current idx, use this to generate atomgroup and pdb
					atom_group = create_atomgroup(data[last_idx:i])
					atomgroup_list.append((atom_group, resid, data.loc[last_idx, 'C_score_bb_ind'], data.loc[last_idx, 'centroid'], data.loc[last_idx, 'pdb_name'], data.loc[last_idx, 'cluster_number']))
					last_idx = i
				last_chain = data.loc[i, 'chain']

		print('Sorting atom groups...')
		### Sort list of atom groups by C_score
		atomgroup_list.sort(key=lambda a: a[2], reverse=True)

		### Make PDB files and record names to load into pymol:
		print('Writing PDBs...')
		pdb_file_names = []
		for atom_group, resid, c_score, is_centroid, pdb, cnum in atomgroup_list:
			#Here, we're only showing the centroids.
			if is_centroid == True:
				pdb_file_name = generate_filename(pdbdir,resid+'_'+str(c_score)[:5]+'_'+pdb+'_clus'+str(cnum),'_centroid.pdb')
				if os.path.exists(pdb_file_name):
					print("Warning: centroid file already exists. Overwriting!!!")
				prody.writePDB(pdb_file_name, atom_group)
				pdb_file_names.append(pdb_file_name)
			else:
				pdb_file_name = generate_filename(pdbdir,resid+'_'+str(c_score)[:5]+'_'+pdb+'_clus'+str(cnum),'.pdb')
				if os.path.exists(pdb_file_name):
					print("Warning: file already exists. Overwriting!!!")
				prody.writePDB(pdb_file_name, atom_group)
				pdb_file_names.append(pdb_file_name)

if __name__ == '__main__':
	main()