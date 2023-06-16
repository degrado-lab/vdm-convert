import numpy as np
import pandas as pd
import prody
import os
from os.path import join
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
### MAIN FUNCTIONS
########################################

def convert(input_dir, output_dir, out_type='PDB', residues=None):
	'''
	Converts parquet.gzip files to specified file type, and stores in output directory.
	Currently supports PDB, PQR, and XYZ.
	Arguments:
		input_dir: directory containing parquet.gzip files
		output_dir: directory to write output files to
		out_type: output file type (default: PDB)
		residues: list of residues to convert (default: all)
	Returns:
		None
	'''

	if residues is None:
		residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
			'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
			'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
			'SER', 'THR', 'TRP', 'TYR', 'VAL']
	
	### Make Dataframes from parquest files
	df_list = []
	for resid in residues:
		parquet_name = resid+'.parquet.gzip'
		if os.path.exists(join(input_dir, parquet_name)):
			df_list.append((resid, pd.read_parquet(join(input_dir, parquet_name))))
		else:
			print('No file for '+resid)

	#Create output directory if it doesn't exist
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	### Create atom groups, so we can order and output as different filetype later
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

	convert_function_dict = {'PDB': prody.writePDB, 'PQR': prody.writePQR, 'XYZ': prody.writeXYZ}

	### Make output files and record names to load into pymol:
	print('Writing '+out_type+'s...')
	out_file_names = []
	for atom_group, resid, c_score, is_centroid, pdb, cnum in atomgroup_list:
		#Here, we're only showing the centroids.
		if is_centroid == True:
			pdb_file_name = generate_filename(output_dir,resid+'_'+str(c_score)[:5]+'_'+pdb+'_clus'+str(cnum),'_centroid.pdb')
			if os.path.exists(pdb_file_name):
				print("Warning: centroid file already exists. Overwriting!!!")
			convert_function_dict[out_type](pdb_file_name, atom_group)
			#prody.writePDB(pdb_file_name, atom_group)
			out_file_names.append(pdb_file_name)
		else:
			pdb_file_name = generate_filename(output_dir,resid+'_'+str(c_score)[:5]+'_'+pdb+'_clus'+str(cnum),'.pdb')
			if os.path.exists(pdb_file_name):
				print("Warning: file already exists. Overwriting!!!")
			convert_function_dict[out_type](pdb_file_name, atom_group)
			out_file_names.append(pdb_file_name)

def main():
	argp = argparse.ArgumentParser()
	argp.add_argument('--input-dir', '-i', default = ".", help="Parquet input directory, relative to working directory.")
	argp.add_argument('--output-dir', '-o', default="convert", help="Output file directory")
	argp.add_argument('--output-type', '-t', default="PDB", help="Filetype to convert to. Currently supports PDB, PQR, and XYZ.")
	args = argp.parse_args()

	convert(args.input_dir, args.output_dir, out_type='PDB')


if __name__ == '__main__':
	main()