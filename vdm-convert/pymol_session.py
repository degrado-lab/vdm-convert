from pymol import cmd
from os import listdir
from os.path import isfile, join

#######################################
### HELPER FUNCTIONS
#######################################
def show_obj(i):
	current_obj = "all" if i == -1 else cmd.get_object_list('(all)')[i]
	cmd.enable(current_obj)

def hide_obj(i):
	current_obj = "all" if i == -1 else cmd.get_object_list('(all)')[i]
	cmd.disable(current_obj)
	
def next_obj():
	hide_obj(cmd.idx)
	cmd.idx+=1
	if cmd.idx > len(cmd.get_object_list('(all)'))-1:
		cmd.idx = -1
	show_obj(cmd.idx)

def prev_obj():
	hide_obj(cmd.idx)
	cmd.idx-=1
	if cmd.idx < -1:
		cmd.idx = len(cmd.get_object_list('(all)'))-1
	show_obj(cmd.idx)
	

def create_session(input_dir, input_type='PDB', residues=None, score_cutoff=None):
	'''
	Creates a pymol session from a directory of converted vdM files.
	Arguments:
		input_dir: directory containing PDB files
		input_type: file type of input files (default: PDB)
		residues: list of residues to include (default: all)
		score_cutoff: minimum score to include (default: None)
	Returns:
		None
	'''

	if residues is None:
		residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
			'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
			'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
			'SER', 'THR', 'TRP', 'TYR', 'VAL']
	
	input_file_names = [join(input_dir, f) for f in listdir(input_dir) if isfile(join(input_dir, f)) and f.endswith(input_type)]

	#Get scores from each infile name:
	scores = []
	for f in input_file_names:
		scores.append(float(f.split('_')[1]))
	#Get residues from each infile name:
	residues = []
	for f in input_file_names:
		residues.append(f.split('_')[0])
	
	#sort input file names by score:
	input_file_names = [x for _,x in sorted(zip(scores,input_file_names))]
	residues = [x for _,x in sorted(zip(scores,residues))]

	if len(input_file_names) > 1000:
		print("Warning: more than 1000 structures to load. Only showing top structures")
	
	for i, pdb in enumerate(input_file_names):
		if i > 1000:
			break
		if score_cutoff is not None and scores[i] > score_cutoff:
			if residues[i] in residues:
				cmd.load(pdb)

	cmd.hide("all")
	cmd.show_as("licorice", "chain X")
	cmd.color("atomic", "chain X")

	cmd.show_as("spheres", "chain Y")
	cmd.color("atomic", "chain Y")
	cmd.set("sphere_transparency", 0.8, "chain Y")

	### Set up iterating over structures
	cmd.idx = 0

	cmd.set_key("right",next_obj)
	cmd.set_key("left",prev_obj)

	#Hide all but first object, set up scene:
	cmd.disable("(all)")
	show_obj(cmd.idx)