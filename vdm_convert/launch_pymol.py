import subprocess
import sys
import argparse

def main():
	argp = argparse.ArgumentParser()
	argp.add_argument('--input-dir', '-i', default = ".", help="Parquet input directory, relative to working directory.")
	argp.add_argument('--input-type', '-t', default="PDB", help="Filetype of converted vdMs. Currently supports PDB, PQR, and XYZ.")
	argp.add_argument('--score-cutoff', '-s', default=None, help="Minimum score to include. Default: include all.")
	argp.add_argument('--residues', '-r', default=None, nargs='+', help="Residues to include. Default: include all.")
	args = argp.parse_args()

	run_string = "vdm_convert.pymol_session.create_session(\'{}\', \'{}\'".format(args.input_dir, args.input_type)
	if args.residues is not None:
		run_string += ", selected_residues={}".format(args.residues)
	if args.score_cutoff is not None:
		run_string += ", score_cutoff={}".format(args.score_cutoff)
	run_string += ")"

	subprocess.Popen(["pymol", "-d", "import vdm_convert; " + run_string])
	
if __name__ == '__main__':
	main()