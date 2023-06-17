# vdm-convert
Python package for converting van der Mer (vdM) parquet.gzip files to other file types. 

## Installation

First, clone the directory anywhere on your machine using the following command:

`git clone https://github.com/njf042/vdm-convert.git`

Then, install the package to your python or conda environment using:

`pip install -e vdm-convert/`

## Usage

After installation, the following terminal commands are available:

### `vdm_convert`
Converts a vdm parquet.gzip file to a different file type.

Usage: 
` vdm_convert [-h] [--input-dir,-i INPUT_DIR] [--output-dir,-o OUTPUT_DIR]
                   [--output-type,-t OUTPUT_TYPE]`

Example: 
To convert all vdMs in the directory `vdms/` to PDB files, and store in the directory `converted_vdms/`:

`vdm_display -i vdms/ -o converted_vdms/`

### `vdm_display`
Launches a pymol session for easy viewing of converted vdMs.

Usage: ` vdm_display [-h] [--input-dir,-i INPUT_DIR] [--input-type,-t INPUT_TYPE]
                   [--score-cutoff,-s SCORE_CUTOFF]
                   [--residues,-r RESIDUES [RESIDUES ...]]`

Example: 
To display all vdMs with a score of 2.5 or greater, and only for the residues ALA, GLY, and LEU:

`vdm_display -i converted_vdms/ -s 2.5 -r ALA GLY LEU`

## PyMOL Tips:
PyMOL must be installed on your system (or environment) and accessible with the command `pymol` for the `vdm_display` command to work.

Once in the pymol session, the following commands are useful:
To change the look of the protein backbone:
`show_as [spheres, lines, sticks], chain X`

To change the look of the ligand:
`show_as [spheres, lines, sticks], chain Y`

To change the transparency of the ligand:
`set sphere_transparency, 0.75, chain Y`

## License
vdm-convert is released under the MIT License. See LICENSE for details.

## Contact
nicholas.freitas@ucsf.edu