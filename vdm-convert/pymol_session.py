from pymol import cmd

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
	

def create_session(pdb_file_names):
    if len(pdb_file_names) > 1000:
        print("Warning: more than 1000 structures to load. Only showing top structures")
    i = 0
    for pdb in pdb_file_names:
        i+=1
        if i > 1000:
            break
        cmd.load(pdb)
        #cmd.disable("all") ###

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