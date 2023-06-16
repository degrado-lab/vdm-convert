import subprocess
import sys

def main():
    command_list = ["pymol_session.py"]+sys.argv[1:]
    subprocess.Popen(["pymol"]+command_list)

if __name__ == '__main__':
	main()