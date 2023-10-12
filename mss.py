# Dakota Kosiorek

# To run: 
#   python mss.py [file path]

from MySecondaryStructure.Hub import app
import sys

def main(pdb_file_path):
    app.run(pdb_file_path)

if __name__ == "__main__":
    for argument in sys.argv[1:]:
        main(argument)