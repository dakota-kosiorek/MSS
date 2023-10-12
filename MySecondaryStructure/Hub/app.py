from MySecondaryStructure.PDBparser import ReadPDB
from MySecondaryStructure.PDBparser import Validate

def run(pdb_file_path):
    """
    About
    ---------------
    Runs the MySecondaryStructure application
    
    Parameters
    ---------------
    pdb_file_path : str
        Path to a pdb file
    """
    valid_file, e = Validate.check_file(pdb_file_path)
    
    if valid_file == True:
        protein = ReadPDB.read(pdb_file_path)
        protein.cmd_about()
        
    else:
        print(f"ERROR: {e} ({pdb_file_path})")