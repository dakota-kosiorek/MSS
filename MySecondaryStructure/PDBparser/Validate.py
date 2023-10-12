from os.path import exists

def check_file(file_path: str):
    """
    About
    ---------------
    Checks to see if the file is an existing PDB file
    
    Parameters
    ---------------
    file_path : str
        Path to a file
        
    Returns
    ---------------
    bool
        True if the 'file_path' parameter leads to a real PDB file, false if otherwise
    str
        Message containing information on what was wrong with 'file_path'
    """
    
    valid_file = exists(file_path)
    
    if valid_file == False:
        return valid_file, "File does not exist!"
    
    valid_file = file_path.lower().endswith(".pdb")
    
    if valid_file == False:
        return valid_file, "File is not in PDB format!"
    
    return valid_file, "File is all good!"
    