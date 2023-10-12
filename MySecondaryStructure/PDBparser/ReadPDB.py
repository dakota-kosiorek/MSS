import pandas as pd

from MySecondaryStructure.Structures.X3D import Protein

def read(file_path: str) -> Protein:
    """
    About
    ---------------
    Reads a PDB file and formats the data in a 'Protein' object.
    
    Parameters
    ---------------
    file_path : str
        The path to the PDB file
    
    Return
    ---------------
    Protein
        A Protein object containing the data from the PDB file found at filepath
    """
    
    raw_PDB = get_raw_PDB(file_path)
    my_protein = process_raw(raw_PDB)
    
    return my_protein

# Gets all data from a pdb file and puts it into a string       
def get_raw_PDB(file_path: str) -> str:
    """
    About
    ---------------
    Read's a file and returns all contents of that file as a string.
    
    Parameters
    ---------------
    raw_PDB : str
        Entire PDB file in a single variable
    
    Returns
    ---------------
    str
        A string containing the contents of the PDB file
    """
    
    raw_PDB = ""
    with open(file_path, "r") as file:
        for line in file.readlines():
            raw_PDB += line
    
    return raw_PDB

def process_raw(raw_PDB: str) -> Protein:    
    """
    About
    ---------------
    Turns PDB data stored in a string into a 'Protein' object.
    
    Parameters
    ---------------
    raw_PDB : str
        Entire PDB file in a single variable
        
    Returns
    ---------------
    Protein
        A protein object filled with data from the raw PDB file
    """
    
    atom_list = []
    helix_list = []
    sheet_list = []
    ssbond_list = []
    conect_list = []
    het_list = []
    hetnam_list = []
    hetsyn_list = []
    formul_list = []
    link_list = []
    cispep_list = []
    seqres_list = []
    
    name = ""
    
    # Line elements are split according to https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    # Get atom positions and secondary structure information from the PDB file
    for line in raw_PDB.split("\n"):
        # Atom information and Hetatm information
        # Hetatm being atoms that are in non-standard residues
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            atom_type = "ATOM"
            if line[:6] == "HETATM":
                atom_type = "HETATM"
            
            line_dict = {
                "Record Type": atom_type,
                "Atom Serial Number": line[6:11].strip(),
                "Atom Name": line[12:16].strip(), 
                "Alternate Location Indicator": line[16].strip(), 
                "Residue Name": line[17:20].strip(), 
                "Chain Identifier": line[21].strip(), 
                "Residue Sequence Number": line[22:26].strip(), 
                "Code for Insertions of Residues": line[26].strip(),
                "X orthogonal Coordinate": line[30:38].strip(), 
                "Y orthogonal Coordinate": line[38:46].strip(), 
                "Z orthogonal Coordinate": line[46:54].strip(), 
                "Occupancy": line[54:60].strip(), 
                "Temperature Factor": line[60:66].strip(), 
                "Segment Identifier": line[72:76].strip(), 
                "Element Symbol": line[76:78].strip(), 
                "Charge": line[78:80].strip()
            }
            
            atom_list.append(line_dict)
            
        elif line[:5] == "HELIX":
            line_dict = {
                "Record Type": "HELIX",
                "Helix Serial Number": line[7:10].strip(),
                "Helix Identifier": line[11:14].strip(),
                "Initial Residue Name": line[15:18].strip(),
                "1st Chain Identifier": line[19].strip(),
                "1st Residue Sequence Number": line[21:25].strip(),
                "1st Code for Insertions of Residues": line[25].strip(),
                "Terminal Residue Name": line[27:30].strip(),
                "2nd Chain Identifier": line[31].strip(),
                "2nd Residue Sequence Number": line[33:37].strip(),
                "2nd Code for Insertions of Residues": line[37].strip(),
                "Type of Helix": line[38:40].strip(),
                "Comment": line[40:70].strip(),
                "Length of Helix": line[71:76].strip()
            }

            helix_list.append(line_dict)
            
        elif line[:5] == "SHEET":
            line_dict = {
                "Record Type": "SHEET",
                "Strand Number (Current Sheet)": line[7:10].strip(),
                "Sheet Identifier": line[11:14].strip(),
                "Number of Strands (Current Sheet)": line[14:16].strip(),
                "Initial Residue Name": line[17:20],
                "1st Chain Identifier": line[21].strip(),
                "1st Residue Sequence Number": line[22:26].strip(),
                "1st Code for Insertions of Residues": line[26].strip(),
                "Terminal Residue Name": line[28:31].strip(),
                "2nd Chain Identifier": line[32].strip(),
                "2nd Residue Sequence Number": line[33:37].strip(),
                "2nd Code for Insertions of Residues": line[37].strip(),
                "Strand Sense with Respect to Previous": line [38:40].strip(),
                # The following fields identify two atoms invovled in a 
                # hydrogen bond
                "HB 1st Atom Name": line[41:45].strip(),
                "HB 1st Residue Name": line[45:48].strip(),
                "HB 1st Chain Identifier": line[49].strip(),
                "HB 1st Residue Sequence Number": line[50:54].strip(),
                "HB 1st Code for Insertions of Residues": line[54].strip(),
                "HB 2nd Atom Name": line[56:60].strip(),
                "HB 2nd Residue Name": line[60:63].strip(),
                "HB 2nd Chain Identifier": line[64].strip(),
                "HB 2nd Residue Sequence Number": line[65:69].strip(),
                "HB 2nd Code for Insertions of Residues": line[69].strip()
            }
            
            sheet_list.append(line_dict)
            
        elif line[:6] == "SSBOND":
            line_dict = {
                "Record Type": "SSBOND",
                "Serial Number": line[7:10].strip(),
                "1st Residue Name": line[11:14].strip(),
                "1st Chain Identifier": line[15].strip(),
                "1st Residue Sequence Number": line[17:21].strip(),
                "1st Code for Insertions of Residues": line[21].strip(),
                "2nd Residue Name": line[25:28].strip(),
                "2nd Chain Identifier": line[29].strip(),
                "2nd Residue Sequence Number": line[31:35].strip(),
                "2nd Code for Insertions of Residues": line[35].strip(),
                "Symmetry Operator for 1st Residue": line[59:65].strip(),
                "Symmetry Operator for 2nd Residue": line[66:72].strip(),
                "Length of Disulfide Bond": line[73:78].strip()
            }

            ssbond_list.append(line_dict)
            
        elif line[:6] == "CONECT":
            line_dict = {
                "Record Type": "CONECT",
                "Atom Serial Number": line[6:11].strip(),
                "Bonded Atom 1 Serial Number": line[11:16].strip(),
                "Bonded Atom 2 Serial Number": line[16:21].strip(),
                "Bonded Atom 3 Serial Number": line[21:26].strip(),
                "Bonded Atom 4 Serial Number": line[26:31].strip(),
                "Hydrogen Bonded Atom 1 Serial Number": line[31:36].strip(),
                "Hydrogen Bonded Atom 2 Serial Number": line[36:41].strip(),
                "Salt Bridged Atom 1 Serial Number": line[41:46].strip(),
                "Hydrogen Bonded Atom 3 Serial Number": line[46:51].strip(),
                "Hydrogen Bonded Atom 4 Serial Number": line[51:56].strip(),
                "Salt Bridged Atom 2 Serial Number": line[56:61].strip()
            }
            
            conect_list.append(line_dict)
    
        elif line[:6] == "HET   ":
            line_dict = {
                "Record Type": "HET",
                "Het ID": line[7:10].strip(),
                "Chain ID": line[12].strip(),
                "Sequence Number": line[13:17].strip(),
                "Insertion Code": line[17].strip(),
                "Number of HETATM Atoms": line[20:25].strip(),
                "Description": line[30:70].strip()
            }
            
            het_list.append(line_dict)
       
        elif line[:6] == "HETNAM":
            line_dict = {
                "Record Type": "HETNAM",
                "Continuation": line[8:10].strip(),
                "Het ID": line[11:14].strip(),
                "Chemical Name": line[15:70].strip()
            } 
            
            hetnam_list.append(line_dict)
        
        elif line[:6] == "HETSYN":
            line_dict = {
                "Record Type": "HETSYN",
                "Continuation": line[8:10].strip(),
                "Het ID": line[11:14].strip(),
                "Het Synonyms": line[15:70].strip()
            }
            
            hetsyn_list.append(line_dict)
            
        elif line[:6] == "FORMUL":
            line_dict = {
                "Record Type": "FORMUL",
                "Component Number": line[8:10].strip(),
                "Het ID": line[12:15].strip(),
                "Continuation Number": line[16:18].strip(),
                "Asterisk (Water)": line[18].strip(),
                "Chemical Formula": line[19:70].strip()
            }
            
            formul_list.append(line_dict)
                
        elif line[:4] == "LINK":
            line_dict = {
                "Record Type": "LINK",
                "1st Atom Name": line[12:16].strip(),
                "1st Alternate location indicator": line[16].strip(),
                "1st Residue Name": line[17:20].strip(),
                "1st Chain Identifier": line[21].strip(),
                "1st Residue Sequence Number": line[22:26].strip(),
                "1st Insertion Code": line[26].strip(),
                "2nd Atom Name": line[42:46].strip(),
                "2nd Alternate location indicator": line[46].strip(),
                "2nd Residue Name": line[47:50].strip(),
                "2nd Chain Identifier": line[51].strip(),
                "2nd Residue Sequence Number": line[52:56].strip(),
                "2nd Insertion Code": line[56].strip(),
                "1st Atom Symmetry Operator": line[59:65].strip(),
                "2nd Atom Symmetry Operator": line[66:72].strip(),
                "Link Distance": line[73:78].strip()
            }
            
            link_list.append(line_dict)
            
        elif line[:6] == "CISPEP":
            line_dict = {
                "Record Type": "CISPEP",
                "Record Serial Number": line[7:10].strip(),
                "1st Residue Name": line[11:14].strip(),
                "1st Chain Identifier": line[15].strip(),
                "1st Residue Sequence Number": line[17:21].strip(),
                "1st Insertion Code": line[21].strip(),
                "2nd Residue Name": line[25:28].strip(),
                "2nd Chain Identifier": line[29].strip(),
                "2nd Residue Sequence Number": line[31:35].strip(),
                "2nd Insertion Code": line[35].strip(),
                "Specific Model": line[43:46].strip(),
                "Angle Measurment (Deg)": line[53:59].strip()
            }
            
            cispep_list.append(line_dict)
            
        elif line[:6] == "SEQRES":
            line_dict = {
                "Record Type": "SEQRES",
                "Serial Number": line[7:10].strip(),
                "Chain Identifier": line[11].strip(),
                "Number of Residuals": line[13:17].strip(),
                "1st Residual Name": line[21:22].strip(),
                "2nd Residual Name": line[23:26].strip(),
                "3rd Residual Name": line[27:30].strip(),
                "4th Residual Name": line[31:34].strip(),
                "5th Residual Name": line[35:38].strip(),
                "6th Residual Name": line[39:42].strip(),
                "7th Residual Name": line[43:46].strip(),
                "8th Residual Name": line[47:50].strip(),
                "9th Residual Name": line[51:54].strip(),
                "10th Residual Name": line[55:58].strip(),
                "11th Residual Name": line[59:62].strip(),
                "12th Residual Name": line[63:66].strip(),
                "13th Residual Name": line[67:70].strip(),
            }
            
            seqres_list.append(line_dict)
            
        elif line[:6] == "HEADER":
            name = line[62:66]
            
    atom_df_col_types = {
        # Everything else is str
        "Atom Serial Number":                       "Int",
        "Residue Sequence Number":                  "Int", 
        "X orthogonal Coordinate":                  "Float", 
        "Y orthogonal Coordinate":                  "Float", 
        "Z orthogonal Coordinate":                  "Float", 
        "Occupancy":                                "Float", 
        "Temperature Factor":                       "Float", 
    }
    
    helix_df_col_types = {
        # Everything else is str
        "Helix Serial Number":                      "Int",
        "1st Residue Sequence Number":              "Int",
        "2nd Residue Sequence Number":              "Int",
        "Type of Helix":                            "Int",
        "Length of Helix":                          "Int",
    }
    
    sheet_df_col_types = {
        # Everything else is str
        "Strand Number (Current Sheet)":            "Int",
        "Number of Strands (Current Sheet)":        "Int",
        "1st Residue Sequence Number":              "Int",
        "2nd Residue Sequence Number":              "Int",
        "Strand Sense with Respect to Previous":    "Int",
        "HB 1st Residue Sequence Number":           "Int",
        "HB 2nd Residue Sequence Number":           "Int",
    }
    
    ssbond_df_col_types = {
        # Everything else is str
        "Serial Number":                            "Int",
        "1st Residue Sequence Number":              "Int",
        "2nd Residue Sequence Number":              "Int",
        "Symmetry Operator for 1st Residue":        "Int",
        "Symmetry Operator for 2nd Residue":        "Int",
        "Length of Disulfide Bond":                 "Float"
    }
    
    conect_df_col_types = {
        # Everything else is str
        "Atom Serial Number":                       "Int",
        "Bonded Atom 1 Serial Number":              "Int",
        "Bonded Atom 2 Serial Number":              "Int",
        "Bonded Atom 3 Serial Number":              "Int",
        "Bonded Atom 4 Serial Number":              "Int",
        "Hydrogen Bonded Atom 1 Serial Number":     "Int",
        "Hydrogen Bonded Atom 2 Serial Number":     "Int",
        "Salt Bridged Atom 1 Serial Number":        "Int",
        "Hydrogen Bonded Atom 3 Serial Number":     "Int",
        "Hydrogen Bonded Atom 4 Serial Number":     "Int",
        "Salt Bridged Atom 2 Serial Number":        "Int"
    }
    
    het_df_col_types = {
        # Everything else is str
        "Sequence Number":                          "Int",
        "Number of HETATM Atoms":                   "Int"
    }
    
    hetnam_df_col_types = {
        # Everything else is str
    }
    
    hetsyn_df_col_types = {
        # Everything else is str
    }
    
    formul_df_col_types = {
        # Everything else is str
        "Component Number":                         "Int",
        "Continuation Number":                      "Int"
    }
    
    link_df_col_types = {
        # Everything else is str
        "1st Residue Sequence Number":              "Int",
        "2nd Residue Sequence Number":              "Int",
        #"1st Atom Symmetry Operator":               "Int",
        #"2nd Atom Symmetry Operator":               "Int",
        "Link Distance":                            "Float"
    }
    
    cispep_df_col_types = {
        # Everything else is str
        "Record Serial Number":                     "Int",
        "1st Residue Sequence Number":              "Int",
        "2nd Residue Sequence Number":              "Int",
        "Angle Measurment (Deg)":                   "Float"
    }
    
    seqres_df_col_types = {
        # Everything else is str
        "Serial Number":                            "Int",
        "Number of Residuals":                      "Int"
    }
    
    atom_df = pd.DataFrame.from_dict(atom_list)
    helix_df = pd.DataFrame.from_dict(helix_list)
    sheet_df = pd.DataFrame.from_dict(sheet_list)
    ssbond_df = pd.DataFrame.from_dict(ssbond_list)
    conect_df = pd.DataFrame.from_dict(conect_list)
    het_df = pd.DataFrame.from_dict(het_list)
    hetnam_df = pd.DataFrame.from_dict(hetnam_list)
    hetsyn_df = pd.DataFrame.from_dict(hetsyn_list)
    formul_df = pd.DataFrame.from_dict(formul_list)
    link_df = pd.DataFrame.from_dict(link_list)
    cispep_df = pd.DataFrame.from_dict(cispep_list)
    seqres_df = pd.DataFrame.from_dict(seqres_list)
    
    if len(atom_df) > 0:
        atom_df = change_col_type(atom_df, atom_df_col_types)
    if len(helix_df) > 0:
        helix_df = change_col_type(helix_df, helix_df_col_types)
    if len(sheet_df) > 0:
        sheet_df = change_col_type(sheet_df, sheet_df_col_types)
    if len(ssbond_df) > 0:
        ssbond_df = change_col_type(ssbond_df, ssbond_df_col_types)
    if len(conect_df) > 0:
        conect_df = change_col_type(conect_df, conect_df_col_types)
    if len(het_df) > 0:
        het_df = change_col_type(het_df, het_df_col_types)
    if len(hetnam_df) > 0:
        hetnam_df = change_col_type(hetnam_df, hetnam_df_col_types)
    if len(hetsyn_df) > 0:
        hetsyn_df = change_col_type(hetsyn_df, hetsyn_df_col_types)
    if len(formul_df) > 0:
        formul_df = change_col_type(formul_df, formul_df_col_types)
    if len(link_df) > 0:
        link_df = change_col_type(link_df, link_df_col_types)
    if len(cispep_df) > 0:
        cispep_df = change_col_type(cispep_df, cispep_df_col_types)
    if len(seqres_df) > 0:
        seqres_df = change_col_type(seqres_df, seqres_df_col_types)
            
    return Protein(atom_df, helix_df, sheet_df, ssbond_df, 
                   conect_df, het_df, hetnam_df, hetsyn_df, 
                   formul_df, link_df, cispep_df, seqres_df,
                   name
    )

def change_col_type(df: pd.DataFrame, col_dict: dict) -> pd.DataFrame:
    """
    About
    ---------------
    Changes the column types of a dataframe
    
    Parameters
    ---------------
    df : pd.DataFrame
        A pandas dataframe
    col_dict
        A dictionary with df columns as the key and 'Int' or 'Float' as the value
        
    Returns
    ---------------
    pd.DataFrame
        A dataframe with changed column types
    """
    
    for col_name, col_type in col_dict.items():
        if col_type == "Int" or col_type == "Float":
            df[col_name] = pd.to_numeric(df[col_name])
            
    return df