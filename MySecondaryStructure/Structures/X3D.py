import pandas as pd

class Protein():
    """
    About
    ---------------
    This class contains all relevant data gathered from a PDB file.
    """
    def __init__(self, atom_df: pd.DataFrame, helix_df: pd.DataFrame, sheet_df: pd.DataFrame, ssbond_df: pd.DataFrame, 
                   conect_df: pd.DataFrame, het_df: pd.DataFrame, hetnam_df: pd.DataFrame, hetsyn_df: pd.DataFrame, 
                   formul_df: pd.DataFrame, link_df: pd.DataFrame, cispep_df: pd.DataFrame, seqres_df: pd.DataFrame,
                   name: str):
        self.atom_df = atom_df
        self.helix_df = helix_df
        self.sheet_df = sheet_df
        self.ssbond_df = ssbond_df
        self.conect_df = conect_df
        self.het_df = het_df
        self.hetnam_df = hetnam_df
        self.hetsyn_df = hetsyn_df
        self.formul_df = formul_df
        self.link_df = link_df
        self.cispep_df = cispep_df
        self.seqres_df = seqres_df
        self.name = name
        
        # Number of atoms in protein residues
        self.num_atoms = len(self.atom_df.loc[self.atom_df["Record Type"] == "ATOM"]["Atom Serial Number"].unique())
        # Number of alpha helices
        self.num_helices = 0
        if len(self.helix_df) > 0:
            self.num_helices = len(self.helix_df["Helix Identifier"].unique())
        # Number of Beta sheets
        self.num_sheets = 0
        if len(self.sheet_df) > 0:
            # Adds togehter the number of sheets found in each chain that have a unique starting residue number
            for c in self.sheet_df["1st Chain Identifier"].unique():
                self.num_sheets += len(self.sheet_df.loc[self.sheet_df["1st Chain Identifier"] == c]["1st Residue Sequence Number"].unique())
        # Number of loops
        self.num_loops, self.chain_loops = self.get_loops()
        
        # ---------------------------------------------------------------------
    
    def get_loops(self):
        num_loops = 0
        chain_loops = list()
        
        # This list contains all of the chains in the protein that has the 
        # largest number of protein secondary structures
        highest_chain = list()
        highest_chain_df = pd.DataFrame()
        lowest_chain = list()
        lowest_chain_df = pd.DataFrame()
        
        helix_chains = list()
        sheet_chains = list()
        
        # If there are alpha helicies in the protein
        if len(self.helix_df) > 0:
            helix_chains = list(self.helix_df["1st Chain Identifier"].unique())
            highest_chain = helix_chains
            highest_chain_df = self.helix_df
        # If there are beta sheets in the protein
        if len(self.sheet_df) > 0:
            sheet_chains = list(self.sheet_df["1st Chain Identifier"].unique())
            lowest_chain = sheet_chains
            lowest_chain_df = self.sheet_df
        
        # If there are more beta sheets than alpha helicies
        # Switch around highest and lowest chain lists and data frames
        if len(helix_chains) < len(sheet_chains):
            highest_chain = sheet_chains
            highest_chain_df = self.sheet_df
            lowest_chain = helix_chains
            lowest_chain_df = self.helix_df

        # Loops through all the protein chains that contain secondary structures (that are not just loops)
        for c in highest_chain:
            num_chain_loops = 0
            
            # Start and end positions of secondary structures (alpha helicies and beta sheets)
            highest_chain_secondary_structure_landmarks = list()
            lowest_chain_secondary_structure_landmarks = list()
            secondary_structure_landmarks = list()
            
            structure_start = highest_chain_df.loc[highest_chain_df["1st Chain Identifier"] == c]["1st Residue Sequence Number"].to_list()
            structure_end = highest_chain_df.loc[highest_chain_df["1st Chain Identifier"] == c]["2nd Residue Sequence Number"].to_list()
            
            # A list of the start and end positions of a secondary structure with 
            highest_chain_secondary_structure_landmarks = sorted(list(set(structure_start + structure_end)))

            if c in lowest_chain:
                structure_start = lowest_chain_df.loc[lowest_chain_df["1st Chain Identifier"] == c]["1st Residue Sequence Number"].to_list()
                structure_end = lowest_chain_df.loc[lowest_chain_df["1st Chain Identifier"] == c]["2nd Residue Sequence Number"].to_list()
                
                lowest_chain_secondary_structure_landmarks = sorted(list(set(structure_start + structure_end)))
            
            secondary_structure_landmarks = sorted(list(set(highest_chain_secondary_structure_landmarks + lowest_chain_secondary_structure_landmarks)))
            
            total_chain_amino_acids = self.seqres_df.loc[self.seqres_df["Chain Identifier"] == c]["Number of Residuals"].unique()[0]
            
            if len(secondary_structure_landmarks) > 0:
                if secondary_structure_landmarks[0] > 0:
                    num_chain_loops += 1
                if secondary_structure_landmarks[-1] < total_chain_amino_acids:
                    num_chain_loops += 1
                
                """    
                for i in range(0, len(secondary_structure_landmarks), 2):
                    print(f"{c}\t{secondary_structure_landmarks[i]}\t-->\t{secondary_structure_landmarks[i + 1]}")
                
                """    
                for i in range(1, len(secondary_structure_landmarks)-1, 2):
                    if secondary_structure_landmarks[i] != secondary_structure_landmarks[i+1]-1:
                        num_chain_loops += 1
            
            chain_loops.append(num_chain_loops)
        
        # If there are beta sheets of alpha helicies then each chain is comprised of one loop
        if len(self.sheet_df) == 0 or len(self.helix_df) == 0:
            for c in self.atom_df["Chain Identifier"].unique():
                chain_loops.append(1)
        
        num_loops = sum(chain_loops)
        return num_loops, chain_loops
    
    def cmd_about(self):
        total_num_secondary_structures = self.num_loops + self.num_helices + self.num_sheets
        
        # Protein name
        print(f"------------------- {self.name} -------------------")
        # Total number of amino acids
        amino_acids_per_chain = {}
        total_amino_acids = 0
        
        for c in self.seqres_df["Chain Identifier"].unique():
            chain_amino_acids = self.seqres_df.loc[self.seqres_df["Chain Identifier"] == c]["Number of Residuals"].unique()[0]
            amino_acids_per_chain[c] = chain_amino_acids
            total_amino_acids += chain_amino_acids
        
        data = {
            "Count": [self.num_loops, self.num_helices, self.num_sheets],
            "% of Structures": [round((self.num_loops/total_num_secondary_structures) * 100, 2), 
                                          round((self.num_helices/total_num_secondary_structures) * 100, 2),
                                          round((self.num_sheets/total_num_secondary_structures) * 100, 2)
                                        ]
        }
        
        df = pd.DataFrame(data, index=["Loop", "Alpha Helix", "Beta Sheet"])
        
        print(f"Total amino aicds: {total_amino_acids}")
        print(f"Number of Chains: {len(self.seqres_df['Chain Identifier'].unique())}\n")
        print(df,"\n")
        
        # All of this info broken down per chain
        #print(amino_acids_per_chain)
        