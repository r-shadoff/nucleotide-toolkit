class Sequence:
    def __init__(self, sequence: str, identifier: str = None):
        self.identifier = identifier or "Unnamed"
        self.sequence = sequence.upper()
        self.seq_type = self._detect_type()

        if not self._is_valid():
            raise ValueError("Invalid characters in sequence!")
        
    def _detect_type(self):
        if "U" in self.sequence and "T" in self.sequence:
            raise ValueError("Mixed T and U bases in sequence!")
        elif "T" in self.sequence:
            return "DNA"
        elif "U" in self.sequence:
            return "RNA"
        else:
            raise ValueError("Sequence type unable to be determined.")
        
    def _is_valid(self):
        valid_dna = set("ACGT")
        valid_rna = set("ACGU")
        if self.seq_type == "DNA":
            return set(self.sequence).issubset(valid_dna)
        elif self.seq_type == "RNA":
            return set(self.sequence).issubset(valid_rna)
        else:
            return False
    
    def count_nucleotides(self):
        pass
    
    def gc_content(self):
        pass
    
    def dna_to_rna(self):
        pass
    
    def rna_to_dna(self):
        pass
    
    def rna_to_protein(self):
        pass
    
    def complement(self):
        pass
    

        