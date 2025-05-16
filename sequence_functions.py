from collections import defaultdict
class Sequence:
    def __init__(self, sequence: str, identifier: str = None):
        self.identifier = identifier or "Unnamed"
        self.sequence = sequence.upper()
        self.seq_type = self._detect_type()

        if not self._is_valid():
            raise ValueError("Invalid characters in sequence!")
    
    #TO DO: Add protein sequence type
    def _detect_type(self):
        """
        Function: Assigns sequence type based on nucleotides in sequence.
        Returns: DNA or RNA as seq_type.
        """
        if "U" in self.sequence and "T" in self.sequence:
            raise ValueError("Mixed T and U bases in sequence!")
        elif "T" in self.sequence:
            return "DNA"
        elif "U" in self.sequence:
            return "RNA"
        else:
            raise ValueError("Sequence type unable to be determined.")  
    
    def _is_valid(self):
        """
        Function: Checks if the input sequence is valid.
        Returns: Boolean
        """
        valid_dna = set("ACGTRYSWKMBDHVN-")
        valid_rna = set("ACGURYSWKMBDHVN-")
        if self.seq_type == "DNA":
            return set(self.sequence).issubset(valid_dna)
        elif self.seq_type == "RNA":
            return set(self.sequence).issubset(valid_rna)
        else:
            return False
    
    def count_nucleotides(self):
        """ 
        Function: Counts the number of times IUPAC nucleotides appear in a DNA sequence
        Returns: a dictionary of nucleotide counts 
        """
        if self.seq_type == "DNA":
            iupac = list("ACGTYRWSKMBDHVN-")
        elif self.seq_type == "RNA":
            iupac = list("ACGUYRWSKMBDHVN-")
        else:
            print("Unsupported sequence type for this function.")
    
        count = defaultdict(int)
        for nucleotide in self.sequence:
            if nucleotide in iupac:
                count[nucleotide] += 1
            else:
                raise ValueError(f"Unknown nucleotide: {nucleotide}")
        
        return dict(count)
    
    def gc_content(self):
        pass
    
    def dna_to_rna(self):
        """
        Function: Transcribes DNA to RNA
        Returns: RNA Sequence 
        """
        rna_dict = {
            "a":"A",
            "c":"C",
            "t":"U",
            "g":"G"
        }
        sequence_lower = self.sequence.lower()
        rna_table = str.maketrans(rna_dict)
        rna = sequence_lower.translate(rna_table)
        rna_upper = rna.upper()
        return(rna_upper)

    
    def rna_to_dna(self):
        pass
    
    def rna_to_protein(self):
        pass
    
    def complement(self):
        pass
    

        