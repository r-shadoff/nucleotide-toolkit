from collections import defaultdict

class Sequence:
    def __init__(self, sequence: str, identifier: str = None):
        self.identifier = identifier or "Unnamed"
        self.sequence = sequence.upper()
        self.seq_type =  self._detect_type()

        if not self._is_valid():
            raise ValueError("Invalid characters in sequence!")

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
 
    def dna_to_rna(self):
        """
        Function: Transcribes DNA to RNA
        Returns: RNA Sequence 
        """
        if self.seq_type != "DNA":
            raise ValueError("Input sequence must be DNA.")
        
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
        """
        Function: Transcribes RNA to cDNA
        Returns: cDNA sequence 
        """
        if self.seq_type != "RNA":
            raise ValueError("Input sequence must be RNA.")
        
        dna_dict = {
            "a":"A",
            "c":"C",
            "u":"T",
            "g":"G"
        }
        sequence_lower = self.sequence.lower()
        dna_table = str.maketrans(dna_dict)
        dna = sequence_lower.translate(dna_table)
        dna_upper = dna.upper()
        return(dna_upper)
    
    def rev_complement(self, sequence):
        """
        Function: Generates the reverse complement of the sequence
        Returns: Reverse complement
        """
        if self.seq_type != "DNA":
            raise ValueError("Input sequence must be DNA.")
        

        complement_dict = {
            "A":"T",
            "C":"G",
            "T":"A",
            "G":"C"
        }

        complement = sequence.upper().translate(str.maketrans(complement_dict))
        rev_complement = complement[::-1]
        return(rev_complement)  

    def gc_content(self):
        """
        Function: Computes GC content. Only considers G and C nucleotides. 
          Ambiguous nucleotides not included in GC content, but are included in the length of the sequence. 
          Ambiguous nucleotides: Y, S, K, M, B, D, H, V, N
        Returns: GC content of the sequence as a percentage.
        """
        counts = {
        "C" : 0,
        "G" : 0
        }
        for nucleotide in self.sequence:
            if nucleotide == "C":
                counts["C"] = counts["C"] +1
            elif nucleotide == "G":
                counts["G"] = counts["G"] +1
            else:
                pass
        
        gc_content = (counts["G"] + counts["C"]) / len(self.sequence) * 100
        return(gc_content)
    
    def rna_to_protein(self):
        if self.seq_type != "RNA":
            raise ValueError("Input sequence must be RNA.")
            return None
        
        codon_dict = {
            "UUU":"F",
            "CUU":"L",
            "AUU":"I",
            "GUU":"V",
            "UUC":"F",
            "CUC":"L",
            "AUC":"I",
            "GUC":"V",
            "UUA":"L",
            "CUA":"L",
            "AUA":"I",
            "GUA":"V",
            "UUG":"L",
            "CUG":"L",
            "AUG":"M",
            "GUG":"V",
            "UCU":"S",
            "CCU":"P",
            "ACU":"T",
            "GCU":"A",
            "UCC":"S",
            "CCC":"P",
            "ACC":"T",
            "GCC":"A",
            "UCA":"S",
            "CCA":"P",
            "ACA":"T",
            "GCA":"A",
            "UCG":"S",
            "CCG":"P",
            "ACG":"T",
            "GCG":"A",
            "UAU":"Y",
            "CAU":"H",
            "AAU":"N",
            "GAU":"D",
            "UAC":"Y",      
            "CAC":"H",      
            "AAC":"N",      
            "GAC":"D",
            "UAA":"Stop",
            "CAA":"Q",
            "AAA":"K",
            "GAA":"E",
            "UAG":"Stop",
            "CAG":"Q",
            "AAG":"K",
            "GAG":"E",
            "UGU":"C",
            "CGU":"R",
            "AGU":"S",
            "GGU":"G",
            "UGC":"C",
            "CGC":"R",
            "AGC":"S",
            "GGC":"G",
            "UGA":"Stop",   
            "CGA":"R",
            "AGA":"R",      
            "GGA":"G",
            "UGG":"W",      
            "CGG":"R",      
            "AGG":"R",
            "GGG":"G"
            }

        protein = ""
        for i in range(0, len(self.sequence), 3):
            codon = self.sequence[i:i+3]
            aa = codon_dict.get(codon, '')
            if aa == "Stop":
                # protein += aa
                break
            else:
                protein += aa
        return(protein)

    def count_snps(self, other):
        len1 = len(self.sequence)
        len2 = len(other.sequence)
        snps = 0
        if len1 == len2:
            for a, b in zip(self.sequence, other.sequence): 
                if a != b:
                    snps += 1
        else:
            print(f"Sequences must be the same length! Sequence 1 length: {len1} Sequence 2 length: {len2}")
            return None

        return snps

    def find_motif(self, motif):
        match_positions = []
        start = 0

        while start <= len(self.sequence) - len(motif):
            index = self.sequence.find(motif, start)
            if index == -1: # indicates no match
                break
            match_positions.append(index + 1) 
            start = index + 1  

        return match_positions
                
    def restriction_sites(self):
        """
        Function: Locates substrings within the sequence and evaluates if the substrings constitute restriction enzyme sites
        (Restriction enzymes cut at sites that are reverse palindromes by binding the enzyme homodimer to the parent and complement strand)
        Returns: List of the position and length of each potential restriction site
        """        
        sites = [] # tuple of starting position and length for each site
        
        for number in range(4,13,1):
            for start_pos in range(len(self.sequence) - number +1):
                subseq = self.sequence[start_pos:start_pos + number]
                rev_comp = self.rev_complement(subseq)
                if subseq == rev_comp:
                    sites.append((start_pos+1, number))
                
        for start,length in sites:
            print(start, length)
        return sites
