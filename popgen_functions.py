def offspring_probability(k: int, m: int, n: int, type="dominant allele"):
    """
    This function assumes a biallelic mendelian locus and returns the probability that offspring of the defined 
    type will be produced by the mating of any two random individuals in the population.
    
    Possible genotypes are:
    k: number of homozygous dominant individuals (AA)
    m: number of heterozygous individuals (Aa)
    n: number of homozygous recessive individuals (aa)

    Types of calculations available: 
    dominant allele: Returns the probability that the offspring will possess at least 1 dominant allele.
    recessive allele: Returns the probability that the offspring will possess at least 1 recessive allele.
    """
    allowed_types = ["dominant allele", "recessive allele"]
    total = k + m + n
    total_pairs = total * (total -1)

    prob_aa_aa = n * (n-1)
    prob_aa_Aa = 2*n*m
    prob_Aa_Aa = m * (m-1)

    recessive_prob = (prob_aa_aa*1 + prob_Aa_Aa*0.25 + prob_aa_Aa*0.5)/total_pairs
    dominant_prob = 1 - recessive_prob
    if type == "dominant allele":
        print(dominant_prob)
        return dominant_prob
    elif type == "recessive allele":
        print(recessive_prob)
        return recessive_prob

offspring_probability(2, 2, 2, type="dominant allele")