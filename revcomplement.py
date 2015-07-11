def revcomplement(s):
    """reverse complement a nucleotide string"""
    complement = {"A": "T", "T":"A", "G":"C", "C":"G"}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

if __name__=="__main__":
   print(revcomplement("TAACCCTAACCCTAACCCT"))
