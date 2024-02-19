genCode = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
           'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
           'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
           'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
           'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
           'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
           'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
           'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
           'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
           'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
           'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
           'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
           'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
           'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
           'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
           'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

def translate(nucSeq, genCode = genCode): #note optional
    nucSeq = nucSeq.upper()
    aaSeq = ''
    for i in range(len(nucSeq) // 3):  # Floor division, returns the integer part of the length divided by three.
        codon = nucSeq[i*3:i*3+3]
        aaSeq = aaSeq + genCode.get(codon,'X')
    return(aaSeq)

def readFasta(seqfile):
    '''Read sequences from fasta file and create a list of tuples (name, seq)'''
    # We'll read the first line and check if it starts with '>' if it
    # does we assign this to the first name and continue. If not we return
    # an error message.
    filein = open(seqfile,'r')
    line = filein.readline()
    if line[0] != '>':
        return "File should start with a sequence in the first line.\n"
    nucseqs = []    # create an empty list
    aaseqs = []
    while line and line[0] == '>':
        name = line.lstrip(">").strip()
        # the first strip removes the '>' in beginning of sequence could
        # also be: name = line[1:]. The second strip removes newlines and
        # whitespaces in the ends of the string.
        line = filein.readline()
        # we have the name of the first sequence, let's get the sequence
        # we can read the file, line-by-line until we reach a line that
        # starts with '>' that signifies the end of the sequence
        seq = ''
        while line and line[0] != '>':
            seq = seq + line.strip()
            line = filein.readline()

        aaseq = translate(seq)
        nucseqs.append((name,seq))
        aaseqs.append((name,aaseq))
    filein.close()
    return nucseqs, aaseqs

def writeFasta(seqfile, seqs):
    fileout = open(seqfile, 'w')
    for seq in seqs:
        fileout.write(">{name}\n{sequence}\n".format(name=seq[0],sequence=seq[1]))
    fileout.close()
    return 0
