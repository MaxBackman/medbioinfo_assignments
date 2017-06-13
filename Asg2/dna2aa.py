codon2aa = {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",  # Alanine
            "TGT": "C", "TGC": "C",  # Cysteine
            "GAT": "D", "GAC": "D",  # Aspartic acid
            "GAA": "E", "GAG": "E",  # Glutamic acid
            "TTT": "F", "TTC": "F",  # Phenylalanine
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",  # Glycine
            "CAT": "H", "CAC": "H",  # Histidine
            "ATT": "I", "ATC": "I", "ATA": "I",  # Isoleucine
            "AAA": "K", "AAG": "K",  # Lysine
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "TTA": "L", "TTG": "L",  # Luecine
            "ATG": "M",  # Methionine
            "AAT": "N", "AAC": "N",  # Asparagine
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",  # Proline
            "CAA": "Q", "CAG": "Q",  # Glutamine
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",  # Arginine
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",  # Serine
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",  # Threonine
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",  # Valine
            "TGG": "W",  # Tryptophan
            "TAT": "Y", "TAC": "Y", # Tyrosine

            "TAA": "*",
            "TAG": "*",
            "TGA": "*",
            }

# Stop codons: TAA, TAG, TGA
"""" A dictionary for converting codons to single-letter amino acid representations """


def remove_comments(fasta):

    """ Takes a fasta string and removes comments on every line """

    line_array_no_comments = []
    for i in list(range(0, int(len(fasta)))):  # remove comments from lines
        line_array_no_comments.append(fasta[i].split(" ")[0])  # the last '[0]' is important to access the first element in the split

    return line_array_no_comments


def convert2aa(sequence):

    """ Takes a line of codons and converts them to single-digit aa-representations in a joined string"""

    # sequence = "".join([x.upper() for x in sequence])  # converts lowercase to uppercase

    number_of_codons = len(sequence)/3
    aa_seq = []

    for nmbr in list(range(1, int(number_of_codons)+1)):  # goes through each codon converting it to an aa

        if "".join([x.upper() for x in sequence])[nmbr*3-3:nmbr*3] in codon2aa:
            aa_seq.append(codon2aa["".join([x.upper() for x in sequence])[nmbr*3-3:nmbr*3]])
        else:
            aa_seq.append("XXX")

    return "".join(aa_seq)


def reverse_complement(sequence):

    """ Takes a line of codons and returns the reverse complement """

    reverse_complement_sequence = ""

    for i in list(range(0, len(sequence))):
        if sequence[i] == "A":
            reverse_complement_sequence += "T"
        elif sequence[i] == "T":
            reverse_complement_sequence += "A"
        elif sequence[i] == "G":
            reverse_complement_sequence += "C"
        elif sequence[i] == "C":
            reverse_complement_sequence += "G"
        else:
            reverse_complement_sequence += sequence[i]

    return reverse_complement_sequence[::-1]


def remove_empty_strings(lines):

    """ Removes empty lines from the fasta-file """

    lines_formatted = []

    for i in list(range(0, len(lines))):

        if not lines[i]:
            continue

        else:
            lines_formatted.append(lines[i])

    return lines_formatted


def join_sequences(lines):

    """ Takes a fasta string and joins all the lines with sequences into one-line sequences """

    j = 0
    lines_formatted = []

    for i in list(range(0, len(lines))):

        if lines[i][0] == ">":
            lines_formatted.append(lines[i])  # lets it be
            lines_formatted.append("")
            j = len(lines_formatted)-1

        else:
            lines_formatted[j] += lines[i]

    return lines_formatted


def longest_orf(sequence):

    """ Function that takes a sequence and finds the longest orf (also looking in the reverse complement) """

    sequence = "".join([x.upper() for x in sequence])  # makes the characters uppercase

    length_1 = 0  #store length
    length_1_pos = 0  #store position
    length_2 = 0  #-''-
    length_2_pos = 0  #-''-

    for i in list(range(0, int(len(sequence)))):  # for i in list(range(0, int(len(sequence)/3))):

        #print(i)
        length1 = distance_to_stop(sequence, i)  # starting at each codon

        if length1 > length_1:
            length_1 = length1
            length_1_pos = i
        else:
            pass

        length2 = distance_to_stop(reverse_complement(sequence), i) # starting at each codon in the reverse complement

        if length2 > length_2:
            length_2 = length2
            length_2_pos = i
        else:
            pass

    # print(convert2aa(sequence[length_1_pos:length_1*3]), convert2aa(reverse_complement(sequence)[length_2_pos:length_2*3]))

    return max(convert2aa(sequence[length_1_pos:length_1*3]), convert2aa(reverse_complement(sequence)[length_2_pos:length_2*3]), key=len)  # has '*3' because of the input from distance_to_stop()



def distance_to_stop(sequence, start):

    """ Function that takes a sequence and a starting position and finds the distance to the first stop codon
     the return value is the number of codons from the start and to right before the stop codon """

    counter = 0
    for i in list(range(start, int(len(sequence)/3))):  #

        if codon2aa[sequence[i*3:(i+1)*3]] == "*":  # does it handle ambiguities?
            break

        else:
            counter += 1
            continue

    return counter


def final(lines):

    """ Output """

    lines_formatted = []

    for i in list(range(0, len(lines))):

        if lines[i][0] == ">":
            lines_formatted.append(lines[i])
        elif "X" in convert2aa(lines[i]):  # handle ambiguities
            lines_formatted.append(convert2aa(lines[i]))

        else:
            lines_formatted.append(longest_orf(lines[i]))

    return lines_formatted


# File handling
#

my_path = str(input("Which file? "))  # ask for a file

my_file = open(my_path, "r")  # open that file

file_contents = my_file.read()  # read it in to a string

line_array = file_contents.split('\n')  # read each line into a list

line_array = remove_comments(line_array) # removes comments, only takes the first string on each line

line_array = remove_empty_strings(line_array) # removes "''" strings, i.e. empty lines

line_array = join_sequences(line_array) # joins sequences written with line-breaks onto the same line

line_array = final(line_array) # does the last step with conversion of DNA into the longest aa-sequences/ORFs

print('\n'.join(line_array))

# Closing file
#

my_file.close()

# testdata/an_exon.fa
# sequence = "ATCAATTCGATTCCAAAACTGTTTTTAACTGTATTGCTAAGCTGCCTTCCACTAGGTTGGTCAGGAAGGATAGTAGAGACAGAGGAAATGGCAGAAGGGACTGGGGGGGTGGGGACAAGATGCTGGCTCTGTTGCATTTTGAAGGCCACATCTGCCAGCCTGGGAGTGGCCTGTGGTGGCCAGGCAACCCAGCTTGATGTTGCCTCTTCTGTCCTGAGGATGGAACAGAGGCAGGTGAGAAGCTTTCTGTGGCTGCTGCAGCAAAACCAACCTATGACCAGGGGCTTTGGCTGCCATTGCCCTCCCTCCAGCCAGAAGTAAGGAGAAGGGAGCCATGCTTTTGGTGGGTAAGCTCCGCGAACAAGAACTGTGGCCCTGGGAGCTTTACTCTCAGTCCTAGCTTCTCCAGGAGAGGGTCAGCCACTGGATGCTTGGAGTCCACACTAACAGCACAACTGATTCTCGGCAGCTGACGACGGGAGTGGAGTATATAGTCACTGTGAAGATTGGCTGGACCAAATGCAAGAGGAATGACACGAGCAATTCTTCCTGCCCCCTGCAAAGCAAGAAGCTGAGAAAGGTGTGTAGTGGAAGTCATCCCAAAGGGACTTGTGAAGTGAATATTTTAAAAAGCAACAGCTAGAAGTTGGCAGTTTTGCCACTTGTGTGGATTGGGGCCAGCTGGGGCCCAGCCCTGGGGCTGTCTGGCACCTCCTGAGGCACCTGGGTGAAGGAGGGAGGCTTTATATTTGGGCTGGTACAGCCCTGCACCCAGGCCTCTGGTTGCCTGGACCCTTCTCCAACCTGGCCTGCTCCCCAGACCCTCTGCGTGAAAGCAGTGGTTCCCTTTCCTTGGGGTCCTTCTTCTGTCTTCTTTCTCCCTCTCCAGGTCTCTAAGGTTGTCCAATTTCCTCCTTGTTCCCTTCTATTTTTCCTGGTCCATGAGTGGCACCCCTTGAAAGCTGGCAGATACAAGGGGCTCCCCCATGCCTGGGGAGGTGCAGGAGGAGGTAACGATGCTATTGCTCCTTGCTGGCCTGCCTGGCCCACATGACAGGGTCCGGCCAGGCCATAGACTCCAGGTTCACCTGCTGTGGGTCCCAAGTGGCTGGGTAGGGGGAAGACAGGGTGACATTCCTGTTCCAGACTGGCACTCCTTTCTGTTTCTACCCTTGACTGTGGGCAGCTCCCCTCTCTGCCCATGGGGTTCATGAGCAGCCACAGGGAGAGACCTGGGGAACAGGTGCCTGTTGCTTGCCTTTGAAAAATACTTCCCCCAAAGTCTGTTTTCTTTCTTCTTCTACAGAGTTTAATTTGCGAGTCTTTGATATACACCATGCCCTGGATAAACTATTTCCAGCTCTGGAACAATTCCTGTCTGGAGGCCGAGCATGTGGGCAGAAACCTCAGATGAGGGCTCATATGATTGAGTTGTGCACTGGCTGTTATTAAACTGTAAAGGATCA"
