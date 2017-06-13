
import random # to be able to generate random seeds for our nucleotide sequence

bases = {0: "A", 1: "C", 2: "G", 3: "T"} # a conversion table for random numbers to nucleotide bases
myrandomsequence = [] # declare an empty list


length = int(input("Length: ")) # asks for an integer
if type(length) == int and length > 0: # checks if the input is a natural number
    print("Thank you!\n\n>myrandomsequence")

    for i in list(range(0, length)): # iterates through a sequence of "length"-length
        myrandomsequence.append(bases[random.randint(0,3)]) # generates a random number
    print("".join(myrandomsequence)) # prints the list of generated nucleotides without spacing

else:
    raise ValueError("Input must be a non-negative integer.")  # returns an error message






