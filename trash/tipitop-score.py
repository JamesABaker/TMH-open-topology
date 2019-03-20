# This is nowhere near complete.


'''
The input sequence should be a single TMH with 5 flanking residues either side.
However, flanks are often fuzzy in biology, so determining the exact TMH
boundary at a single residue position is difficult, if not, impossible.
Users are encouraged to submit several versions of the TMH with different
boundaries.
'''


# Replace the below AA sequence with your own.
def pinot(input_seq, tmh_coordinates)


    def characters(input_sequence):
        amino_acids = ["I", "V", "L", "F", "C", "M", "A", "G", "T",
                       "S", "W", "Y", "P", "H", "E", "Q", "D", "N", "K", "R"]
        character_check = True
        for i in list(input_sequence):
            if str(i) not in amino_acids:
                print("Character detected that do not represent an amino acid. Please remove any non-AA characters or convert the characters to upper case.")
                character_check = False
        if character_check == True:
            return(True)
        else:
            return(False)


    def length(input_sequence):
        lower_length_cutoff = 20
        higher_length_cutoff = 40
        if len(input_sequence) > lower_length_cutoff and len(input_sequence) < higher_length_cutoff:
            return(True)
        else:
            print("Sequence must be between", lower_length_cutoff,
                  "and", higher_length_cutoff, "residues long.")
            return(False)


    def topologyscore(input_sequence):
        print("\nInside to outside score")
        for position, residue in enumerate(list(input_sequence)):
            print(residue, "at position", position, "scores", )
        print("\nOutside to inside score")
        for position, residue in enumerate(list(input_sequence[::-1])):
            print(residue, "at position", position, "scores", )
        return("Topology-score, Likelihood score")


    # The sequence validity is set to false, then the checks are run to see if it can be set to true.
    sequence_integrity = False

    if characters(str(input_sequence)) == True and length(str(input_sequence)) == True:
        print("Sequence valid")
        sequence_integrity = True


    if sequence_integrity == True:
        result = topologyscore(input_sequence)

    return result
