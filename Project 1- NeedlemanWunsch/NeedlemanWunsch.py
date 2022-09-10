import sys
import csv
import os
import numpy


# input.csv Reader
# used to store the data from the input file
# @name- load_input_file()
# @return- input_file - returns the file that stores the sequences, as strings

def load_input_file():
    try:
        file = open(sys.argv[1], 'r')
        input_file = csv.reader(file)
    except:
        print("invalid file")
        exit()
    return input_file


# Needleman-Wunsch methpd
# @name- needleman_wunsch()
# @param- seq_1 - string representation of the 1st sequence
# @param- seq_2 - string representation of the 2nd sequence
def needleman_wunsch(seq_1: str, seq_2: str):
    # Create Matrices
    matrix = numpy.zeros((len(seq_1) + 1, len(seq_2) + 1))
    checker_matrix = numpy.zeros((len(seq_1), len(seq_2)))

    # Providing the fixed scores for match, mismatch and gap penalty
    match_reward = 1
    mismatch_penalty = -1
    gap_penalty = -2

    # Fill the match checker matrix accrording to if there was a match or a mismatch
    for i in range(len(seq_1)):
        for j in range(len(seq_2)):
            if seq_1[i] == seq_2[j]:
                checker_matrix[i][j] = match_reward
            else:
                checker_matrix[i][j] = mismatch_penalty

    # Filling up the matrix using Needleman_Wunsch algorithm
    # First the initialization of some matrix values
    for i in range(len(seq_1) + 1):
        matrix[i][0] = i * gap_penalty
    for j in range(len(seq_2) + 1):
        matrix[0][j] = j * gap_penalty

# Filling the matrix with dynamic programming techniques
    for i in range(1, len(seq_1) + 1):
        for j in range(1, len(seq_2) + 1):
            matrix[i][j] = max(matrix[i - 1][j - 1] + checker_matrix[i - 1][j - 1],
                                matrix[i - 1][j] + gap_penalty,
                                matrix[i][j - 1] + gap_penalty)
    score = matrix[-1][-1]

    # Backtracing to write the sequences
    # Define our
    alignment_text_1 = ""
    alignment_text_2 = ""

    len_seq_1 = len(seq_1)
    len_seq_2 = len(seq_2)

    while len_seq_1 > 0 and len_seq_2 > 0:

        if len_seq_1 > 0 and len_seq_2 > 0 and matrix[len_seq_1][len_seq_2] == matrix[len_seq_1 - 1][len_seq_2 - 1] + \
                checker_matrix[len_seq_1 - 1][len_seq_2 - 1]:

            alignment_text_1 = seq_1[len_seq_1 - 1] + alignment_text_1
            alignment_text_2 = seq_2[len_seq_2 - 1] + alignment_text_2

            len_seq_1 = len_seq_1 - 1
            len_seq_2 = len_seq_2 - 1

        elif len_seq_1 > 0 and matrix[len_seq_1][len_seq_2] == matrix[len_seq_1 - 1][len_seq_2] + gap_penalty:
            alignment_text_1 = seq_1[len_seq_1 - 1] + alignment_text_1
            alignment_text_2 = "-" + alignment_text_2

            len_seq_1 = len_seq_1 - 1
        else:
            alignment_text_1 = "-" + alignment_text_1
            alignment_text_2 = seq_2[len_seq_2 - 1] + alignment_text_2

            len_seq_2 = len_seq_2 - 1
    return ["".join(i for i in [alignment_text_1, " ", alignment_text_2]), int(score)]


# Couldn't figure out how to extract the info and write results that well
# with numpy, so more robust way
load = load_input_file()

with open('./results.csv', 'w', newline='') as output_csv:
    result = csv.writer(output_csv, delimiter=",")
    r = next(load) + ['alignment text', 'alignment score']
    result.writerow(r)

    for row in load:
        row += needleman_wunsch(row[0], row[1])
        result.writerow(r)

    output_csv.close()