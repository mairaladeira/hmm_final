import os
import sys
from data_structures import *
from algorithms.baum_welch_algorithm import *

__author__ = 'Maira'


def manual_input():
    print('--------------------------------------------------------------')
    print('--------------------------------------------------------------')
    print('Welcome to the HMM implementation')
    print('This project implements the following algorithms:')
    print('- Forward Backward')
    print('- Viterbi')
    print('- Baum Welch')
    print('')
    print('The required inputs should be  provided in three different file format:')
    print('- <file_name>.trans: containing the initialization and transition matrix')
    print('- <file_name>.emit: containing the emission probabilities matrix')
    print('- <file_name>.input: containing the observations')
    print('For more information about this files, please consult the readme file')
    print('')
    trans_file = input('Transition matrix file: ')
    while '.trans' not in trans_file:
        print('The transition matrix file should have the format .trans')
        trans_file = input('Transition matrix file: ')
    while not os.path.isfile(trans_file):
        print('File not found!')
        trans_file = input('Transition matrix file: ')
    emission_file = input('Emission matrix file: ')
    while '.emit' not in emission_file:
        print('The emission matrix file should have the format .emit')
        emission_file = input('Emission matrix file: ')
    while not os.path.isfile(emission_file):
        print('File not found!')
        emission_file = input('Emission matrix file: ')
    input_file = input('Observations file: ')
    while '.input' not in input_file:
        print('The observations file should have the format .input')
        input_file = input('Observations file: ')
    while not os.path.isfile(input_file):
        print('File not found!')
        input_file = input('Observations file: ')
    max_iterations = input('Max number of iterations:')
    print('--------------------------------------------------------------')
    print('--------------------------------------------------------------')
    return trans_file, emission_file, input_file, max_iterations


def automatic_input(args):
    if args[1] != '-t' or args[3] != '-e' or args[5] != '-i':
        print('Wrong arguments instantiation.')
        print('Correct use: python main.py -t <trans_file>.trans -e <emission_file>.emit -i <input_file>.input -bw_i <max_iterations>')
        sys.exit()
    trans_file = args[2]
    emit_file = args[4]
    input_file = args[6]
    if args[7] == '-bw_i':
        max_iterations = int(args[8])
    if '.trans' not in trans_file:
        print('Wrong file extension for transition matrix file.')
        print('Correct extension: .trans')
        sys.exit()
    if '.emit' not in emit_file:
        print('Wrong file extension for emission matrix file.')
        print('Correct extension: .emit')
        sys.exit()
    if '.input' not in input_file:
        print('Wrong file extension for observations file.')
        print('Correct extension: .input')
        sys.exit()
    if not os.path.isfile(trans_file):
        print('Transition matrix file not found. Please check the path')
        sys.exit()
    if not os.path.isfile(emit_file):
        print('Emission matrix file not found. Please check the path')
        sys.exit()
    if not os.path.isfile(input_file):
        print('Observation file not found. Please check the path')
        sys.exit()
    return trans_file, emit_file, input_file, max_iterations


if __name__ == "__main__":
    if len(sys.argv) == 1:
        t_file, e_file, i_file, max_it = manual_input()
    else:
        args = sys.argv
        t_file, e_file, i_file, max_it = automatic_input(args)
    trans_matrix = TransitionMatrix(t_file)
    emission_matrix = EmissionProbability(e_file)
    obs_obj = Observations(i_file)
    states = trans_matrix.states
    for obs in obs_obj.all:
        bw = BaumWelchAlgorithm(
                trans_matrix,
                emission_matrix,
                trans_matrix.states,
                trans_matrix.init_p,
                trans_matrix.finals,
                obs,
                max_it)
        bw.print_fw()
        bw.print_viterbi()
        bw.baum_welch()
        bw.print_bw()
    print('------------------------------------------')
