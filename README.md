# Hidden Markov Model
DMKM project that aims to implement the algorithms to solve Rabiner's [1] three classical problems:
* Given an observation sequence O1, O2, ..., OT and the model λ = (A, B, π),
how to compute P(O|λ), the probability of the observation sequence? Solved with the Forward-Backward algorithm
* Given the observation sequence O = O1, O2, ..., OT , how to choose a state
sequence I = i1, i2, ..., i3 that maximizes the probability of this observation? Solved with the Viterbi algorithm
* How to adjust the model parameters λ = (A, B, π) to maximize P r(O|λ)? Solved with the Baum-Welch algorithm

**Usage:**
python main.py -t <trans_file>.trans -e <emission_file>.emit -i <input_file>.input -bw_i <max_iterations>

the files folder contains some examples of files that can be used for the execution of this program.

**Group:** 
Maira Ladeira and Gabriela Hernandez

[1] Rabiner, Lawrence R. "A tutorial on hidden Markov models and selected applications in
speech recognition." Proceedings of the IEEE 77.2 (1989): 257-286.

