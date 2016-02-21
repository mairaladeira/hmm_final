from algorithms.forward_backward_algorithm import *
__author__ = 'Maira'


class ViterbiAlgorithm(ForwardBackwardAlgorithm):
    path = None
    prob_path = None

    def __init__(self, a, b, states, pi, finals, obs):
        ForwardBackwardAlgorithm.__int__(self, a, b, states, pi, finals, obs)

    def viterbi(self):
        viterbi = []
        i = 0
        prev_viterbi = {}
        self.path = []
        most_likely = ('', 0)
        for obs in self.obs:
            cur_viterbi = {}
            total_prob = 0.0
            for s in self.states:
                if i == 0:
                    cur_viterbi[s] = self.pi[s] * \
                                     self.get_emission_prob(s, obs)
                    total_prob += cur_viterbi[s]
                else:
                    max_prob = 0
                    for s2 in self.states:
                        new_prob = prev_viterbi[s2]*self.get_transition_prob(s2, s)
                        if new_prob > max_prob:
                            max_prob = new_prob
                    cur_viterbi[s] = max_prob * self.get_emission_prob(s, obs)
                    total_prob += cur_viterbi[s]
            prev_viterbi = cur_viterbi
            viterbi.append(cur_viterbi)
            most_likely = ('', 0)
            for s in cur_viterbi:
                if cur_viterbi[s] > most_likely[1]:
                    most_likely = (s, cur_viterbi[s])

            self.path.append(str(obs)+" "+str(most_likely[0]))
            i += 1
        self.prob_path = most_likely[1]/self.prob_obs

    def __getitem__(self, item):
        if item == 'prob_path':
            return self.prob_path
        if item == 'path':
            return self.path
        return str(item)+' not found in the class'

    def __str__(self):
        string = '\tP(path):'+str(self.prob_path)+'\n'
        string += '\tMost likely path:\n'
        for elem in self.path:
            string += '\t'+str(elem) + '\n'
        return string

    def print_viterbi(self):
        print('\tResults from Viterbi algorithm:')
        print('\tP(path):'+str(self.prob_path))
        print('\tMost likely path:')
        for elem in self.path:
            print('\t'+str(elem))
        print('\t--------------')

