from algorithms.forward_backward_algorithm import *
__author__ = 'Gabriela and Maira'


class ViterbiAlgorithm(ForwardBackwardAlgorithm):
    path = None
    prob_path = None

    def __init__(self, a, b, states, pi, finals, obs):
        ForwardBackwardAlgorithm.__int__(self, a, b, states, pi, finals, obs)

    def viterbi(self):
        deltas = []
        psis = []
        self.path = []
        for t in range(len(self.obs)):
            delta = {}
            psi = {}
            for i in self.states:
                if t == 0:
                    delta[i] = self.pi[i]*self.get_emission_prob(i, self.obs[t])
                    psi[i] = 0
                else:
                    max_delta = 0
                    max_psi = None
                    for j in self.states:
                        p = deltas[t-1][j]*self.get_transition_prob(j, i)*self.get_emission_prob(i, self.obs[t])
                        if p > max_delta:
                            max_delta = p
                            max_psi = j
                    delta[i] = max_delta
                    psi[i] = max_psi
            deltas.append(delta)
            psis.append(psi)
        max_prob = 0
        max_state = None
        for s in deltas[len(deltas)-1]:
            if deltas[len(deltas)-1][s] > max_prob:
                max_prob = deltas[len(deltas)-1][s]
                max_state = s
        for t in range(len(self.obs)-1):
            mp = 0
            ms = None
            for s in deltas[t+1]:
                if deltas[t+1][s] > mp:
                    mp = deltas[t+1][s]
                    ms = s
            self.path.append(str(self.obs[t])+" "+str(psis[t+1][ms]))
        self.path.append(str(self.obs[len(self.obs)-1])+" "+str(max_state))
        self.prob_path = max_prob/self.prob_obs

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

