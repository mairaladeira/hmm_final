from algorithms.viterbi_algorithm import *
__author__ = 'Maira'


class BaumWelchAlgorithm(ViterbiAlgorithm):
    init_Gamma = None
    sum_Gammas = None
    sum_Gammas_end = None
    sum_Xis = None
    sum_observed_Gammas = None
    sum_observed_Gammas_end = None
    iterations = None

    def __init__(self, a, b, states, pi, finals, obs, iterations):
        ViterbiAlgorithm.__int__(self, a, b, states, pi, finals, obs)
        ForwardBackwardAlgorithm.forward(self)
        ForwardBackwardAlgorithm.backward(self)
        ViterbiAlgorithm.viterbi(self)
        self.iterations = iterations

    def get_possible_obs_list(self):
        possible_obs = []
        for s in self.states:
            for o in self.B[s]:
                if o not in possible_obs:
                    possible_obs.append(o)
        return possible_obs

    def set_vars(self):
        self.sum_Gammas = {}
        self.sum_Xis = {}
        self.sum_observed_Gammas = {}
        for i in range(len(self.obs)-1):
            cur_gamma = {}
            obs = self.obs[i]
            for s in self.states:
                const = self.constants[i]
                #if const == 0:
                 #   const = 1
                gamma = (self.alphas[i][s]*self.betas[i][s])/const
                cur_gamma[s] = gamma
                if s not in self.sum_Gammas:
                    self.sum_Gammas[s] = 0
                    self.sum_observed_Gammas[s] = {}
                    self.sum_Xis[s] = {}
                self.sum_Gammas[s] += gamma
                if obs not in self.sum_observed_Gammas[s]:
                    self.sum_observed_Gammas[s][obs] = 0
                self.sum_observed_Gammas[s][obs] += gamma
                for ns in self.states:
                    xi = self.alphas[i][s] * \
                         self.get_transition_prob(s, ns) * \
                         self.get_emission_prob(ns, self.obs[i+1]) * \
                         self.betas[i+1][ns]
                    if ns not in self.sum_Xis[s]:
                        self.sum_Xis[s][ns] = 0
                    self.sum_Xis[s][ns] += xi
            if i == 0:
                self.init_Gamma = cur_gamma
            self.sum_Gammas_end = self.sum_Gammas.copy()
            self.sum_observed_Gammas_end = self.sum_observed_Gammas.copy()
            last_elem = len(self.obs)-1
            obs = self.obs[last_elem]
            for s in self.states:
                const = self.constants[last_elem]
                #if const == 0:
                 #   const = 1
                gamma = (self.alphas[last_elem][s]*self.betas[last_elem][s])/const
                if s not in self.sum_Gammas_end:
                    self.sum_Gammas_end[s] = 0
                    self.sum_observed_Gammas_end[s] = {}
                self.sum_Gammas_end[s] += gamma
                if obs not in self.sum_observed_Gammas_end[s]:
                    self.sum_observed_Gammas_end[s][obs] = 0
                self.sum_observed_Gammas_end[s][obs] += gamma

    def update_params(self, possible_obs):
        n_pi = {}
        n_a = {}
        n_b = {}
        for s in self.states:
            # All the comparisons with 0 are to avoid sparse initializations
            if s in self.init_Gamma and self.init_Gamma != 0:
                n_pi[s] = self.init_Gamma[s]
            if s not in n_a:
                n_a[s] = {}
                n_b[s] = {}
            for s2 in self.states:
                if self.sum_Xis[s][s2] != 0:
                    n_a[s][s2] = self.sum_Xis[s][s2]/self.sum_Gammas[s]
            for obs in possible_obs:
                if obs in self.sum_observed_Gammas[s] and self.sum_observed_Gammas[s][obs] != 0:
                    n_b[s][obs] = self.sum_observed_Gammas_end[s][obs]/self.sum_Gammas_end[s]
        for s in self.states:
            if n_a[s] == {}:
                del n_a[s]
            if n_b[s] == {}:
                del n_b[s]
        return n_pi, n_a, n_b

    def baum_welch(self):
        cur_a = self.A
        cur_b = self.B
        cur_pi = self.pi
        i = 0
        possible_obs = self.get_possible_obs_list()
        while (cur_a != self.A or cur_b != self.B or cur_pi != self.pi) and i < self.iterations:
            if i != 0:
                self.A = cur_a
                self.B = cur_b
                self.pi = cur_pi
                ForwardBackwardAlgorithm.forward(self)
                ForwardBackwardAlgorithm.backward(self)
            self.set_vars()
            cur_pi, cur_a, cur_b = self.update_params(possible_obs)
        self.A = cur_a
        self.B = cur_b
        self.pi = cur_pi

    def __str__(self):
        string = '\tNew initialization probabilities:\n'
        string += '\t'+str(self.pi)+'\n'
        string += '\tNew transition probabilities:\n'
        string += '\t'+str(self.A)+'\n'
        string += '\tNew emission probabilities:\n'
        string += '\t'+str(self.B)+'\n'
        return string

    def print_bw(self):
        print('\tResults from Baum-Welch algorithm:')
        print('\t'+str(self.pi))
        print('\tNew initialization probabilities:')
        print('\tNew transition probabilities:')
        print('\t'+str(self.A))
        print('\tNew emission probabilities:')
        print('\t'+str(self.B))
        print('------------------------------------------')