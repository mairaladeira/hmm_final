from algorithms.viterbi_algorithm import *
__author__ = 'Gabriela and Maira'


class BaumWelchAlgorithm(ViterbiAlgorithm):
    init_Gamma = None
    sum_Gammas = None
    sum_Gammas_end = None
    sum_Xis = None
    sum_observed_Gammas = None
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
                gamma = (self.alphas[i][s]*self.betas[i][s])*const
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
            last_elem = len(self.obs)-1
            obs = self.obs[last_elem]
            for s in self.states:
                const = self.constants[last_elem]
                #if const == 0:
                 #   const = 1
                gamma = (self.alphas[last_elem][s]*self.betas[last_elem][s])*const
                if s not in self.sum_Gammas_end:
                    self.sum_Gammas_end[s] = 0
                    self.sum_observed_Gammas[s] = {}
                self.sum_Gammas_end[s] += gamma
                if obs not in self.sum_observed_Gammas[s]:
                    self.sum_observed_Gammas[s][obs] = 0
                self.sum_observed_Gammas[s][obs] += gamma

    def normalize_params(self, n_a, n_b, n_pi):
        const_pi = 0
        for s in self.states:
            const_a = 0
            const_b = 0
            if s in n_pi:
                const_pi += n_pi[s]
            if s in n_a:
                for s2 in n_a[s]:
                    const_a += n_a[s][s2]
                for s2 in n_a[s]:
                    n_a[s][s2] /= const_a
            if s in n_b:
                for obs in n_b[s]:
                    const_b += n_b[s][obs]
                for obs in n_b[s]:
                    n_b[s][obs] /= const_b
        for s in self.states:
            if s in n_pi:
                n_pi[s] /= const_pi
        return n_a, n_b, n_pi

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
                    n_b[s][obs] = self.sum_observed_Gammas[s][obs]/self.sum_Gammas_end[s]
        for s in self.states:
            if n_a[s] == {}:
                del n_a[s]
            if n_b[s] == {}:
                del n_b[s]
        n_a, n_b, n_pi = self.normalize_params(n_a, n_b, n_pi)
        return n_pi, n_a, n_b

    def baum_welch(self):
        i = 1
        possible_obs = self.get_possible_obs_list()
        self.set_vars()
        cur_pi, cur_a, cur_b = self.update_params(possible_obs)
        while (cur_a != self.A or cur_b != self.B or cur_pi != self.pi) and i < self.iterations:
            self.A = cur_a
            self.B = cur_b
            self.pi = cur_pi
            ForwardBackwardAlgorithm.forward(self)
            ForwardBackwardAlgorithm.backward(self)
            self.set_vars()
            cur_pi, cur_a, cur_b = self.update_params(possible_obs)
            i += 1
        self.A = cur_a
        self.B = cur_b
        self.pi = cur_pi
        return self.sum_Gammas, self.sum_Xis, self.sum_observed_Gammas, self.sum_Gammas_end

    def multiple_observations(self, list_sum_gammas, list_sum_xis, list_sum_gammas_observed, list_sum_gammas_end):
        n_a = {}
        n_b = {}
        for i in self.states:
            if i not in n_a:
                n_a[i] = {}
                n_b[i] = {}
            for j in self.states:
                n_a_top = 0
                n_a_bottom = 0
                for k in range(len(list_sum_gammas)):
                    for t in range(len(list_sum_gammas[k])):
                        if j in list_sum_xis[k][t][i] and j in list_sum_gammas[k][t][i]:
                            n_a_top += list_sum_xis[k][t][i][j]
                            n_a_bottom += list_sum_gammas[k][t][i][j]
                if n_a_bottom != 0:
                    n_a[i][j] = n_a_top/n_a_bottom
            possible_obs = self.get_possible_obs_list()
            for obs in possible_obs:
                n_b_top = 0
                n_b_bottom = 0
                for k in range(len(list_sum_gammas_end)):
                    for t in range(len(list_sum_gammas_end[k])):
                        if obs in list_sum_gammas_observed[k][t][i]:
                            n_b_top += list_sum_gammas_observed[k][t][i][obs]
                        if obs in list_sum_gammas_end:
                            n_b_bottom += list_sum_gammas_end[k][t][i][obs]
                if n_b_bottom != 0:
                    n_b[i][obs] = n_b_top/n_b_bottom
        return n_a, n_b

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
        print('\tNew initialization probabilities:')
        print('\t'+str(self.pi))
        print('\tNew transition probabilities:')
        print('\t'+str(self.A))
        print('\tNew emission probabilities:')
        print('\t'+str(self.B))
        print('------------------------------------------')