from algorithms.viterbi_algorithm import *
__author__ = 'Gabriela and Maira'


class BaumWelchAlgorithm(ViterbiAlgorithm):
    gammas = None
    xis = None
    iterations = None

    def __init__(self, a, b, states, pi, finals, obs, iterations):
        ViterbiAlgorithm.__init__(self, a, b, states, pi, finals, obs)
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
        self.gammas = []
        self.xis = []
        for i in range(len(self.obs)):
            cur_gamma = {}
            cur_xi = {}
            const = self.constants[i]
            for s in self.states:
                if s not in cur_xi:
                    cur_xi[s] = {}
                cur_gamma[s] = (self.alphas[i][s]*self.betas[i][s])*const
                for ns in self.states:
                    if i < len(self.obs)-1:
                        xi = self.alphas[i][s] * \
                             self.get_transition_prob(s, ns) * \
                             self.get_emission_prob(ns, self.obs[i+1]) * \
                             self.betas[i+1][ns]
                        cur_xi[s][ns] = xi
            self.gammas.append(cur_gamma)
            self.xis.append(cur_xi)

    def normalize_params(self, n_a, n_b, n_pi=None):
        const_pi = 0
        for s in self.states:
            const_a = 0
            const_b = 0
            if n_pi is not None and s in n_pi:
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
        if n_pi is not None:
            for s in self.states:
                if s in n_pi:
                    n_pi[s] /= const_pi
        return n_a, n_b, n_pi

    def update_params(self, possible_obs):
        n_pi = self.gammas[0]
        n_a = {}
        n_b = {}
        T = len(self.obs)
        for i in self.states:
            for j in self.states:
                n_a_top = 0
                n_a_bottom = 0
                for t in range(T-1):
                    n_a_top += self.xis[t][i][j]
                    n_a_bottom += self.gammas[t][i]
                if n_a_top != 0 and n_a_bottom != 0:
                    if i not in n_a:
                        n_a[i] = {}
                    n_a[i][j] = n_a_top/n_a_bottom
            for k in possible_obs:
                n_b_top = 0
                n_b_bottom = 0
                for t in range(T):
                    if k == self.obs[t]:
                            n_b_top += self.gammas[t][i]
                    n_b_bottom += self.gammas[t][i]
                if n_b_top != 0 and n_b_bottom != 0:
                    if i not in n_b:
                        n_b[i] = {}
                    n_b[i][k] = n_b_top/n_b_bottom
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
        return self.alphas, self.betas, self.constants
        #return self.sum_Gammas, self.sum_Xis, self.sum_observed_Gammas, self.sum_Gammas_end

    def multiple_observations(self, list_alphas, list_betas, constants, observations):
        n_a = {}
        n_b = {}
        for i in self.states:
            if i not in n_a:
                n_a[i] = {}
                n_b[i] = {}
            for j in self.states:
                n_a_top = 0
                n_a_bottom = 0
                for k in range(len(list_alphas)):
                    for t in range(len(list_alphas[k])-1):
                        n_a_top += list_alphas[k][t][i]\
                                   * self.get_transition_prob(i, j)\
                                   * list_betas[k][t+1][j]\
                                   * self.get_emission_prob(j, observations[k][t+1])
                        n_a_bottom += list_alphas[k][t][i]*list_betas[k][t][i]/constants[k][t]
                if n_a_top != 0 and n_a_bottom != 0:
                    n_a[i][j] = n_a_top/n_a_bottom
            possible_obs = self.get_possible_obs_list()
            for obs in possible_obs:
                n_b_top = 0
                n_b_bottom = 0
                for k in range(len(list_alphas)):
                    for t in range(len(list_alphas[k])-1):
                        if observations[k][t] == obs:
                            n_b_top += list_alphas[k][t][i]*list_betas[k][t][i]/constants[k][t]
                        n_b_bottom += list_alphas[k][t][i]*list_betas[k][t][i]/constants[k][t]
                if n_b_top != 0 and n_b_bottom != 0:
                    n_b[i][obs] = n_b_top/n_b_bottom

            n_a, n_b, n_pi = self.normalize_params(n_a, n_b)
        return n_a, n_b, self.pi

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
        for a in self.A:
            print('\t\t'+str(a)+': '+str(self.A[a]))
        print('\tNew emission probabilities:')
        for b in self.B:
            print('\t\t'+str(b)+': '+str(self.B[b]))
        print('------------------------------------------')