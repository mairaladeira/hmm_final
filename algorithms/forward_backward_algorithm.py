from algorithms.algorithm import *
__author__ = 'Gabriela and Maira'


class ForwardBackwardAlgorithm(Algorithm):
    alphas = None
    betas = None
    constants = None
    prob_obs = None

    def __init__(self, a, b, states, pi, finals, obs):
        Algorithm.__int__(self, a, b, states, pi, finals, obs)

    def forward(self):
        alphas = []
        constants = []
        previous_alpha_estimated = {}
        previous_alpha = {}
        constants.append(0)
        for s in self.states:
            previous_alpha[s] = self.pi[s] * self.get_emission_prob(s, self.obs[0])
            constants[0] += previous_alpha[s]
        for s in self.states:
            previous_alpha_estimated[s] = previous_alpha[s]/constants[0]
        alphas.append(previous_alpha_estimated)

        for index in range(1, len(self.obs)):
            current_alpha = {}
            current_alpha_estimated = {}
            constants.append(0)
            for j in self.states:
                if j not in current_alpha:
                    current_alpha[j] = 0
                for i in self.states:
                    if i in previous_alpha_estimated:
                        current_alpha[j] += previous_alpha_estimated[i]*self.get_transition_prob(i, j)
                current_alpha[j] *= self.get_emission_prob(j,self.obs[index])
                constants[index] += current_alpha[j]
            for s in self.states:
                if constants[index] != 0:
                    current_alpha_estimated[s] = current_alpha[s]/constants[index]
                else:
                    current_alpha_estimated[s] = 0
            previous_alpha_estimated = current_alpha_estimated.copy()
            alphas.append(current_alpha_estimated)
        f_const = 1
        for c in constants:
            f_const *= 1/c
        prob = 1/f_const
        self.alphas = alphas
        self.constants = constants
        self.prob_obs = prob

    def backward(self):
        betas = []
        index = 1
        prev_beta = {}
        const = self.constants[::-1]
        for s in self.states:
            prev_beta[s] = 1/const[0]
        betas.append(prev_beta)
        for obs in reversed(self.obs):
            cur_beta = {}
            if index < len(self.obs):
                for s in self.states:
                    cur_beta[s] = 0.0
                    for s2 in self.states:
                        val = prev_beta[s2] * \
                                        self.get_transition_prob(s, s2) * \
                                        self.get_emission_prob(s2, obs)
                        cur_beta[s] += val
                for s in cur_beta:
                    cur_beta[s] /= const[index]
                prev_beta = cur_beta.copy()
                betas.append(cur_beta)
                index += 1
        betas = betas[::-1]
        self.betas = betas

    def __getitem__(self, item):
        if item == 'alphas':
            return self.alphas
        if item == 'betas':
            return self.betas
        if item == 'constants':
            return self.constants
        if item == 'prob_obs':
            return self.prob_obs

    def __str__(self):
        string = '\tObservation:'+str(self.obs)+'\t\n'
        string += '\tP(Observation):'+str(self.prob_obs)
        return string

    def print_fw(self):
        print('------------------------------------------')
        print('\tResults from Forward-Backward algorithm:')
        print('\tAlphas:'+str(self.alphas)+'\t')
        print('\tBetas:'+str(self.betas)+'\t')
        print('\tObservation:'+str(self.obs)+'\t')
        print('\tP(Observation):'+str(self.prob_obs))
        print('\t--------------')