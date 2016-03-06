__author__ = 'Gabriela and Maira'


class Algorithm(object):

    A = None
    B = None
    pi = None
    finals = None
    obs = None
    states = None

    def __int__(self, a, b, states, pi, finals, obs):
        """
        Define all data for the algorithms
        :param a: Transition matrix
        :param b: Emission matrix
        :param states: possible states
        :param pi: Init probabilities
        :param finals: final probabilities
        :param obs: Observations
        """
        self.A = a
        self.B = b
        self.states = states
        self.pi = pi
        self.finals = finals
        self.obs = obs

    def get_transition_prob(self, s1, s2):
        try:
            if s2 in self.A[s1]:
                return self.A[s1][s2]
        except KeyError:
            return 0
        return 0

    def get_emission_prob(self, state, obs):
        try:
            if obs in self.B[state]:
                return self.B[state][obs]
        except KeyError:
            return 0
        return 0
