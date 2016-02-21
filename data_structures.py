__author__ = 'Gabriela and Maira'

import copy


class TransitionMatrix(object):

    init_p = None
    trans_mat = None
    final = None
    _states = None

    def __init__(self, file):
        f = open(file, 'r')
        i = 0
        first_state = ""
        self.trans_mat = {}
        self.final = {}
        self._states = []
        for line in f:
            line = line.replace("\n", "")
            if i == 0:
                first_state = line
            else:
                data = line.split("\t")
                if len(data) != 3:
                    continue
                if data[0] not in self.trans_mat:
                    self.trans_mat[data[0]] = {}
                self.trans_mat[data[0]][data[1]] = float(data[2])
                if data[1] == 'FINAL':
                    self.final[data[0]] = float(data[2])
            i += 1
        self.init_p = copy.deepcopy(self.trans_mat[first_state])
        del self.trans_mat[first_state]
        self._states = list(self.trans_mat.keys())
        for s in self._states:
            if s not in self.init_p:
                self.init_p[s] = 0

    def __getitem__(self, item):
        if item in self.trans_mat:
            return self.trans_mat[item]
        return {}

    @property
    def initial_probabilities(self):
        return self.init_p

    @property
    def states(self):
        return self._states

    @property
    def finals(self):
        return self.final

    def __str__(self):
        return str(self.trans_mat)


class EmissionProbability:
    def __init__(self, file):
        f = open(file, 'r')
        self.emission_probability = {}
        for line in f:
            line = line.replace("\n", "")
            data = line.split("\t")
            if len(data) != 3:
                continue
            if data[0] not in self.emission_probability:
                self.emission_probability[data[0]] = {}
            self.emission_probability[data[0]][data[1]] = float(data[2])

    def __getitem__(self, index):
        if index in self.emission_probability:
            return self.emission_probability[index]
        return {}

    def __str__(self):
        return str(self.emission_probability)


class Observations:
    def __init__(self, file):
        f = open(file, 'r')
        self.paths = []
        for line in f:
            line = line.replace("\n", "")
            line = line.strip()
            data = line.split(" ")
            if len(data) > 1:
                self.paths.append(data)

    @property
    def all(self):
        return self.paths

    def __getitem__(self, path):
        return self.paths[path]
