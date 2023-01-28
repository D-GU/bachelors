import sympy as sp
import numpy as np
import sage


class GeomZeta:
    def __init__(
            self, a=None, b=None, c=None, A=None, B=None, C=None, P=None, Q=None, R=None, S=None
    ):
        self.a = a
        self.b = b
        self.c = c
        self.A = A
        self.B = B
        self.C = C
        self.P = P
        self.Q = Q
        self.R = R
        self.S = S

    @property
    def collinear(self):
        if sp.factor(self.a * sp.conjugate(self.b) - sp.conjugate(self.b) - self.a) == 0:
            return True
        return False

    @property
    def concurrent(self):
        if sp.factor(
                sp.intersection(self.P, self.a, self.Q, self.b) - sp.intersection(self.Q, self.b, self.R, self.c)) == 0:
            return True
        return False

    @property
    def concyclic(self):
        if sp.factor(
                (self.P - self.R) * (self.Q - self.S) / (self.Q - self.R) / (self.P - self.S) - sp.conjugate(
                    (self.P - self.R) * (self.Q - self.S) / (self.Q - self.R) / (self.P - self.S))) == 0:
            return True
        return False

    @property
    def equidistant(self):
        if sp.factor(sp.sqf_norm(self.P - self.Q) - sp.sqf_norm(self.R - self.S)) == 0:
            return True
        return False
