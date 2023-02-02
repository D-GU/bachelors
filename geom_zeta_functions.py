import sympy as sp
import complex as cp
import numpy as np
import sage


class GeomZeta:
    def __init__(
            self, a=None, b=None, c=None, A=None, B=None, C=None, P=None, Q=None, R=None, S=None, Z=0
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
        self.Z = Z

    @property
    def equidistant_construct(self):
        eq = self.projection
        if sp.factor(sp.sqf_norm(self.P - self.Q) - sp.sqf_norm(self.R - eq)) == 0:
            return True
        return False

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

    @property
    def perpendicular(self):
        if sp.factor(self.a * sp.conjugate(self.b) + sp.conjugate(self.a) * self.b) == 0:
            return True
        return False

    @property
    def similar(self):
        if sp.factor(
                self.A * sp.conjugate(self.Q) +
                self.B * sp.conjugate(self.R) +
                self.C * sp.conjugate(self.P) -
                self.B * sp.conjugate(self.P) -
                self.C * sp.conjugate(self.Q) -
                self.A * sp.conjugate(self.R)
        ) == 0:
            return True

        return False

    @property
    def radical(self):
        return sp.factor(
            (self.Q * sp.conjugate(self.Q) -
             self.Q * sp.conjugate(self.P) -
             sp.conjugate(self.Q) * self.P -
             self.S * sp.conjugate(self.S) +
             self.S * sp.conjugate(self.R) +
             sp.conjugate(self.S) * self.R -
             sp.conjugate(self.P) * self.R +
             self.P * sp.conjugate(self.R)) / 2 / sp.conjugate(self.R - self.P))

    @property
    def area(self):
        return sp.factor(
            ((sp.conjugate(self.B) - sp.conjugate(self.A)) *
             (self.C - self.A) - (self.B - self.A) *
             (sp.conjugate(self.C) - sp.conjugate(self.A))) / 4 / sp.I)

    @property
    def projection(self):
        return sp.intersection(self.Z, sp.I * self.a, self.P, self.a)

    @property
    def center(self):
        return sp.intersection(
            (self.P + self.Q) / 2, sp.I * (self.P - self.Q),
            (self.Q + self.R) / 2, sp.I * (self.Q - self.R))

    @property
    def centroid(self):
        return sp.factor((self.A + self.B + self.C) / 3)

    @property
    def orthocenter(self):
        return sp.intersection(self.A, sp.I * (self.B - self.C), self.B, sp.I * (self.A - self.C))

    def symmetry(self, Z, P):
        return sp.factor(2 * self.P - self.Z)

    def reflection(self, Z, P, a):
        return self.symmetry(self.Z, self.projection)

    def second_point(self):
        return self.reflection(self.P, self.Q, sp.I * self.a)

    def scalar_product(self):
        return sp.factor((self.a * sp.conjugate(self.b) + sp.conjugate(self.a) * self.b) / 2)
