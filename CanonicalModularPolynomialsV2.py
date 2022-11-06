import PolynomialsV2
from gmpy2 import *
import math

#given p, l, loads the lth Canonical modular polynomial from stored text file and reduces modulo p. 
#the text files contain data as a triple a,b,c which represents the monomial in x and y:  c*x^a y^b
#various methods below evaluate the polynomial at partial derivatives.

class CanonicalModularPolynomial:
    def __init__(self, p, l):
        phi = open("CanonicalModularPolynomials/Phi" + str(l) + ".txt", "r")
        self.coefficients = []
        self.l = l
        self.P = PolynomialsV2.polynomials_over(p)
        self.p = p
        for line in phi:
            data = line.split()
            pow_x = int(data[0])
            pow_y = int(data[1])
            coefficient = mpz(data[2]) % self.p
            self.coefficients.append([pow_x, pow_y, coefficient])
        phi.close()

    def evaluate(self, x, y):
        z = 0
        for monomial in self.coeffs:
            deg_x, deg_y, c = monomial
            z += (pow(x, deg_x, self.p) * (pow(y, deg_y, self.p)) * c)
        return z % self.p

    def eval_at_y(self, y):
        coefficients = [0] * (self.l + 2)
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            coefficients[deg_x] += (c * pow(y, deg_y, self.p)) % self.p
        return self.P(coefficients)

    def eval_at_x(self, x):
        coefficients = [0] * (self.l + 2)
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            coeffslist[deg_y] += c * (pow(x, deg_x, self.p)) % self.p
        return self.P(coefficients)

    def partial_x(self, x, y):
        val = 0
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            if deg_x > 0:
                val = val + c * deg_x * (pow(x, (deg_x - 1), self.p) * pow(y, deg_y, self.p)) % self.p
        return val % self.p

    def partial_y(self, x, y):
        val = 0
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            if deg_y > 0:
                val = val + (c * deg_y) * pow(y, (deg_y - 1), self.p) * pow(x, deg_x, self.p) % self.p
        return val % self.p

    def partial_xx(self, x, y):
        val = 0
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            if deg_x > 1:
                val = val + (c * deg_x * (deg_x - 1) * pow(x, (deg_x - 2), self.p) * pow(y, deg_y, self.p)) % self.p
        return val % self.p

    def partial_xy(self, x, y):
        val = 0
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            if deg_x > 0 and deg_y > 0:
                val = val + (c * deg_x * deg_y * pow(y, (deg_y - 1), self.p) * pow(x, (deg_x - 1), self.p)) % self.p
        return val % self.p

    def partial_yy(self, x, y):
        val = 0
        for monomial in self.coefficients:
            deg_x, deg_y, c = monomial
            if deg_y > 1:
                val = val + (c * deg_y * (deg_y - 1) * pow(y, deg_y - 2, self.p) * pow(x, deg_x, self.p)) % self.p
        return val % self.p
