
import PolynomialsV2
import json
from gmpy2 import *

#loads the lth classical modular polynomial modulo p from text file. Data in text file is stored as
#(a,b) c representing the monomial c * x^a y^b
#class has various methods to evaluate at specific values and partial derivatives
class ModularPolynomial:
    def __init__(self, p, l):
        phi = open('Modular Polynomials/ClassicalModularPolynomials/phi_j_' + str(l) + '.txt', 'r')
        self.coefficients = []
        self.l = l
        self.P = PolynomialsV2.polynomials_over(p)
        self.p = p
        for line in phi:
            data = line.split(' ')
            powers = json.loads(data[0])
            coefficient = mpz(data[1].split('\\')[0]) % self.p
            powers = [int(powers[0]), int(powers[1])]
            self.coefficients.append([powers, coefficient])
        phi.close()

    def sub(self, x, y):
        z = 0
        for term in self.coefficients:
            z += pow(x, term[0][0], self.p) * pow(y, term[0][1], self.p) * term[1] % self.p
            if term[0][1] > term[0][0]:
                z += pow(y, term[0][0], self.p) * pow(x, term[0][1], self.p) * term[1] % self.p
        return z % self.p

    def eval_at_y(self, k):
        coefficients = [0] * (self.l + 2)
        for term in self.coefficients:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            coefficients[i] += c * pow(k, j, self.p)
            if i > j:
                coefficients[j] += c * pow(k, i, self.p)
        return self.P(coefficients)

    def partial_x(self, a, b):
        res = 0
        for term in self.coefficients:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i > 0:
                res = res + i * c * (pow(a, (i - 1), self.p) * (pow(b, j, self.p))) % self.p
            if i > j > 0:
                res = res + j * c * (pow(a, (j - 1), self.p) * pow(b, i, self.p))
        return res % self.p

    def partial_y(self, a, b):
        res = 0
        for term in self.coefficients:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i > 0:
                res = res + i * c * (pow(a, j, self.p) * pow(b, (i - 1), self.p)) % self.p
            if i > j > 0:
                res = res + j * c * pow(a, i, self.p) * pow(b, (j - 1), self.p) % self.p
        return res % self.p

    def partial_xx(self, a, b):
        res = 0
        for term in self.coefficients:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i > 1:
                res = res + (i - 1) * i * c * pow(a, (i - 2), self.p) * pow(b, j, self.p) % self.p
            if i > j > 1:
                res = res + (j - 1) * j * c * pow(a, (j - 2), self.p) * pow(b, i, self.p) % self.p
        return res % self.p

    def partial_yy(self, a, b):
        res = 0
        for term in self.coefficients:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i > 1:
                res = res + i * (i - 1) * c * pow(a, j, self.p) * pow(b, (i - 2), self.p) % self.p
            if i > j > 1:
                res = res + j * (j - 1) * c * pow(a, i, self.p) * pow(b, (j - 2), self.p) % self.p
        return res % self.p

    def partial_xy(self, a, b):
        res = 0
        for term in self.coefficients:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i > 0 and j > 0:
                res = res + i * j * c * pow(a, (i - 1), self.p) * pow(b, (j - 1), self.p) % self.p
            if i > j > 0:
                res = res + j * i * c * pow(a, (j - 1), self.p) * pow(b, (i - 1), self.p) % self.p
        return res % self.p


