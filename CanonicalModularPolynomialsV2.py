#import FiniteFieldZP
import PolynomialsV2
from ClassicalModularPolynomialsV2 import*
from QuotientRingV2 import *
from Euclidean import*
#import math
#import scipy.special
from gmpy2 import*
import math
from FracV2 import*
from PointPolynomialArithmeticV2 import*
class CanonicalModularPolynomial():
    def __init__(self,p,l):
        Phi = open("CanonicalModularPolynomials/Phi"+str(l)+".txt", "r")
        self.coeffs = []
        self.l =l
        self.P = PolynomialsV2.PolynomialsOver(p)
        self.p = p
        for line in Phi:
            inf = line.split()
            powx = int(inf[0])
            powy = int(inf[1])
            coeff = mpz(inf[2])%self.p
            self.coeffs.append([powx, powy, coeff])
            
            
    def evaluate(self,x,y):
        z= 0
        for monomial in self.coeffs:
            degx,degy,c = monomial
            z += (pow(x,degx,self.p)*(pow(y, degy,self.p))*c)
        return z%self.p
    
    
    def eval_at_y(self,y):
        coeffslist = [0]*(self.l+2)
        for monomial in self.coeffs:
            degx,degy,c = monomial
            coeffslist[degx] += (c*pow(y, degy,self.p))%self.p
        return self.P(coeffslist)
    
    def eval_at_x(self,x):
        coeffslist = [0]*(self.l+2)
        for monomial in self.coeffs:
            degx,degy,c = monomial
            coeffslist[degy] += c*(pow(x,degx, self.p))%self.p
        return self.P(coeffslist)
            
    def partial_x(self,x,y):
        val =0
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degx>0:
                val = val+ c*degx*(pow(x, (degx-1), self.p)*pow(y, degy, self.p))%self.p
        return val%self.p
    
    def partial_y(self,x,y):
        val = 0
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degy>0:
                val = val+ (c*degy)*pow(y, (degy-1), self.p)*pow(x, degx, self.p)%self.p
        return val%self.p
    
    def partial_xx(self,x,y):
        val = 0
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degx>1:
                val = val+ (c*degx*(degx-1)*pow(x, (degx-2), self.p)*pow(y, degy, self.p))%self.p
        return val%self.p
        
        
    def partial_xy(self,x,y):
        val =0
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degx>0 and degy>0:
                val = val+ (c*degx*degy*pow(y, (degy-1), self.p)*pow(x, (degx-1), self.p))%self.p
        return val%self.p
    
    
    def partial_yy(self,x,y):
        val = 0
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degy>1:
                val = val+ (c*degy*(degy-1)*pow(y, degy-2, self.p)*pow(x, degx, self.p))%self.p
        return val%self.p

    def Eisensteins(self,a,b):
        E_4 = -48 *a
        E_6 = 864 * b
        return E_4 % self.p, E_6 % self.p

    def j_inv_prime(self,a,b,j):
        E_4, E_6 = self.Eisensteins(a,b)
        jprime = (-E_6 * pow(E_4, self.p - 2, self.p)) * j
        return jprime % self.p

    def PrimeCheck(self, l, ModularPolynomial):
        Q = GenQuotientRing(ModularPolynomial, self.P)
        x = Q(self.P.x)
        Xq = (x ** self.p).polynomial
        FieldPoly = Xq - (x.polynomial)
        GCD = ModularPolynomial.gcd(FieldPoly)
        lead = GCD.coefficients[-1]
        GCD = GCD * self.P([pow(lead, self.p - 2, self.p)])
        return GCD

    def IsogenousCurve(self, l, ModularPolynomial, gcd,a,b,j):
        # print(gcd.deg())

        q = self.p
        if gcd.deg() == 1:
            jtilde = -gcd.coefficients[0] % q
        if gcd.deg() == 2:
            Delta = gcd.coefficients[1] ** 2 - (4) * gcd.coefficients[0] * gcd.coefficients[2] % q
            disc = modular_sqrt(Delta, q) % q
            jtilde = (-gcd.coefficients[1] + disc) * pow(2 * gcd.coefficients[2], -1, q) % q

        jprime = self.j_inv_prime(a,b,j)
        Phi_x = ModularPolynomial.partial_x(j, jtilde) % q
        Phi_y = ModularPolynomial.partial_y(j, jtilde) % q
        jtilde_prime = -(jprime * Phi_x) * pow(l * Phi_y, -1, q) % q
        atilde = -(jtilde_prime ** 2) * pow(48 * jtilde * (jtilde - 1728), -1, q) % q
        btilde = -(jtilde_prime ** 3) * pow(864 * (jtilde ** 2) * (jtilde - (1728)), -1, q) % q
        return jprime, jtilde, jtilde_prime, atilde, btilde

    def ElkiesProcedure(self, l, Phi_l, gcd,a,b,j):
        q = self.p
        jprime, jtilde, jtilde_prime, atilde, btilde = self.IsogenousCurve(l, Phi_l, gcd,a,b,j)
        print("Classical jtilde,", jtilde)
        t1 = jprime ** 2 * Phi_l.partial_xx(j, jtilde) % q
        t2 = 2 * l * jprime * jtilde_prime * Phi_l.partial_xy(j, jtilde) % q
        t3 = (l ** 2) * (jtilde_prime ** 2) * Phi_l.partial_yy(j, jtilde) % q

        E4_ql = -48 * atilde % q
        E6_ql = 864 * btilde % q
        E4, E6 = self.Eisensteins(a,b)
        vii14 = -(t1 + t2 + t3) * pow(jprime * Phi_l.partial_x(j, jtilde), -1, q) % q

        T1 = pow(2, -1, q) * l * vii14 % q
        T2 = pow(4, -1, q) * l * ((E4 ** 2) * pow(E6, -1, q) - l * ((E4_ql ** 2) * pow(E6_ql, -1, q))) % q
        T3 = pow(3, -1, q) * l * (E6 * pow(E4, -1, q) - l * (E6_ql * pow(E4_ql, -1, q))) % q
        '''coefficient of x^((l-1)/2-1)'''
        p1 = (T1 + T2 + T3) % q

        atilde = (l ** 4) * atilde % q
        btilde = (l ** 6) * btilde % q
        d = int((l - 1) / 2)

        C_k, C_ktilde = self.WeierstrassCoeffecients(d, a,b, atilde, btilde)

        coeffs = self.KernelPolyCoefficients(p1, d, C_k, C_ktilde)
        return self.P(coeffs)

    def ElkiesProcedureV2(self, l, Phi, gcd,a,b,j):
        q = self.p
        E_4 = -a * pow(3, -1, q) % q
        E_6 = -b * pow(2, -1, q) % q
        Delta = (E_4 ** 3 - E_6 ** 2) * pow(1728, -1, q) % q

        if gcd.deg() == 1:
            # print('deg = 1')
            g = -gcd.coefficients[0]
        if gcd.deg() == 2:
            #  print('2 roots')
            disc = gcd.coefficients[1] ** 2 - 4 * gcd.coefficients[0] * gcd.coefficients[2] % q
            disc = (modular_sqrt(disc, q)) % q
            g = (-gcd.coefficients[1] - disc) * pow(2 * gcd.coefficients[2], -1, q)
        # if gcd.deg()==l+1:
        # print('l+1 roots')
        #  g = gcd.findroot()
        Dg = g * Phi.partial_x(g, j) % q
        Dj = j * Phi.partial_y(g, j) % q
        s = int(12 / math.gcd(l - 1, 12))
        Delta_l = (g ** (int(12 / s))) * Delta * pow(l ** 12, -1, q) % q
        if Dj == 0:
            print("Dj==0")
            E4_l = E_4 * pow(l ** 2, -1, q) % q
            atilde = -3 * (l ** 4) * E4_l % q
            jtilde = (E4_l ** 3) * pow(Delta_l, -1, q) % q
            btilde = (2 * (l ** 6)) * (modular_sqrt((((jtilde - (1728)) * Delta_l)).n, q)) % q
            p1 = 0
        else:
            E_2 = -12 * E_6 * Dj * pow(s * E_4 * Dg, -1, q) % q
            gprime = -s * pow(12, -1, q) * E_2 * g % q
            jprime = -(E_4 ** 2) * E_6 * pow(Delta, -1, q) % q
            E_0 = E_6 * pow(E_4 * E_2, -1, q) % q

            Dgprime = gprime * Phi.partial_x(g, j) + g * (
                        gprime * Phi.partial_xx(g, j) + jprime * Phi.partial_xy(g, j)) % q
            Djprime = jprime * Phi.partial_y(g, j) + j * (
                        jprime * Phi.partial_yy(g, j) + gprime * Phi.partial_xy(g, j)) % q

            E0_prime = ((-s * pow(12, -1, q)) * Dgprime - E_0 * Djprime) * pow(Dj, -1, q) % q
            E4_l = pow(l ** 2, -1, q) * (E_4 - E_2 * (
                        12 * E0_prime * pow(E_0, -1, q) + (6 * (E_4 ** 2) * pow(E_6, -1, q)) - (
                            4 * E_6 * pow(E_4, -1, q))) + E_2 ** 2) % q
            jtilde = (E4_l ** 3) * pow(Delta_l, -1, q) % q
            print("Cannonical Jtilde", jtilde)
            f = (l ** s) * pow(g, -1, q)
            fprime = s * E_2 * f * pow(12, -1, q) % q
            Dgstar = Phi.partial_x(f, jtilde) % q
            Djstar = Phi.partial_y(f, jtilde) % q
            jtilde_prime = -fprime * Dgstar * pow(l * Djstar, -1, q) % q
            E6_l = -E4_l * jtilde_prime * pow(jtilde, -1, q) % q
            atilde = -3 * (l ** 4) * E4_l % q
            btilde = -2 * (l ** 6) * E6_l % q
            p1 = -l * E_2 % q
        d = int((l - 1) / 2)
        Ck, Ck_tilde = self.WeierstrassCoeffecients(d, a,b, atilde, btilde)
        coeffs = self.KernelPolyCoefficients2(p1, d, Ck, Ck_tilde)
        return self.P(coeffs)

    def WeierstrassCoeffecients(self, d, a, b, atilde, btilde):
        q = self.p
        '''compute coefficients of weierstrass p function for E and E-tilde'''
        Ck = [-a * pow(5, -1, q) % q, -b * pow(7, -1, q) % q]
        Ck_tilde = [-atilde * pow(5, -1, q) % q, -btilde * pow(7, -1, q) % q]
        for k in range(3, d + 1):
            coeffk = 3 * pow((k - 2) * (2 * k + 3), -1, q) % q
            c_k = coeffk * sum(Ck[j - 1] * Ck[k - 2 - j] for j in range(1, k - 1)) % q
            c_ktilde = coeffk * sum(Ck_tilde[j - 1] * Ck_tilde[k - 2 - j] for j in range(1, k - 1)) % q
            Ck.append(c_k)
            Ck_tilde.append(c_ktilde)
        return Ck, Ck_tilde

    def KernelPolyCoefficients(self, p1, d, Ck, Ck_tilde):
        '''Compute Coefficients for factor of division polynomial that
        vanhishes on eigenspace of the frobenius endomorphism'''
        Ak = [-p1 * pow(2, -1, self.p)]
        for j in range(1, d):
            Ak.append(
                -(Ck_tilde[j - 1] - (2 * d + 1) * Ck[j - 1]) * pow(((2 * j + 1) * (2 * j + 2)), -1, self.p) % self.p)
        C_series = [self.P.one]

        f = self.P(Ak)
        A = self.P([1] + Ak)
        s = self.P(Ck)
        C = self.P([0] + Ck)

        C_series.append(C)
        j = 2
        g = f
        h = self.P(Ck)
        while j <= d:
            g = g.mul_upto(f, d - j)
            q = self.P([pow((mpz(math.factorial(j))), -1, self.p)]) * g.shift(j)
            A = A + q
            h = h.mul_upto(s, d - j)
            C_series.append(h.shift(j))
            j = j + 1
        Coeffs = [0] * (d - 1) + [-p1 * pow(2, -1, self.p)] + [1]
        for i in range(1, d + 1):
            fl = A.coefficients[i]
            for k in range(1, i + 1):
                c = sum(mpz(math.factorial(d - i + k) / (math.factorial(k - j) * math.factorial(d - i + j))) % self.p *
                        C_series[k - j].coefficients[j] for j in range(0, k)) % self.p
                fl = fl - c * Coeffs[d - i + k] % self.p
            Coeffs[d - i] = fl
        return Coeffs

    def KernelPolyCoefficients2(self, p1, d, Ck, Ck_tilde):
        '''Compute Coefficients for factor of division polynomial that
        vanhishes on eigenspace of the frobenius endomorphism'''
        Ak = [-p1 * pow(2, -1, self.p)]
        Ak = Ak + [
            (-(Ck_tilde[j - 1] - (2 * d + 1) * Ck[j - 1]) * pow(((2 * j + 1) * (2 * j + 2)), -1, self.p) % self.p)
            for j in range(1, d)]
        C_series = [self.P.one]
        multiplier = self.P(Ak)
        A = self.P.one + self.P(Ak).shift(1)
        s = self.P(Ck)
        C = self.P([0] + Ck)
        C_series.append(C)
        j = 2
        g = multiplier
        h = s
        mod = self.P.x.shift(d)
        while j <= d:
            modj = self.P.x.shift(d - j)
            g = g * multiplier % modj
            q = self.P([pow((mpz(math.factorial(j))), -1, self.p)]) * g.shift(j)
            A = A + q
            h = h * s % modj
            C_series.append(h.shift(j))
            j = j + 1
        Coeffs = [0] * (d - 1) + [-p1 * pow(2, -1, self.p) % self.p] + [1]
        for i in range(1, d + 1):
            fl = A.coefficients[i]
            for k in range(1, i + 1):

                    vec = [mpz(math.factorial(d - i + k) // (math.factorial(k - j) * math.factorial(d - i + j))) % self.p *
                    C_series[k - j].coefficients[j] for j in range(0, k)]
                    c = sum(vec)%self.p
                    fl = fl - c * Coeffs[d - i + k] % self.p
            Coeffs[d - i] = fl
        return Coeffs


    def DivisionPoly(self, n,A,B):
        y2 = self.P([B, A, 0, 1])
        y4 = y2 * y2
        phi0 = self.P([0])
        phi1 = self.P([1])
        phi2 = self.P([2])
        phi3 = self.P([-(A ** 2), 12 * B, 6 * A, 0, 3])
        phi4 = self.P([-4 * (A ** 3) - 32 * (B ** 2), -16 * A * B, -20 * (A ** 2), 80 * B, 20 * A, 0, 4])
        Phi = [phi0, phi1, phi2, phi3, phi4]

        if n < 5:
            return Phi[0:n + 1]
        i = 5
        while i <= n:
            if i % 4 == 1:
                m = int((i - 1) / 2)

                new = (y4 * (Phi[m + 2] * (Phi[m] ** 3))) - (Phi[m - 1] * (Phi[m + 1] ** 3))
                Phi.append(new)

            if i % 4 == 3:
                m = int((i - 1) / 2)
                new = (Phi[m + 2] * (Phi[m] ** 3)) - (y4 * (Phi[m - 1] * (Phi[m + 1] ** 3)))
                Phi.append(new)

            if i % 2 == 0:
                m = i // 2
                r = self.P([1*pow(2,-1,self.p)]) * Phi[m] * (
                            (Phi[m + 2] * (Phi[m - 1] ** 2)) - (Phi[m - 2] * (Phi[m + 1] ** 2)))
                Phi.append(r)

            i = i + 1
        return Phi
