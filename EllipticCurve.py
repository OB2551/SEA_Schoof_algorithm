from PolynomialsV2 import *
from QuotientRingV2 import *
from Euclidean import *
from ClassicalModularPolynomialsV2 import *
from PointPolynomialArithmeticV2 import *
from CanonicalModularPolynomialsV2 import *
from gmpy2.gmpy2 import mpz
from FracV2 import *
import math
import itertools
import random as rand
from gmpy2 import *

primes_list = get_primes(1000)


def point_factory(a, b, p):
    '''returns a class that generates random points on the curve y^2=x^3+ax+b over the field F_p'''
    class CurvePoint:
        def __init__(self, *args):
            if not args:
                while True:
                    x = rand.randint(1, p - 1)
                    m = (pow(x, 3, p) + a * x + b) % p
                    if m == 0:
                        self.x, self.y, self.z = x, 0, 1
                        break
                    y = modular_sqrt(m, p)
                    if y != 0:
                        self.x, self.y, self.z = x, y, 1
                        break
            else:
                self.x, self.y, self.z = args

        def double(self):
            if self.z == 0:
                return self
            if self.y == 0:
                return CurvePoint(0, 1, 0)
            else:
                m = ((3 * (self.x ** 2) + a) * pow((2 * self.y), p - 2, p)) % p
                x = (m ** 2 - self.x - self.x) % p
                y = (-m * (x - self.x) - self.y) % p
                return CurvePoint(x, y, 1)

        def add(self, other):

            if self.z == 0:
                return other
            if other.z == 0:
                return self
            if self.x == other.x and self.y == other.y:
                return self.double()
            if self.x == other.x:
                return CurvePoint(0, 1, 0)
            else:
                m = ((self.y - other.y) * pow((self.x - other.x), -1, p)) % p
                x = m ** 2 - other.x - self.x
                y = m * (other.x - x) - other.y
                return CurvePoint(x % p, y % p, 1)

        def scalar(self, m):
            point = self
            n = abs(m)
            res = CurvePoint(0, 1, 0)
            while n > 0:
                if n % 2 == 1:
                    res = res.add(point)
                point = point.double()
                n = n // 2
            if m < 0 and res.z != 0:
                return CurvePoint(res.x % p, -res.y % p, 1)
            return res

        def sub(self, other):
            if other.z == 0:
                return self
            r = CurvePoint(other.x % p, -other.y % p, 1)
            return self.add(r)

    return CurvePoint


class EllipticCurve:
    def __init__(self, p1, a, b):
        p = mpz(p1)
        self.P = polynomials_over(p)
        self.a = mpz(a) % p
        self.b = mpz(b) % p
        self.EC = self.P([self.b, self.a, mpz(0), mpz(1)])
        self.p = mpz(p)
        if self.discriminant == 0:
            print("Singular Curve")

    def j_inv(self):
        """compute j-invariant of Elliptic Curve"""
        a = self.a
        b = self.b
        j = 1728 * 4 * (a ** 3) * pow((4 * (a ** 3) + 27 * (b ** 2)), -1, self.p)
        return j % self.p
    
    def discriminant(self):
        """compute discriminant of elliptic curve"""
        return (4*(self.a**3)+27*(self.b**2))%self.p

    def division_poly(self, n):
        """Compute up to the nth division polynomial required for Schoof algorithm"""
        a = self.a
        b = self.b
        y2 = self.P([b, a, 0, 1])
        y4 = y2 * y2
        phi0 = self.P([0])
        phi1 = self.P([1])
        phi2 = self.P([2])
        phi3 = self.P([-(a ** 2), 12 * b, 6 * a, 0, 3])
        phi4 = self.P([-4 * (a ** 3) - 32 * (b ** 2), -16 * a * b, -20 * (a ** 2), 80 * b, 20 * a, 0, 4])
        phi = [phi0, phi1, phi2, phi3, phi4]

        if n < 5:
            return phi[0:n + 1]
        i = 5
        while i <= n:
            if i % 4 == 1:
                m = int((i - 1) / 2)
                new = (y4 * (phi[m + 2] * (phi[m] ** 3))) - (phi[m - 1] * (phi[m + 1] ** 3))
                phi.append(new)

            if i % 4 == 3:
                m = int((i - 1) / 2)
                new = (phi[m + 2] * (phi[m] ** 3)) - (y4 * (phi[m - 1] * (phi[m + 1] ** 3)))
                phi.append(new)

            if i % 2 == 0:
                m = i // 2
                r = self.P([pow(2, -1, self.p)]) * phi[m] * (
                            (phi[m + 2] * (phi[m - 1] ** 2)) - (phi[m - 2] * (phi[m + 1] ** 2)))
                phi.append(r)

            i = i + 1
        return phi

    def odd_l_subroutine(self, phi, l):
        """compute the trace of Frobenius modulo l for odd primes l"""
        phi_l = phi[l]
        #construct field of fractions of F_p[x]/<psi_l(x)> and a base point with 
        #coordinates in the fraction field
        quotient = gen_quotient_ring(phi_l, self.P)
        ltor = gen_fraction_field(quotient)
        base = gen_base_point(self.EC, self.a)
        point = base(ltor.x, ltor.one)
        
        q = self.p
        ql = q % l
        if ql > l // 2:
            ql = ql - l
        q1, q2 = point.frobenius(q ** 2), point.multiply_v2(ql)
        
        #catch edge case when x coordinates of q1,q2 are equal
        if q1.x == q2.x:
            w = modular_sqrt(ql, l)
            if w == 0:
                return 0
            else:
                pw = point.multiply_v2(w)
                if (point.x ** q - pw.x).num.polynomial.gcd(phi_l).deg() == 0:
                    return 0
                else:
                    if (point.y ** q - pw.y).num.polynomial.gcd(phi_l).deg() == 0:
                        return (-2 * w) % l
                    else:
                        return (2 * w) % l

        point_q = q1 + q2
        bound = int((l - 1) / 2)
        j = 1
        j_point = point.frobenius(q)
        adder = j_point
        while j <= bound:
            if j_point.x == point_q.x:
                break
            j += 1
            j_point += adder
        if j_point.y == point_q.y:
            return j
        if j_point.y == -point_q.y:
            return -j
        else:
            print("fail")

    def a_mod2(self):
        """determine the trace of Frobenius moduluo 2"""
        q = self.p
        gcd_mod2 = (self.P.x.powmod(q, self.EC) - self.P.x).gcd(self.EC)
        if gcd_mod2.deg() == 0:
            return 1
        return 0

    def schoof(self):
        """schoofs algorithm"""
        #get primes required and division polynomials
        primes_list = primes(self.p)
        max_l = primes_list[-1]
        phi = self.division_poly(max_l)
        
        c1, c2 = [], []
        x, y = self.a_mod2(), 2
        c1.append(x)
        c2.append(2)
        
        #determine trace modulo l
        for l in primes_list[1:]:
            x = self.odd_l_subroutine(phi, l)
            c1.append(x)
            c2.append(l)
            
        #use chinese remainder theorem to recover trace
        a, modulus = crt(c1, c2)
        if a > 2 * math.sqrt(self.p):
            a = a - modulus
            
        return "#E(F_p) = " + str(self.p + 1 - a)

    def eisensteins(self):
        e_4 = -48 * self.a
        e_6 = 864 * self.b
        return e_4 % self.p, e_6 % self.p

    def j_inv_prime(self):
        e_4, e_6 = self.eisensteins()
        j_prime = (-e_6 * pow(e_4, self.p - 2, self.p)) * self.j_inv()
        return j_prime % self.p

    def prime_check(self, modular_polynomial):
        """Check whether l is an Elkies or Atkin prime"""
        pol = (self.P.x.powmod(self.p, modular_polynomial) - self.P.x)
        divisor = pol.gcd(modular_polynomial)
        lead = divisor.coefficients[-1]
        divisor = divisor * self.P([pow(lead, -1, self.p)])
        return divisor

    def isogenous_curve(self, l, modular_polynomial, divisor):
        """For an Elkies prime l, determine equaiton of an l-isogenous curve"""
        j = self.j_inv()
        q: mpz = self.p
        if divisor.deg() == 1:
            j_tilde = -divisor.coefficients[0] % q
        else:
            delta = divisor.coefficients[1] ** 2 - 4 * divisor.coefficients[0] * divisor.coefficients[2] % q
            disc = modular_sqrt(delta, q) % q
            j_tilde = (-divisor.coefficients[1] - disc) * pow(2 * divisor.coefficients[2], -1, q) % q
        j_prime = self.j_inv_prime()
        phi_x = modular_polynomial.partial_x(j, j_tilde) % q
        phi_y = modular_polynomial.partial_y(j, j_tilde) % q
        j_tilde_prime = -(j_prime * phi_x) * pow(l * phi_y, -1, q) % q
        a_tilde = -(j_tilde_prime ** 2) * pow(48 * j_tilde * (j_tilde - 1728), -1, q) % q
        b_tilde = -(j_tilde_prime ** 3) * pow(864 * (j_tilde ** 2) * (j_tilde - 1728), -1, q) % q
        return j, j_prime, j_tilde, j_tilde_prime, a_tilde, b_tilde

    def elkies_procedure(self, l, phi_l, divisor):
        """Elkies method to recover divisor f_l(x) of division polynomial psi_l(x)"""
        q = self.p
        j, j_prime, j_tilde, j_tilde_prime, a_tilde, b_tilde = self.isogenous_curve(l, phi_l, divisor)
        t1 = j_prime ** 2 * phi_l.partial_xx(j, j_tilde) % q
        t2 = 2 * l * j_prime * j_tilde_prime * phi_l.partial_xy(j, j_tilde) % q
        t3 = (l ** 2) * (j_tilde_prime ** 2) * phi_l.partial_yy(j, j_tilde) % q
        e4_ql = -48 * a_tilde % q
        e6_ql = 864 * b_tilde % q
        e4, e6 = self.eisensteins()
        vii14 = -(t1 + t2 + t3) * pow(j_prime * phi_l.partial_x(j, j_tilde), -1, q) % q
        t1 = pow(2, -1, q) * l * vii14 % q
        t2 = pow(4, -1, q) * l * ((e4 ** 2) * pow(e6, -1, q) - l * ((e4_ql ** 2) * pow(e6_ql, -1, q))) % q
        t3 = pow(3, -1, q) * l * (e6 * pow(e4, -1, q) - l * (e6_ql * pow(e4_ql, -1, q))) % q
       
        #coefficient of x^((l-1)/2-1)
        p1 = (t1 + t2 + t3) % q
        a_tilde = (l ** 4) * a_tilde % q
        b_tilde = (l ** 6) * b_tilde % q
        d = int((l - 1) / 2)
        c_k, c_k_tilde = self.weierstrass_coefficients(d, self.a, self.b, a_tilde, b_tilde)
        coefficients = self.kernel_poly_coefficients2(p1, d, c_k, c_k_tilde)
        return self.P(coefficients)

    def elkies_procedure_v2(self, l, phi, divisor):
        """Second version of method to recover divisor f_l(x) of lth division polynomial psi_l(x)
        , this method is used for canonical modular polynomials"""
        q = self.p
        e_4 = -self.a * pow(3, -1, q) % q
        e_6 = -self.b * pow(2, -1, q) % q
        delta = (e_4 ** 3 - e_6 ** 2) * pow(1728, -1, q) % q
        j = self.j_inv()
        
        if divisor.deg() == 1:
            g = -divisor.coefficients[0]
        else:
            disc = divisor.coefficients[1] ** 2 - 4 * divisor.coefficients[0] * divisor.coefficients[2] % q
            disc = (modular_sqrt(disc, q)) % q
            g = (-divisor.coefficients[1] - disc) * pow(2 * divisor.coefficients[2], -1, q)
        dg = g * phi.partial_x(g, j) % q
        dj = j * phi.partial_y(g, j) % q
        s = int(12 / math.gcd(l - 1, 12))
        delta_l = (g ** (int(12 / s))) * delta * pow(l ** 12, -1, q) % q
        
        if dj == 0:
            e4_l = e_4 * pow(l ** 2, -1, q) % q
            a_tilde = -3 * (l ** 4) * e4_l % q
            j_tilde = (e4_l ** 3) * pow(Delta_l, -1, q) % q
            b_tilde = (2 * (l ** 6)) * (modular_sqrt(((j_tilde - 1728) * delta_l).n, q)) % q
            p1 = 0
        else:
            e_2 = -12 * e_6 * dj * pow(s * e_4 * dg, -1, q) % q
            g_prime = -s * pow(12, -1, q) * e_2 * g % q
            j_prime = -(e_4 ** 2) * e_6 * pow(delta, -1, q) % q
            e_0 = e_6 * pow(e_4 * e_2, -1, q) % q

            dg_prime = g_prime * phi.partial_x(g, j) + g * (
                        g_prime * phi.partial_xx(g, j) + j_prime * phi.partial_xy(g, j)) % q
            dj_prime = j_prime * phi.partial_y(g, j) + j * (
                        j_prime * phi.partial_yy(g, j) + g_prime * phi.partial_xy(g, j)) % q

            e0_prime = ((-s * pow(12, -1, q)) * dg_prime - e_0 * dj_prime) * pow(dj, -1, q) % q
            e4_l = pow(l ** 2, -1, q) * (e_4 - e_2 * (
                        12 * e0_prime * pow(e_0, -1, q) + (6 * (e_4 ** 2) * pow(e_6, -1, q)) - (
                            4 * e_6 * pow(e_4, -1, q))) + e_2 ** 2) % q
            
            j_tilde = (e4_l ** 3) * pow(delta_l, -1, q) % q
            f = (l ** s) * pow(g, -1, q)
            f_prime = s * e_2 * f * pow(12, -1, q) % q
            dg_star = phi.partial_x(f, j_tilde) % q
            dj_star = phi.partial_y(f, j_tilde) % q
            j_tilde_prime = -f_prime * dg_star * pow(l * dj_star, -1, q) % q
            e6_l = -e4_l * j_tilde_prime * pow(j_tilde, -1, q) % q
            a_tilde = -3 * (l ** 4) * e4_l % q
            b_tilde = -2 * (l ** 6) * e6_l % q
            p1 = -l * e_2 % q
            
        d = int((l - 1) / 2)
        ck, ck_tilde = self.weierstrass_coefficients(d, self.a, self.b, a_tilde, b_tilde)
        coefficients = self.kernel_poly_coefficients2(p1, d, ck, ck_tilde)
        return self.P(coefficients)

    def weierstrass_coefficients(self, d, a, b, a_tilde, b_tilde):
        """compute weierstrass coefficients of the weierstrass p function for an elliptic curve E, 
        and l-isogenous elliptic curve E' """
        q = self.p
        '''compute coefficients of weierstrass p function for E and E-tilde'''
        ck = [-a * pow(5, -1, q) % q, -b * pow(7, -1, q) % q]
        ck_tilde = [-a_tilde * pow(5, -1, q) % q, -b_tilde * pow(7, -1, q) % q]
        for k in range(3, d + 1):
            coefficient_k = 3 * pow((k - 2) * (2 * k + 3), -1, q) % q
            c_k = coefficient_k * sum(ck[j - 1] * ck[k - 2 - j] for j in range(1, k - 1)) % q
            c_k_tilde = coefficient_k * sum(ck_tilde[j - 1] * ck_tilde[k - 2 - j] for j in range(1, k - 1)) % q
            ck.append(c_k)
            ck_tilde.append(c_k_tilde)
        return ck, ck_tilde

    def kernel_poly_coefficients2(self, p1, d, ck, ck_tilde):
        """Compute Coefficients for divisor of division polynomial that
            vanishes on eigenspace of the frobenius endomorphism"""
        ak = [-p1 * pow(2, -1, self.p)]
        ak = ak + [
            (-(ck_tilde[j - 1] - (2 * d + 1) * ck[j - 1]) * pow(((2 * j + 1) * (2 * j + 2)), -1, self.p) % self.p)
            for j in range(1, d)]
        c_series = [self.P.one]
        multiplier = self.P(ak)
        a_series = self.P.one + self.P(ak).shift(1)
        s = self.P(ck)
        c_prime = self.P([0] + ck)
        c_series.append(c_prime)
        j = 2
        g = multiplier
        h = s
        while j <= d:
            modj = self.P.x.shift(d - j)
            g = g * multiplier % modj
            q = self.P([pow((mpz(math.factorial(j))), -1, self.p)]) * g.shift(j)
            a_series += q
            h = h * s % modj
            c_series.append(h.shift(j))
            j = j + 1
        coefficients = [0] * (d - 1) + [-p1 * pow(2, -1, self.p) % self.p] + [1]
        for i in range(1, d + 1):
            fl = a_series.coefficients[i]
            for k in range(1, i + 1):
                vec = [
                    mpz(math.factorial(d - i + k) // (math.factorial(k - j) * math.factorial(d - i + j))) % self.p *
                    c_series[k - j].coefficients[j] for j in range(0, k)]
                c = sum(vec) % self.p
                fl = fl - c * coefficients[d - i + k] % self.p
            coefficients[d - i] = fl
        return coefficients

    def frobenius_order(self, l, phi):
        """compute the order of frobenius r in the case of an atkin prime l so that
            the lth modular polynomial splits over F_p^r"""
        b = modular_sqrt(self.p, l)
        if b == 0:
            b = -1
        else:
            b = 1
        k = 0
        h = self.P.x
        
        #special edge cases where r can be determined easily
        if l in [37, 61, 73]:
            if ((h.powmod(self.p ** 2, phi)) - self.P.x).gcd(phi).deg() == l + 1:
                return 2
            else:
                if b == 1:
                    return (l + 1) // 2
                if b == -1:
                    return l + 1
                
        for i in range(1, l + 2):
            if (l + 1) % i == 0 and ((-1) ** ((l + 1) // i) == b):
                for c in range(1, i - k + 1):
                    h = h.powmod(self.p, phi)
                k = i
                g = (h - self.P.x).gcd(phi)
                if g.deg() == l + 1:
                    return i

    def atkin(self, l, r):
        """Process an Atkin prime"""
        pl = polynomials_over(l)
        d = get_non_square(l)
        irreducible_modulus = pl([-d, 0, 1])
        ql = gen_quotient_ring(irreducible_modulus, pl)
        divs = list(factors(l ** 2 - 1))
        divs.sort()
        divs = divs[1:-1]
        while True:
            # find a generator of multiplicative group F_l^2
            a = rand.randint(0, l)
            g = ql(pl([a, 1]))  # sqrt{d}-a
            t = 1
            for d in divs:
                if g ** ((l ** 2 - 1) / d) == ql.one:
                    t = 0
                    break
            if t == 1:
                generator = g
                break

        # compute candidates for rth roots of unity in F_l
        g = generator ** ((l ** 2 - 1) / r)
        s = ql.one
        unities = []
        k = 0
        rel_primes = get_rel_primes(r)
        for i in rel_primes:
            s = s * (g ** (i - k))
            k = i
            unities.append(s)
        traces = []
        for gamma_r in unities:
            g1 = gamma_r.polynomial.coefficients[0]
            x1_sq = (self.p * (g1 + 1) * pow(2, -1, l)) % l
            if x1_sq == 0:
                traces.append(0)
            x1 = modular_sqrt(x1_sq, l)
            if x1 != 0:
                traces.append(2 * x1 % l)
                traces.append(-2 * x1 % l)
        return list(set(traces))

    def atkin_traces(self, atkin, number_of_traces):
        """split list of Atkin primes and their possible traces into two approximately equal sets"""
        l = len(atkin[0])
        i = 0
        traces_count = 0
        s1 = []
        s1_traces = []
        if l == 1:
            s1 = atkin[1]
            s1_traces = atkin[0]
            t1 = [(trace, s1[0]) for trace in s1_traces]
            return t1, []
        while i < l:
            if traces_count + len(atkin[0][i]) > number_of_traces / 2:
                s2 = atkin[1][i:]
                s2_traces = atkin[0][i:]
                break
            traces_count += len(atkin[0][i])
            s1.append(atkin[1][i])
            s1_traces.append(atkin[0][i])
            i = i + 1
        s1_traces_product = itertools.product(*s1_traces)
        s2_traces_product = itertools.product(*s2_traces)

        t1 = [crt(trace, s1) for trace in s1_traces_product]
        t2 = [crt(trace, s2) for trace in s2_traces_product]

        return t1, t2

    def baby_step_giant_step(self, elkies, atkin1, atkin2):
        """Algorithm to find a match for trace after sorting possible trace values"""
        t3, m3 = crt(elkies[0], elkies[1])
        points = point_factory(self.a, self.b, self.p)
        traces1, traces2 = atkin1, atkin2
        if not traces2:
            m1 = traces1[0][1]
            point = points()
            for t1 in traces1[0][0]:
                t, _ = crt([t1, t3], [m1, m3])
                if abs(t) > 2 * math.sqrt(self.p):
                    t = t - _
                if point.scalar(self.p + 1 - t).z == 0:
                    return t
        m1, m2 = traces1[0][1], traces2[0][1]

        r1, r2 = [], []
        r1, r2 = sorted(r1, key=abs), sorted(r2, key=abs)

        for t1, m1 in traces1:
            r = ((t1 - t3) * pow(m2 * m3, -1, m1)) % m1
            if r > m1 / 2:
                r1.append(r - m1)
            else:
                r1.append(r)
        for t2, m2 in traces2:
            r = ((t2 - t3) * pow(m1 * m3, -1, m2)) % m2
            r2.append(r)
            r2.append(r - m2)
        i = 1
        #keep iterating unitl we find a match. we expect to find a match quickly, so if no match is found then something
        #has gone wrong in the previous steps, so break out of loop after an arbitrary 40 iterations
        while i < 40:
            point = points()
            point_q = point.scalar(self.p + 1 - t3)
            m3point = point.scalar(m3)
            m1m3point = m3point.scalar(m1)
            m2m3point = m3point.scalar(m2)
            point_q_r1 = {}

            for r1_i in r1:
                q_r1_i = point_q.sub(m2m3point.scalar(r1_i))
                if q_r1_i.z != 0:
                    point_q_r1[q_r1_i.x] = q_r1_i.y, r1_i
                if q_r1_i.z == 0:
                    point_q_r1["O"] = r1_i

            for r2_i in r2:
                r2_im1m3point = m1m3point.scalar(r2_i)
                if r2_im1m3point.z != 0:
                    if r2_im1m3point.x % self.p in point_q_r1:
                        y_r_1_i, r_1_i = point_q_r1[r2_im1m3point.x]
                        
                        if (y_r_1_i - r2_im1m3point.y) % self.p == 0:
                            t = t3 + m3 * (m2 * r_1_i + m1 * r2_i)
                            return t
                        
                        else:
                            t = t3 + m3 * (m2 * r_1_i - m1 * r2_i)
                            return t
                        
                if r2_im1m3point.z == 0:
                    if "O" in point_q_r1:
                        r_1_i = point_q_r1["O"]
                        t = t3 + m3 * (m2 * r_1_i + m1 * r2_i)
                        return t
            i += 1
        print("Failed to find match at BSGS step")

    def sea(self, canonical=True, atkin_filter=True):
        """SEA algorithm. Options availabl to filter Atkin primes and use canonical or classical moudlar polynomials"""
        x, y = self.a_mod2(), 2
        prod = 2
        e1 = [x]
        e2 = [y]
        a1 = []
        a2 = []
        number_of_traces = 0
        counter = 1
        elkies_count = 0
        atkin_count = 0
        val = 4 * math.sqrt(self.p)
        while prod <= val:
            l = primes_list[counter]
            
            #load modular polynomimal
            if canonical and l <= 197:
                phi_l = CanonicalModularPolynomial(self.p, l)
            else:
                phi_l = ModularPolynomial(self.p, l)
                
            j = self.j_inv()
            phi = phi_l.eval_at_y(j)
            divisor = self.prime_check(phi)
            
            #pathalogical Elkies prime, which we ignore
            if divisor.deg() == l + 1:
                counter += 1
                continue
                
            #Atkin prime    
            if divisor.deg() == 0: 
                if atkin_filter:
                    if l > 30 and (l != 73 or l != 61):
                        counter += 1
                        print(str(l) + ' rejected')
                        continue
                atkin_count += 1
                r = self.frobenius_order(l, phi)
                traces = self.atkin(l, r)
                number_of_traces += len(traces)
                a1.append(traces)
                a2.append(l)
                counter += 1
                prod = l * prod
                
            #Elkies prime
            if divisor.deg() > 0:
                elkies_count += 1
                
                #find divisor f_l(x) of division polynomial psi_l(x)
                if canonical and l <= 197:
                    f_l = self.elkies_procedure_v2(l, phi_l, divisor)
                else:
                    f_l = self.elkies_procedure(l, phi_l, divisor)

                #generate fraction field of F_p[x]/<f_(x)> and base point on elliptic curve
                #with coordinates in this fraction field
                ltorring = gen_quotient_ring(f_l, self.P)
                ltor = gen_fraction_field(ltorring)
                base = gen_base_point(self.EC, self.a)
                point = base(ltor.x, ltor.one)
                frob_point = point.frobenius(self.p)
                
                
                if divisor.deg() == 2:
                    bound = (l - 1) // 2
                    j = 1
                    point_j = point
                    light = 0
                    
                    #find eigenvalue of frobenius endomorphism
                    while j <= bound:
                        if (point_j.x - frob_point.x).num.polynomial.gcd(f_l).deg() > 0:
                            light = 1
                            break
                        point_j += point
                        j = j + 1
                    if light == 0:
                        print('no eigenvalue found, skipping')
                        counter += 1
                        continue
                    if (point_j.y - frob_point.y).num.polynomial.gcd(f_l).deg() > 0:
                        t = (j + self.p * pow(j, l - 2, l)) % l
                       
                    else:
                        t = (-j + self.p * pow(-j, l - 2, l)) % l
                        
                        
                #catch edge case
                else:
                    w = modular_sqrt(self.p, l)
                    wp = point.multiply(w)
                    if (wp.y - frob_point.y).num.polynomial.gcd(f_l).deg() > 0:
                        t = (2 * w % l)
                        
                    else:

                        t = (-2 * w % l)
                      

                e1.append(t)
                e2.append(l)
                counter += 1
                prod = l * prod
                
        elkies = e1, e2
        
        #if no atkin primes used
        if not a1:
            t, m = crt(e1, e2)
            if abs(t) > 2 * math.sqrt(self.p):
                return "E(F_p)=" + str(self.p + 1 - (t - m))
            return "E(F_p)=" + str(self.p + 1 - t)
        
        #otherwise use match sort algorithm to find trace t
        atkin = a1, a2
        atkin1, atkin2 = self.atkin_traces(atkin, number_of_traces)
      
        t = self.baby_step_giant_step(elkies, atkin1, atkin2)
        
        return "#E(F_p) = " + str(self.p + 1 - t)

    
    
    def baby_step(self):
        """baby step giant step algorithm to compute the number of points on the curve 
        by finding divisors of #E(F_p)"""
        
        m = rand.randint(math.floor(math.sqrt(math.sqrt(self.p))), math.ceil(math.sqrt(math.sqrt(self.p)) + 10))
        points = point_factory(self.a, self.b, self.p)
        l = 1
        thresh = 4 * math.sqrt(self.p)
        
        while l < thresh:
            point = points()
            #find the order of the random point by a baby step giant step collision algorithm
            j_ps = [(point.scalar(j), j) for j in range(0, m + 1)]
            j_psx = [point[0].x for point in j_ps]
            j_psy = [(point[0].y, point[1]) for point in j_ps]
            dic = dict(zip(j_psx, j_psy))
            point_q = point.scalar(self.p + 1)
            two_m_point = point.scalar(2 * m)
            k = 0
            order = 0
            
            while k < self.p + 2 * math.sqrt(self.p):
                z = point_q.add(two_m_point.scalar(k))
                if z.x in dic.keys():
                    if z.y == dic[z.x][0]:
                        j = dic[z.x][1]
                    else:
                        j = - dic[z.x][1]
                    order = self.p + 1 + (2 * m * k) - j
                    break
                k += 1
                
            if order == 0:
                continue
            primes = prime_factors(order)
            for prime in primes:
                if point.scalar(order // prime).z % self.p == 0:
                    order = order // prime
                    
            l = (order * l) // math.gcd(order, l)

        #once l is large enough there will be only one value in this interval with n%l == 0
        for n in range(int(self.p + 1 - 2 * math.sqrt(self.p)), int(self.p + 1 + 2 * math.sqrt(self.p))):
            if n % l == 0:
                return "#E(F_p) = " + str(n)

            
    def find_point(self, divs):
        """finds a point of maximal prime order given a sorted list of prime divisors of #E(F_p)"""
        points = point_factory(self.a, self.b, self.p)
        while True:
            indicator = 0
            point = points()
            for d in divs[:-1]:
                if point.scalar(d).z == 0:
                    indicator = 1
                    break
            if indicator == 0:
                if point.scalar(divs[-1]).z == 0:
                    return P
            else:
                continue



