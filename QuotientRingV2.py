from Euclidean import *


def gen_quotient_ring(modulus, ring):
    """returns a Quotient ring class  a polynomial ring modulo the specified modulus, ie something of the form Z[x]/<f(x)>.
    Operators overloaded."""
    
    class QuotientRing:
        def __init__(self, polynomial):
            self.polynomial = polynomial % modulus

        def __add__(self, other):
            return QuotientRing(self.polynomial + other.polynomial)

        def __sub__(self, other):
            return QuotientRing(self.polynomial - other.polynomial)

        def __mul__(self, other):
            return QuotientRing(self.polynomial * other.polynomial)

        def __truediv__(self, other):
            return QuotientRing(self.polynomial * other.inverse())

        def inverse(self):
            x, y, z = extended_euclidean_algorithm2(self.polynomial, self.modulus)
            if z.deg() == 0:
                return x
            else:
                return self.ring.zero

        def __neg__(self):
            return QuotientRing(-self.polynomial)

        def __repr__(self):
            return str(self.polynomial)

        def __pow__(self, n):
            q = self
            r = QuotientRing.one
            while n > 0:
                if n % 2 == 1:
                    r = q * r
                q = (q * q)
                n = n // 2
            return r

        def __eq__(self, other):
            if self.polynomial == other.polynomial:
                return True
            return False

        def gcd(self, other):
            z = self.polynomial.gcd(other.polynomial)
            return QuotientRing(z)

    QuotientRing.ring = ring
    QuotientRing.modulus = modulus
    QuotientRing.one = QuotientRing(ring.one)
    QuotientRing.zero = QuotientRing(ring.zero)
    QuotientRing.x = QuotientRing(ring.x)
    return QuotientRing
