from Euclidean import *
from gmpy2 import *


def gen_base_point(curve, a):
    """returns a class Point for creating points and doring arithmetic  on an elliptic curve, where the point coordinates are 
    quotients of polynomials"""
    class Point:
        def __init__(self, x, y, z=1):
            self.x = x
            self.y = y
            self.z = z

        def __add__(self, other):
            if other.z == 0:
                return self
            if self.x == other.x:
                if self.y == other.y:
                    return self.double()
            return self.add(other)

        def frobenius(self, q):
            fraction = self.x.__class__
            quotient = fraction.quotient
            x_new = self.x ** q
            y_new = (fraction(quotient(curve)) ** ((q - 1) // 2)) * self.y
            return Point(x_new, y_new)

        def __eq__(self, other):
            if self.x == other.x and self.y == other.y:
                return True
            return False

        def add(self, other):
            fraction = self.x.__class__
            quotient = fraction.quotient
            m_num = self.y - other.y
            m_den = self.x - other.x
            m = m_num / m_den
            m2 = (m ** 2) * fraction(quotient(curve))
            x_new = m2 - other.x - self.x
            y_new = -m * (x_new - other.x) - other.y
            return Point(x_new, y_new)

        def double(self):
            fraction = self.x.__class__
            quotient = fraction.quotient
            ring = quotient.ring
            three = fraction(quotient(ring([3])))
            a_img = fraction(quotient(ring([a])))
            m_num = three * (self.x ** 2) + a_img
            two = fraction(quotient(ring([2])))
            m = (m_num / (two * self.y))
            m2 = m ** 2 / fraction(quotient(curve))
            x_new = m2 - two * self.x
            y_new = (-m / fraction(quotient(curve)) * (x_new - self.x)) - self.y
            return Point(x_new, y_new)

        def __neg__(self):
            return Point(self.x, -self.y)

        def multiply(self, n):
            m = abs(n)
            point = self
            i = 1
            while i < m:
                point = point + self
                i = i + 1
            if n > 0:
                return point
            return -point

        def multiply_v2(self, n):
            m = abs(n)
            point = Point.O
            point_q = self
            while m > 0:
                if m % 2 == 1:
                    point = point_q+point
                m = m // 2
                point_q = point_q.double()
            if n >= 0:
                return point
            return Point(point.x, -point.y)

    Point.O = Point(0, 1, 0)
    return Point
