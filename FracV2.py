from QuotientRingV2 import *


def gen_fraction_field(quotient):
    """generates a fraction field class given a quotient ring class. Stores fractions as a pair (a,b)
    where a is the numerator and b the denominator"""
    class Fraction:
        def __init__(self, num, dem=None):
            self.num = num
            if dem is None:
                self.dem = quotient.one
            else:
                self.dem = dem

        def __add__(self, other):
            new_num = self.num * other.dem + other.num * self.dem
            new_dem = self.dem * other.dem
            return Fraction(new_num, new_dem)

        def __sub__(self, other):
            new_num = self.num * other.dem - other.num * self.dem
            new_dem = self.dem * other.dem
            return Fraction(new_num, new_dem)

        def __mul__(self, other):
            new_num = self.num * other.num
            new_dem = self.dem * other.dem
            return Fraction(new_num, new_dem)

        def inverse(self):
            return Fraction(self.dem, self.num)

        def __truediv__(self, other):
            new_num = self.num * other.dem
            new_dem = self.dem * other.num
            return Fraction(new_num, new_dem)

        def __eq__(self, other):
            if self.num * other.dem == self.dem * other.num:
                return True
            return False

        def __neg__(self):
            return Fraction(-self.num, self.dem)

        def __pow__(self, n):
            return Fraction(self.num ** n, self.dem ** n)

        def __repr__(self):
            return str(self.num) + "/" + str(self.dem)

    Fraction.quotient = quotient
    Fraction.zero = Fraction(quotient.zero, quotient.one)
    Fraction.one = Fraction(quotient.one, quotient.one)
    Fraction.x = Fraction(quotient.x)
    return Fraction
