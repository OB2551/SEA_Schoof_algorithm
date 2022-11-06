import numpy as np
import sympy as sy
import math
from functools import reduce
import operator
import gmpy2
import random as rand


#file of various utility functions and operators on euclidean rings
def prod(iterable):
    return reduce(operator.mul, iterable, 1)


def gcd(a, b):
    if abs(a) < abs(b):
        return Gcd(b, a)

    while abs(b) > 0:
        _, r = divmod(a, b)
        a, b = b, r

    return a


def random_n_bit_prime(n):
    """get a random n bit prime"""
    lower, upper = 2 ** (n - 1), (2 ** n) - 1
    while True:
        m = rand.randint(lower, upper)
        prime = gmpy2.next_prime(m)
        if int(math.floor(math.log(prime, 2)) + 1) == n:
            return prime


def check_type(s):
    """check if int string contains any illegal characters"""
    alpha = "abcdefghijklmnopqrstuvwxyz%!$#^&*()-+~/'{}="
    for i in s:
        if i in alpha:
            return False
    return True


def polynomial_gcd(u, v):
    """for polynomials with coefficients of class Finite Field"""
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v

    while lesser_remainder != zero:
        quotient, next_remainder = divmod(larger_remainder, lesser_remainder)

        # Shift roles
        larger_remainder = lesser_remainder
        lesser_remainder = next_remainder

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    c = larger_remainder.coefficients[0]
    c = larger_remainder.Field.one / c
    p = larger_remainder.__class__([c])
    return p * larger_remainder


def extended_euclidean_algorithm_ints(d, e):
    """for integers"""

    a = d
    b = e
    if abs(b) == 0:
        return 1, 0, a

    x1, x2 = 0, 1
    while abs(b) > 0:
        q, r = divmod(a, b)

        x = x2 - q * x1

        a, b, x2, x1 = b, r, x1, x
    c = (a - x2 * d) // e
    return x2, c, a


def extended_euclidean_algorithm(u, v):
    """for polynomials"""
    one = u.__class__.one
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v
    larger_scalar, lesser_scalar = one, zero
    while lesser_remainder != zero:
        quotient, next_remainder = divmod(larger_remainder, lesser_remainder)
        next_scalar = larger_scalar - quotient * lesser_scalar

        # Shift roles
        larger_remainder, larger_scalar = lesser_remainder, lesser_scalar
        lesser_remainder, lesser_scalar = next_remainder, next_scalar

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    other_scalar = (larger_remainder - u * larger_scalar) // v
    c = larger_remainder.coefficients[0]
    c = larger_remainder.Field.one / c
    p = larger_remainder.__class__([c])
    return p * larger_scalar, p * other_scalar, p * larger_remainder


def extended_euclidean_algorithm2(u, v):
    """for polynomials"""
    one = u.__class__.one
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v
    larger_scalar, lesser_scalar = one, zero
    while lesser_remainder != zero:
        quotient, next_remainder = divmod(larger_remainder, lesser_remainder)
        next_scalar = larger_scalar - quotient * lesser_scalar

        # Shift roles
        larger_remainder, larger_scalar = lesser_remainder, lesser_scalar
        lesser_remainder, lesser_scalar = next_remainder, next_scalar

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    other_scalar = (larger_remainder - u * larger_scalar) // v
    c = larger_remainder.coefficients[0]
    p = u.__class__.modulus
    c = pow(c, -1, p)
    p = larger_remainder.__class__([c])
    return p * larger_scalar, p * other_scalar, p * larger_remainder


def crt(congruences1, congruences2):
    x = 0
    mod = prod(congruences2)

    for i in range(len(congruences1)):
        n_i = congruences1[i]
        m_i = congruences2[i]

        x = x + (mod // m_i) * pow(mod // m_i, -1, m_i) * n_i

    return int(x % mod), mod


def get_primes(n):
    numbers = set(range(n, 1, -1))
    primes_list = []
    while numbers:
        p = numbers.pop()
        primes_list.append(p)
        numbers.difference_update(set(range(p * 2, n + 1, p)))
    return primes_list


def primes(p):
    primes_list = get_primes(1000)
    product = 1
    i = 0
    j = -1
    while product < 4 * math.sqrt(p):
        if primes_list[i] != p:
            product = product * primes_list[i]
        else:
            j = i
        i += 1
    if j != -1:
        return primes_list[0:j] + primes_list[j + 1:i]
    return primes_list[0:i]


def get_non_square(p):
    """get a non quadratic residue mod p"""
    while True:
        a = np.random.randint(2, p)
        if modular_sqrt(a, p) == 0:
            return a


def phi(n):
    """euler totient function, count number of rel primes to n"""
    amount = 0
    for k in range(1, n + 1):
        if math.gcd(n, k) == 1:
            amount += 1
    return amount


def get_rel_primes(n):
    """Get relatively prime integers <= n"""
    s = []
    for i in range(1, n):
        if math.gcd(i, n) == 1:
            s.append(i)
    return s


def factors(n):
    """Factors of n"""
    return set(reduce(list.__add__,
                      ([i, n // i] for i in range(1, int(n ** 0.5) + 1) if n % i == 0)))


def prime_factors(n):
    """prime factors of n, with multiplicity"""
    primefactors = []
    d = 2
    while d * d <= n:
        while (n % d) == 0:
            primefactors.append(d)  # supposing you want multiple factors repeated
            n //= d
        d += 1
    if n > 1:
        primefactors.append(n)
    return primefactors


def modular_sqrt(a, p):
    """get modular square root of a mod p"""

    def legendre_symbol(a1, p1):
        """ Compute the Legendre symbol a|p using
            Euler's criterion. p is a prime, a is
            relatively prime to p (if p divides
            a, then a|p = 0)
            Returns 1 if a has a square root modulo
            p, -1 otherwise.
        """

        ls = pow(a1, (p1 - 1) // 2, p1)
        return -1 if ls == p1 - 1 else ls

    """ Find a quadratic residue (mod p) of 'a'. p
        must be an odd prime.
        Solve the congruence of the form:
            x^2 = a (mod p)
        And returns x. Note that p - x is also a root.
        0 is returned is no square root exists for
        these a and p.
        The Tonelli-Shanks algorithm is used (except
        for some simple cases in which the solution
        is known from an identity). This algorithm
        runs in polynomial time (unless the
        generalized Riemann hypothesis is false).
    """
    # Simple cases

    if legendre_symbol(a, p) != 1:

        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return p
    elif p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Partition p-1 to s * 2^e for an odd s (i.e.
    # reduce all the powers of 2 from p-1)
    #
    s = p - 1
    e = 0
    while s % 2 == 0:
        s //= 2
        e += 1

    # Find some 'n' with a legendre symbol n|p = -1.
    # Shouldn't take long.
    #
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1

    # Here be dragons!
    # Read the paper "Square roots from 1; 24, 51,
    # 10 to Dan Shanks" by Ezra Brown for more
    # information
    #

    # x is a guess of the square root that gets better
    # with each iteration.
    # b is the "fudge factor" - by how much we're off
    # with the guess. The invariant x^2 = ab (mod p)
    # is maintained throughout the loop.
    # g is used for successive powers of n to update
    # both a and b
    # r is the exponent - decreases with each update
    #
    x = pow(a, (s + 1) // 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m
