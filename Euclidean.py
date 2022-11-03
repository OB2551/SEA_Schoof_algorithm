
import numpy as np
import sympy as sy
import math
from functools import reduce
import operator
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

def Gcd(a, b):
   if abs(a) < abs(b):
      return Gcd(b, a)             

   while abs(b)>0:
      _,r = divmod(a,b)
      a,b = b,r

   return a


def polynomialGCD(u,v):
    '''for polynomials'''
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v
 
    while lesser_remainder != zero:
        quotient, next_remainder = divmod( larger_remainder, lesser_remainder )
        
        
        # Shift roles
        larger_remainder = lesser_remainder
        lesser_remainder = next_remainder

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    c = larger_remainder.coefficients[0]
    C = larger_remainder.Field.one/c
    p = larger_remainder.__class__([C])
    return p*larger_remainder

def polynomialGCD2(u,v):
    '''for polynomials'''
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v
 
    while lesser_remainder != zero:
        quotient, next_remainder = divmod( larger_remainder, lesser_remainder )
        
        
        # Shift roles
        larger_remainder = lesser_remainder
        lesser_remainder = next_remainder

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    c = larger_remainder.coefficients[0]
    C = pow(c,-1, u.__class__.modulus)
    p = larger_remainder.__class__([C])
    return p*larger_remainder
    
    
def extendedEuclideanAlgorithm(d, e):
   '''for integers'''
  
   a = d
   b = e
   if abs(b) == 0:
      return (1, 0, a)

   x1, x2 = 0, 1
   while abs(b) > 0:
      q, r = divmod(a,b)

      x = x2 - q*x1
     
      a, b, x2, x1= b, r, x1, x
   c = (a-x2*d)//e
   return (x2, c, a)


def extended_euclidean_algorithm(u, v):
    '''for polynomials'''
    one = u.__class__.one
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v
    larger_scalar, lesser_scalar = one, zero
    while lesser_remainder != zero:
        quotient, next_remainder = divmod( larger_remainder, lesser_remainder )
        next_scalar = larger_scalar - quotient * lesser_scalar 
        
        # Shift roles
        larger_remainder, larger_scalar = lesser_remainder, lesser_scalar
        lesser_remainder, lesser_scalar = next_remainder, next_scalar

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    other_scalar = ( larger_remainder  -  u * larger_scalar ) //v  
    c = larger_remainder.coefficients[0]
    C = larger_remainder.Field.one/c
    p = larger_remainder.__class__([C])
    return (p*larger_scalar, p*other_scalar, p*larger_remainder )

def extended_euclidean_algorithm2(u,v):
    '''for polynomials'''
    one = u.__class__.one
    zero = u.__class__.zero
    larger_remainder, lesser_remainder = u, v
    larger_scalar, lesser_scalar = one, zero
    while lesser_remainder != zero:
        quotient, next_remainder = divmod( larger_remainder, lesser_remainder )
        next_scalar = larger_scalar - quotient * lesser_scalar 
        
        # Shift roles
        larger_remainder, larger_scalar = lesser_remainder, lesser_scalar
        lesser_remainder, lesser_scalar = next_remainder, next_scalar

    # larger_remainder is the gcd;
    # larger_scalar is the factor for the linear combination
    other_scalar = ( larger_remainder  -  u * larger_scalar ) //v  
    c = larger_remainder.coefficients[0]
    p = u.__class__.modulus
    C = pow(c, -1,p)
    p = larger_remainder.__class__([C])
    return (p*larger_scalar, p*other_scalar, p*larger_remainder )

def Fastpow(a, n):
      Q = a
      R = 1
      while n > 0:
         if n%2 == 1:
            R = Q * R
         Q = (Q * Q)
         n = n//2
    

      return R

def powermod(p,n,mod,char):
    x = sy.Symbol('x')
    Q = p
    R = sy.Poly(1, x, modulus = char)
    while n>0:
        if n%2 == 1:
            R = ((Q*R)%mod)
        Q = ((Q*Q)%mod)
        n = n//2
    return R

def CRT(congruences1, congruences2):
    x = 0
    M = prod(congruences2)
 
   
    for i in range(len(congruences1)):
        n_i = congruences1[i]
        m_i = congruences2[i]
        
        x = x + (M//m_i)*pow(M//m_i,-1, m_i)*n_i
        
    return int(x%M), M
    
def get_primes(n):
    numbers = set(range(n, 1, -1))
    primes = []
    while numbers:
        p = numbers.pop()
        primes.append(p)
        numbers.difference_update(set(range(p*2, n+1, p)))
    return primes

def primes(p):
    primes = get_primes(1000)
    prod = 1
    i = 0
    j = -1
    while prod < 4*math.sqrt(p):
        if primes[i] != p:
         prod = prod*primes[i]
        else:
            j = i
        i += 1
    if j != -1:
     return primes[0:j]+primes[j+1:i]
    return primes[0:i]

def get_non_square(p):
     while True:
                a = np.random.randint(2, p)
                if modular_sqrt(a,p) ==0:
                    return a
                
         
def phi(n):
    amount = 0        
    for k in range(1, n + 1):
        if math.gcd(n, k) == 1:
            amount += 1
    return amount                
                
def get_relprimes(n):
    s = []
    for i in range(1,n):
        if math.gcd(i,n)==1:
            s.append(i)
    return s

def factors(n):    
    return set(reduce(list.__add__, 
                ([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))             


def primesfac(n):
    primfac = []
    d = 2
    while d*d <= n:
        while (n % d) == 0:
            primfac.append(d)  # supposing you want multiple factors repeated
            n //= d
        d += 1
    if n > 1:
       primfac.append(n)
    return primfac                

def modular_sqrt(a, p):

    def legendre_symbol(a, p):
        """ Compute the Legendre symbol a|p using
            Euler's criterion. p is a prime, a is
            relatively prime to p (if p divides
            a, then a|p = 0)
            Returns 1 if a has a square root modulo
            p, -1 otherwise.
        """
    
        ls = pow(a, (p - 1) // 2, p)
        return -1 if ls == p - 1 else ls

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


    