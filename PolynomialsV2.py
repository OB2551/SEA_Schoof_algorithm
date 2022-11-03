import numpy as np
import random
from Euclidean import *
try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest
import math


      
def PolynomialsOver(M):
    class Polynomial():
        def __init__(self, C):
            self.coefficients = [c%M for c in C]
            c = self.coefficients
            '''remove any leading 0s'''
            while len(c)>1 and c[-1] == 0 :
                    del c[-1]
                    

        def __add__(self, other):
         newCoefficients = [sum(x)%M for x in zip_longest(self.coefficients, other.coefficients, fillvalue=0)]
         return Polynomial(newCoefficients)     
     
        def __sub__(self, other):
            return (self + -other)
        
                
        def __mul__(self, other):
             
              result = [ 0 ] * (self.deg() + other.deg() + 2)
              for i, x in enumerate(self.coefficients):
                 for j, y in enumerate(other.coefficients):
                     result[i + j]  +=  x * y%M
              return Polynomial( result )
        
       
        def __len__(self):
            return len(self.coefficients)
        
        def deg(self):
            return len(self)-1
        
        def __neg__(self): return Polynomial([-a for a in self.coefficients])
            
        def __repr__(self):
          x = ' + '.join(['%s x^%d' % (a,i) if i > 0 else '%s'%a
                              for i,a in enumerate(self.coefficients)])
          return x

      
      
        def __rmul__(self,other):
            return Polynomial(other*self)
        
        def __radd__(self, other):
            return Polynomial(other + self)
        
        def __divmod__(self, other):
            if other.deg()>self.deg():
                return Polynomial.zero, self
            dividend = self.coefficients[:]
            divisor = other.coefficients[:]
            n = other.deg()
           
            quotient = [ 0 ] * (self.deg() - n + 1)
            for k in reversed(range( 0, len(quotient) )):

                quotient[k] = pow(divisor[n], -1, M)*dividend[n + k]
                for j in range(k, n + k):
                   dividend[j] -= quotient[k] * divisor[j - k]
            remainder = dividend[ 0 : n]
            if remainder == []:
                r = Polynomial.zero
                return Polynomial(quotient), r
            return Polynomial(quotient), Polynomial(remainder)
       
    
        def remlead(self):
            return Polynomial(self.coefficients[0:-1])
                
        def __floordiv__(self, other):
            x, _ = divmod(self, other)
            return x
        
             
        def __pow__(self, n):
          
           Q = self
           R = Polynomial.one
           while n > 0:
             if n%2 == 1:
               R = Q * R
             Q = (Q * Q)
             n = n//2

           return R
        
        def __mod__(self, other):
            _,r  = divmod(self, other)
            return r
        
        def __eq__(self, other):
            if self.deg() != other.deg():
             return False
            if self.coefficients == other.coefficients:

             return True
            return False
             
        def gcd(self, other):
            z = polynomialGCD2(self, other)
            return z
        
        def powmod(self, n, mod):
            Q = self
            R = Polynomial.one
            while n > 0:
             if n%2 == 1:
               R = (Q * R)%mod
             Q = (Q * Q)%mod
             n = n//2
            return R
        
        def evaluate(self, x):
            res = 0
            s = 1
            for n in range(0, self.deg()+1):
                
                res = res + self.coefficients[n]*(s)%M
                s = s*x%M
            return res%M
        
        def shift(self, j):
            coeffs = [0]*j+self.coefficients
            return Polynomial(coeffs)
        
        def mul_upto(self, other, d):
            n = self.deg()
            m = other.deg()
            coeffs = [0]*(d+1)
            for i in range(min(d,n)+1):
                for j in range(min(d-i,m)+1):
                    coeffs[i+j] += self.coefficients[i]*other.coefficients[j]
            return Polynomial(coeffs)
                
    Polynomial.x = Polynomial([0,1])   
    Polynomial.zero = Polynomial([0]) 
    Polynomial.one = Polynomial([1])
    Polynomial.modulus = M
    return Polynomial       

