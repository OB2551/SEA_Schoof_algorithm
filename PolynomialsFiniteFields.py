import numpy as np
import FiniteFieldZP
import random
from Euclidean import *
try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest

def PolynomialsOver(Field):
    class Polynomial():
        def __init__(self, C):
            '''check if C same type'''
            if type(C) is Polynomial:
                self.coefficients = C.coefficients
                
            if type(C[0]) is Field:
                self.coefficients = C
               
            else:
                self.coefficients = [Field(c) for c in C] 
            c = self.coefficients
            '''remove any leadinf 0s'''
            while len(c)>1 and c[-1] == Field.zero :
                    del c[-1]
                    
            
          
                
        def __add__(self, other):
         newCoefficients = [sum(x) for x in zip_longest(self.coefficients, other.coefficients, fillvalue=self.Field.zero)]
         return Polynomial(newCoefficients)     
     
        def __sub__(self, other):
            return (self + -other)
        
                
        def __mul__(self, other):
              zero = Field.zero
              result = [ zero ] * (self.deg() + other.deg() + 2)
              for i, x in enumerate(self.coefficients):
                 for j, y in enumerate(other.coefficients):
                     result[i + j]  +=  x * y 
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
            #if other.deg()==0:
              #  return self,self
            if other.deg()>self.deg():
                return Polynomial.zero, self
            dividend = self.coefficients[:]
            divisor = other.coefficients[:]
            n = other.deg()
            zero = Field.zero
            quotient = [ zero ] * (self.deg() - n + 1)
            for k in reversed(range( 0, len(quotient) )):
                quotient[k] = dividend[n + k] / divisor[n]
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
            z = polynomialGCD(self, other)
            return z
        
        def powmod(self,n,mod):
            Q = self
            R = Polynomial.one
            while n > 0:
             if n%2 == 1:
               R = (Q * R)%mod
             Q = (Q * Q)%mod
             n = n//2
            return R
        
        def evaluate(self, x):
            res = self.Field.zero
            s = self.Field.one
            x = self.Field(x)
            for n in range(0, self.deg()+1):
                
                res = res + self.coefficients[n]*(s)
                s = s*x
            return res
        
        def shift(self, j):
            coeffs = [self.Field.zero]*j+self.coefficients
            return Polynomial(coeffs)
        
        def mul_upto(self, other, d):
            n = self.deg()
            m = other.deg()
            coeffs = [self.Field.zero]*(d+1)
            for i in range(min(d,n)+1):
                for j in range(min(d-i,m)+1):
                    coeffs[i+j] += self.coefficients[i]*other.coefficients[j]
            return Polynomial(coeffs)
        
        def findroot(self):
            '''randomized search to find a root of f(x) given that f splits of our field F'''
            G = self
            while G.deg()>1:
                a = random.randrange(self.Field.p)
                G = G.gcd((Polynomial([a, 1]))**(int((self.Field.p-1)/2))-Polynomial.one)
        
            return -G.coefficients[0]/G.coefficients[1]
                
    Polynomial.x = Polynomial([0,1])   
    Polynomial.zero = Polynomial([0]) 
    Polynomial.one = Polynomial([1])
    Polynomial.Field = Field
    return Polynomial       
        
