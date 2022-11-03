
from Euclidean import *
from math import gcd as bltin_gcd
import numpy as np


def IntegersModP(p):
    
    class IntegerModP(int):
        def __init__(self, n):
            self.n = int(n)%IntegerModP.p
           
       
        
        
        def __add__(self, other):
             x = IntegerModP(self.n+other.n)
             return x
         
        def __sub__(self,other):
            return IntegerModP(self.n-other.n)
       
        def __truediv__(self, other):
            return IntegerModP(self.n*other.inverse())
        
        def __mul__(self, other):
            return IntegerModP(self.n*other.n)
        
        def __abs__(self):
            return abs(self.n)
        
        def __repr__(self):
            return str(self.n)
        
        def __str__(self):
           return str(self.n)
        
        def __pow__(self, power):
            return IntegerModP(pow(self.n, power,p))
        
        def inverse(self):
            return IntegerModP(pow(self.n, p-2 ,p))
        
        def __radd__(self, other ):
            return IntegerModP(other + self.n)
        
        def __neg__(self):
            return IntegerModP(-self.n)
        
        def __eq__(self, other):
            if self.n == other.n:
                return True
            return False
        
       
                
        
   
    IntegerModP.p = p
    IntegerModP.zero = IntegerModP(0)
    IntegerModP.one = IntegerModP(1)
    
    return IntegerModP