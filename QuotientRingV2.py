from Euclidean import*
def GenQuotientRing(modulus, ring):
   class QuotientRing():
       def __init__(self, polynomial):
           self.polynomial = polynomial%modulus
           
           
       def __add__(self, other):
           return QuotientRing(self.polynomial+other.polynomial)
       
       def __sub__(self, other):
           return QuotientRing(self.polynomial - other.polynomial)
       
       def __mul__(self, other):
           return QuotientRing(self.polynomial*other.polynomial)
       
       def __truediv__(self, other):
 
             
             return QuotientRing(self.polynomial*other.inverse())
            
         
       def inverse(self):
            x,y,z = extended_euclidean_algorithm2(self.polynomial, self.modulus)
          
            if z.deg() == 0:
             return x
            else:
             return self.ring.zero
        
       def __neg__(self):
           return QuotientRing(-self.polynomial)
       def __repr__(self):
           return str(self.polynomial)
       
       def __pow__(self, n):
           Q = self
           R = QuotientRing.one
           while n > 0:
             if n%2 == 1:
               R = Q * R
             Q = (Q * Q)
             n = n//2

           return R
       
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