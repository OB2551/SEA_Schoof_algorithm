from Euclidean import*
from gmpy2 import*
def GenBasePoint(Curve, A,B):
  class Point():
    def __init__(self, x, y,z=1):
        self.x = x
        self.y = y
        self.z = z
      
    
    def __add__(self, other):
         if other.z == 0:
            return self
         if self.x == other.x:
            if self.y == other.y:
                return self.Double()
         return self.Add(other)

    def Frobenius(self, q):
        F = self.x.__class__
        Q = F.quotient
        x_new = self.x**q
        y_new = (F(Q(Curve))**((q-1)//2))*self.y
        return Point(x_new,y_new)
    
    def __eq__(self,other):
        if self.x == other.x and self.y == other.y:
            return True
        return False
        
    def Add(self, other):
        F = self.x.__class__
        Q = F.quotient
        mnum = (self.y-other.y)
        mden = (self.x-other.x)
        m = mnum/mden
        m2 = (m**2)*F(Q(Curve))
        x_new = m2-other.x-self.x
        y_new = -m*(x_new-other.x)-other.y
        return Point(x_new,y_new)
 
        
    def Double(self):
        F = self.x.__class__
        Q = F.quotient
        ring = Q.ring
        three = F(Q(ring([3])))
        a = F(Q(ring([A])))
        m_num = three*(self.x**2)+a
        two = F(Q(ring([2])))
        m = (m_num/(two*self.y))

        m2 = m**2/F(Q(Curve))

        x_new = m2-two*self.x
        y_new = (-m/F(Q(Curve))*(x_new-self.x))-self.y
       
        return Point(x_new, y_new)
    
    def __neg__(self):
        return Point(self.x, -self.y)
     
    def multiply(self, n):
        m = abs(n)
        P = self
        i = 1
        while i<m:
            P = P+self
            i = i+1
        if n>0:
         return P
        return -P
    
    def mulv2(self, n):
        m = abs(n)
    
        P = Point.O
        Q = self
       
        while m >0:
            if m%2 == 1:
                P = Q+P
            m = m//2
            Q = Q.Double()
        if n>=0:
            return P
        
        return Point(P.x, -P.y)
    
  Point.O = Point(0,1,0)
  return Point