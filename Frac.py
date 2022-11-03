from QuotientRing import *
def GenFracField(quotient):
    class Frac():
        def __init__(self, num, dem=None):
            self.num = num
            if dem is None:
                self.dem = quotient.one 
            else:
                self.dem = dem
            
        def __add__(self, other):
            new_num = self.num*other.dem+other.num*self.dem
            new_dem = self.dem*other.dem
            #if new_dem == quotient.zero:
                #print("error divide by 0 occured in Frac(R)")
            return Frac(new_num, new_dem)
            
        def __sub__(self, other):
            new_num = self.num*other.dem-other.num*self.dem
            new_dem = self.dem*other.dem
           # if new_dem == quotient.zero:
              #  print("error divide by 0 occured in Frac(R)")
            return Frac(new_num, new_dem)
        
        def __mul__(self, other):
            new_num = self.num*other.num
            new_dem = self.dem*other.dem
            return Frac(new_num, new_dem)
        
        def inverse(self):
            #if self == Frac.zero:
                #print("0, cant invert in Frac(R)")
            
            return Frac(self.dem, self.num)
            
        def __truediv__(self, other):
            if other != Frac.zero:
             new_num = self.num*other.dem
             new_dem = self.dem*other.num
             return Frac(new_num, new_dem)
           # print("div by 0 error in Frac(R)")
            
        def __eq__(self, other):
            if self.num*other.dem == self.dem*other.num:
                return True
            return False
        
        def __neg__(self):
            return Frac(-self.num, self.dem)
        
        def __pow__(self,n):
            return Frac(self.num**n, self.dem**n)
        
        def __repr__(self):
            return str(self.num)+"/"+str(self.dem)
        
        
    Frac.quotient = quotient        
    Frac.zero = Frac(quotient.zero, quotient.one)
    Frac.one = Frac(quotient.one, quotient.one)  
    Frac.x = Frac(quotient.x)
    return Frac