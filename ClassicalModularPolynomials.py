import json
import FiniteFieldZP
import PolynomialsFiniteFields


class ModularPolynomial():
    def __init__(self, p, l):
        Phi = open('ClassicalModularPolynomials/phi_j_'+str(l)+'.txt', 'r')
        self.coeffs = []
        self.l = l
        self.F = FiniteFieldZP.IntegersModP(p)
        self.P = PolynomialsFiniteFields.PolynomialsOver(self.F)
        for line in Phi:
            L = line.split(' ')
            pows = json.loads(L[0])
            coeff = int(L[1].split('\\')[0])
        
        
            coeff = self.F(coeff)
            #print(pows, coeff)
            pows = [int(pows[0]), int(pows[1])]
            self.coeffs.append([pows, coeff])
       
    
    def sub(self,x,y):
        z = self.F(0)
        x = self.F(x)
        y = self.F(y)
        for term in self.coeffs:
        
            
            z +=(x**term[0][0])*(y**term[0][1])*term[1]
            if term[0][1]>term[0][0]:
             z+= (y**term[0][0]*(x**term[0][1])*term[1])
            
        return z
    
    def eval_at_y(self, k):
        k = self.F(k)
        coefflist = [self.F(0)]*(self.l+2)
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            coefflist[i]+= c*(k**j)
            if i>j:
               coefflist[j]+= c*(k**i)
        return self.P(coefflist)
    
    
    def partial_x(self,a,b):
    
        res = self.F(0)
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>0:
                res= res+self.F(i)*c*((a**(i-1))*(b**j))
            if i>j and j>0:
                res= res + self.F(j)*c*((a**(j-1))*(b**i))
            
        return res
    
    def partial_y(self, a,b):
        res = self.F(0)
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>0:
                res= res+self.F(i)*c*((a**(j))*(b**(i-1)))
            if i>j and j>0:
                res= res + self.F(j)*c*((a**(i))*(b**(j-1)))
            
        return res
    
    def partial_xx(self,a,b):
        res = self.F(0)
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>1:
                res= res+(self.F(i-1)*self.F(i)*c*((a**(i-2))*(b**j)))
            if i>j and j>1:
                res= res + (self.F(j-1)*self.F(j)*c*((a**(j-2))*(b**i)))
            
        return res
        
    
    
    def partial_yy(self, a,b):
        res = self.F(0)
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>1:
                res= res+(self.F(i)*self.F(i-1)*c*((a**(j))*(b**(i-2))))
            if i>j and j>1:
                res= res + (self.F(j)*self.F(j-1)*c*((a**(i))*(b**(j-2))))
            
        return res
        
        
        

    def partial_xy(self, a, b):
        res = self.F(0)
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>0 and j>0:
                res= res+(self.F(i)*self.F(j)*c*((a**(i-1))*(b**(j-1))))
            if i>j and j>0:
                res= res + (self.F(j)*self.F(i)*c*((a**(j-1))*(b**(i-1))))
            
        return res
            
  #  def j_tilde_prime(self)


            
    