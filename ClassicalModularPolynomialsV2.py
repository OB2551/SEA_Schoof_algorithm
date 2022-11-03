#import FiniteFieldZP
import PolynomialsV2
import json
from gmpy2 import*
class ModularPolynomial():
    def __init__(self, p, l):
        Phi = open('ClassicalModularPolynomials/phi_j_'+str(l)+'.txt', 'r')
        self.coeffs = []
        self.l = l
        self.P = PolynomialsV2.PolynomialsOver(p)
        self.p = p

        for line in Phi:
            L = line.split(' ')
            pows = json.loads(L[0])
            coeff = mpz(L[1].split('\\')[0])%self.p

            pows = [int(pows[0]), int(pows[1])]
            self.coeffs.append([pows, coeff])

    
    def sub(self,x,y):
        z = 0
        for term in self.coeffs:
            z +=pow(x, term[0][0], self.p)*pow(y, term[0][1], self.p)*term[1]%self.p
            if term[0][1]>term[0][0]:
             z+= pow(y, term[0][0],self.p )*pow(x, term[0][1], self.p)*term[1]%self.p
        return z%self.p
    
    def eval_at_y(self, k):
        coefflist = [0]*(self.l+2)
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            coefflist[i]+= c*pow(k, j, self.p)
            if i>j:
               coefflist[j]+= c*pow(k,i, self.p)
        return self.P(coefflist)
    
    
    def partial_x(self,a,b):
    
        res = 0
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>0:
                res= res+i*c*(pow(a, (i-1), self.p)*(pow(b,j, self.p)))%self.p
            if i>j and j>0:
                res= res + j*c*(pow(a, (j-1),self.p)*pow(b, i, self.p))
            
        return res%self.p
    
    def partial_y(self, a,b):
        res = 0
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>0:
                res= res+i*c*(pow(a,j, self.p)*pow(b, (i-1), self.p))%self.p
            if i>j and j>0:
                res= res +j*c*pow(a, i, self.p)*pow(b, (j-1), self.p)%self.p
            
        return res%self.p
    
    def partial_xx(self,a,b):
        res = 0
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>1:
                res= res+(i-1)*i*c*pow(a, (i-2), self.p)*pow(b, j, self.p)%self.p
            if i>j and j>1:
                res= res + (j-1)*j*c*pow(a, (j-2), self.p)*pow(b,i, self.p)%self.p
            
        return res%self.p
        
    
    
    def partial_yy(self, a,b):
        res = 0
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>1:
                res= res+i*(i-1)*c*pow(a, j, self.p)*pow(b, (i-2), self.p)%self.p
            if i>j and j>1:
                res= res + j*(j-1)*c*pow(a, i, self.p)*pow(b, (j-2), self.p)%self.p
            
        return res%self.p
        
        
        

    def partial_xy(self, a, b):
        res = 0
        
        for term in self.coeffs:
            i = term[0][0]
            j = term[0][1]
            c = term[1]
            if i>0 and j>0:
                res= res+i*j*c*pow(a, (i-1), self.p)*pow(b, (j-1), self.p)%self.p
            if i>j and j>0:
                res= res + j*i*c*pow(a, (j-1), self.p)*pow(b, (i-1), self.p)%self.p
            
        return res%self.p
            
  #  def j_tilde_prime(self)


            
    