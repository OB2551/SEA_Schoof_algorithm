import FiniteFieldZP
import PolynomialsFiniteFields
from ClassicalModularPolynomials import*
from QuotientRing import *
from Euclidean import*
import math
import scipy.special

class CanonicalModularPolynomial():
    def __init__(self,p,l):
        Phi = open("CanonicalModularPolynomials/Phi"+str(l)+".txt", "r")
        self.coeffs = []
        self.l =l
        self.F = FiniteFieldZP.IntegersModP(p)
        self.P = PolynomialsFiniteFields.PolynomialsOver(self.F)
        self.p = p
        for line in Phi:
            inf = line.split()
            powx = int(inf[0])
            powy = int(inf[1])
            coeff = self.F(int(inf[2]))
            self.coeffs.append([powx, powy, coeff])
            
            
    def evaluate(self,x,y):
        z= self.F(0)
        x =self.F(x)
        y = self.F(y)
        for monomial in self.coeffs:
            degx,degy,c = monomial
            z += (x**degx)*(y**degy)*c
        return z
    
    
    def eval_at_y(self,y):
        y = self.F(y)
        coeffslist = [self.F(0)]*(self.l+2)
        for monomial in self.coeffs:
            degx,degy,c = monomial
            coeffslist[degx] += c*(y**degy)
        return self.P(coeffslist)
    
    def eval_at_x(self,x):
        x = self.F(x)
        coeffslist = [self.F(0)]*(self.l+2)
        for monomial in self.coeffs:
            degx,degy,c = monomial
            coeffslist[degy] += c*(x**degx)
        return self.P(coeffslist)
            
    def partial_x(self,x,y):
        x,y = self.F(x), self.F(y)
        val = self.F(0)
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degx>0:
                val = val+ (c*self.F(degx)*(x**(degx-1))*(y**degy))
        return val
    
    def partial_y(self,x,y):
        x,y = self.F(x), self.F(y)
        val = self.F(0)
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degy>0:
                val = val+ (c*self.F(degy)*(y**(degy-1))*(x**degx))
        return val
    
    def partial_xx(self,x,y):
        x,y = self.F(x), self.F(y)
        val = self.F(0)
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degx>1:
                val = val+ (c*self.F(degx)*self.F(degx-1)*(x**(degx-2))*(y**degy))
        return val
        
        
    def partial_xy(self,x,y):
        x,y = self.F(x), self.F(y)
        val = self.F(0)
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degx>0 and degy>0:
                val = val+ (c*self.F(degx)*self.F(degy)*(y**(degy-1))*(x**(degx-1)))
        return val
    
    
    def partial_yy(self,x,y):
        x,y = self.F(x), self.F(y)
        val = self.F(0)
        for monomial in self.coeffs:
            degx, degy, c = monomial
            if degy>1:
                val = val+ (c*self.F(degy)*self.F(degy-1)*(y**(degy-2))*(x**degx))
        return val
    
    
    
    


