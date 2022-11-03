from PolynomialsFiniteFields import*
from FiniteFieldZP import *
from QuotientRing import *
from Euclidean import*
#from Prototyping import*
from ClassicalModularPolynomials import*
from PointPolynomialArithmetic import*
from CanonicalModularPolynomials import*
from Frac import*
import scipy.special
import numpy as np
import time
import math
import itertools
import random as rand
import timeit

Primes = get_primes(1000)

first = get_primes(10000)[5:]


def pointfactory(a,b,p,field):
    
 class CurvePoint(object):
     def __init__(self,*args):
        
           if not args:
            while True:
             x = rand.randint(1, p-1)
             m = (x**3+a*x+b)%p
             if m == 0:
                self.x, self.y, self.z = field(x), field.zero,field.one
                break
             y = modular_sqrt(m,p)
             if y!= 0:
                self.x,self.y,self.z=field(x),field(y),field.one
                break
           else:
              self.x,self.y,self.z = args

     def Double(self):
           
            if self.z ==field.zero:
                return self
            if self.y == field.zero:
                return CurvePoint(field.zero, field.one, field.zero)
            else:
                m = (field(3)*(self.x**2)+field(a))/(field(2)*self.y)
                X = m**2-self.x-self.x
                Y = -m*(X-self.x)-self.y
                return CurvePoint(X,Y,field.one)
            
            
     def Add(self,Q):
         
            if self.z == field.zero:
                return Q
            if Q.z == field.zero:
                return self
            if self.x ==Q.x and self.y == Q.y:
                return self.Double()
            if self.x==Q.x:
                return CurvePoint(field.zero, field.one, field.zero)
            else:
                m = (self.y-Q.y)/(self.x-Q.x)
                x = m**2-Q.x-self.x
                y = m*(Q.x-x)-Q.y
                return CurvePoint(x,y,field.one)
        
     def scal(self,m):
            P = self
            n = abs(m)
            res = CurvePoint(field.zero, field.one, field.zero)
            while n>0:
                 if n%2 == 1:
                     res = res.Add(P)
                 P = P.Double()
             
                 n = n//2
            if m<0 and res.z != field.zero:
                return CurvePoint(res.x, -res.y, field.one)
            return res
        
     def sub(self,Q):
            if Q.z == field.zero:
                return self
            R = CurvePoint(Q.x, -Q.y, field.one)
            return self.Add(R) 
        
        
        
 return CurvePoint


class EllipticCurve():
        def __init__(self, p, a,b):
            self.field = IntegersModP(p)
            self.P = PolynomialsOver(self.field)
            self.a = a
            self.b = b
            self.EC = self.P([self.b,self.a,0,1])
            self.p = p
            if self.j_inv()==self.field.zero:
                print("Singular Curve")
                
        def j_inv(self):
            '''compute j-invariant of Elliptic Curve'''
            A = self.field(self.a)
            B = self.field(self.b)
            j = (self.field(1728)*self.field(4)*(A**3))/(self.field(4)*(A**3)+self.field(27)*(B**2))
            return j
                
            
        def DivisionPoly(self,n):
            A = self.a
            B = self.b
            y2 = self.P([B,A,0,1])
            y4 = y2*y2
            phi0 = self.P([0])
            phi1 = self.P([1])
            phi2 = self.P([2])
            phi3 = self.P([-(A**2), 12*B, 6*A, 0, 3])
            phi4 = self.P([-4*(A**3)-32*(B**2), -16*A*B, -20*(A**2), 80*B, 20*A,0, 4])
            Phi = [phi0, phi1, phi2, phi3, phi4]
            
            if n <5:
                return Phi[0:n+1]
            i = 5
            while i <= n:
             if i%4 == 1:
                m = int((i-1)/2)
                
                new = (y4*(Phi[m+2]*(Phi[m]**3)))-(Phi[m-1]*(Phi[m+1]**3))
                Phi.append(new)
                
             if i%4 == 3:
                m = int((i-1)/2)
                new = (Phi[m+2]*(Phi[m]**3))-(y4*(Phi[m-1]*(Phi[m+1]**3)))
                Phi.append(new)
                
             if i%2 == 0:
                m = i//2
                r = self.P([(self.field(1)/self.field(2))])*Phi[m]*((Phi[m+2]*(Phi[m-1]**2))-(Phi[m-2]*(Phi[m+1]**2)))
                Phi.append(r)
                
             i = i +1
            return Phi
        
    
    
        def odd_l_subroutine(self, Phi, l):
           '''compute phi^2(P)-q(P) mod l, psi(l) in polynomial eqn
           and find lambda s.t phi^2(P)-q(P) = lambda(P)'''
           phi_l = Phi[l]
           Quotient = GenQuotientRing(phi_l, self.P)
           Ltor = GenFracField(Quotient)
           Point =  GenBasePoint(self.EC, self.a, self.b)
           P = Point(Ltor.x,Ltor.one)
           q = self.p
           ql = q%l
           if ql>l/2:
               ql = ql-l

           Q1,Q2 = P.Frobenius(q**2),P.mulv2(ql)
           if Q1.x == Q2.x:
               w = modular_sqrt(ql, l)
               if w == 0:
                   return 0
               else:
                   Pw = P.mulv2(w)
                   if (P.x**q-Pw.x).num.polynomial.gcd(phi_l).deg()==0:
                       return 0
                   else:
                       if (P.y**q-Pw.y).num.polynomial.gcd(phi_l).deg() == 0:
                           return -2*w%l
                       else:
                           return 2*w%l
               
           Q = Q1+Q2        
           x_q, y_q = Q.x, Q.y
      
           bound = int((l-1)/2)
           j = 1
           jP = P
           x_j, y_j = jP.x, jP.y
    
           while j <= bound:   
               if (x_j**q-x_q).num.polynomial.gcd(phi_l).deg()>0:
                   break
               j = j+1
               jP = jP+P
               x_j, y_j = jP.x, jP.y 
           
           yj = (y_j**q)*((Ltor(Quotient(self.EC))**((q-1)//2)))
           
           if (yj-y_q).num.polynomial.gcd(phi_l).deg()>0:
               return j
           else:
               return -j
           
              
        def a_mod2(self):
            q = self.p
            Q = GenQuotientRing(self.EC, self.P)
            GCD= (Q.x**q-Q.x).polynomial.gcd(self.EC)
            if  GCD.deg() == 0:
             
                return 1
            return 0
           
    
  
        def Schoof(self):
            t1 = time.time()
            primeslist = primes(self.p)
            print(primeslist)
            t2 = time.time()
            #print('get primes:'+ str(t2-t1))
            L = primeslist[-1]
            Phi = self.DivisionPoly(L)
            #print(Phi)
            c1 = []
            c2 = []
            t1 = time.time()
            x,y = self.a_mod2(), 2
            #print(x,y)
            t2 = time.time()
            #print('check a mod 2: '+str(t2-t1))
            c1.append(x)
            c2.append(2)
            
            for l in primeslist[1:]:
                 t1 = time.time()
                 x =  self.odd_l_subroutine(Phi, l)
                 print(x)
                 t2 = time.time()
                 #print(str(l)+': '+str(t2-t1))
                 c1.append(x)
                 c2.append(l)
            #print(c1)
            a, M = CRT(c1,c2)
           
            if abs(a)>2*math.sqrt(self.p):
                a = a-M
            return self.p+1-a
            
        
        def Eisensteins(self):
            E_4 = self.field(-48)*self.field(self.a)
            E_6 = self.field(864)*self.field(self.b)
            return E_4, E_6
        
        
            
        
        
        def j_inv_prime(self):
            E_4,E_6 = self.Eisensteins()
            jprime = (-E_6/E_4)*self.j_inv()
            return jprime
        
        def PrimeCheck(self, l, ModularPolynomial):
            Q = GenQuotientRing(ModularPolynomial, self.P)
            x = Q(self.P.x)
            Xq =(x**self.p).polynomial
            FieldPoly = Xq-(x.polynomial)
            GCD = ModularPolynomial.gcd(FieldPoly)
            lead = GCD.coefficients[-1]
            GCD = GCD*self.P([self.field(1)/lead])
            return GCD
        
        def IsogenousCurve(self,l, ModularPolynomial, gcd):
            #print(gcd.deg())
            j = self.j_inv()
            if gcd.deg()==1:
             jtilde= -gcd.coefficients[0]
            if gcd.deg() == 2:
             Delta = gcd.coefficients[1]**2-self.field(4)*gcd.coefficients[0]*gcd.coefficients[2]
             disc = self.field(modular_sqrt(Delta.n,self.p))
             jtilde= (-gcd.coefficients[1]-disc)/(self.field(2)*gcd.coefficients[2])
   
                
            j = self.j_inv()
            jprime = self.j_inv_prime()
            Phi_x = ModularPolynomial.partial_x(j, jtilde)
            Phi_y = ModularPolynomial.partial_y(j, jtilde)
            jtilde_prime = -(jprime*Phi_x)/(self.field(l)*Phi_y)
            atilde = -(jtilde_prime**2)/(self.field(48)*jtilde*(jtilde-self.field(1728)))
            btilde = -(jtilde_prime**3)/(self.field(864)*(jtilde**2)*(jtilde-self.field(1728)))
            return j,jprime, jtilde, jtilde_prime, atilde, btilde
            
        def ElkiesProcedure(self, l, Phi_l, gcd):
             j,jprime, jtilde, jtilde_prime, atilde, btilde = self.IsogenousCurve(l,Phi_l,gcd)
             
             t1 = jprime**2*Phi_l.partial_xx(j, jtilde)
             t2 = self.field(2*l)*jprime*jtilde_prime*Phi_l.partial_xy(j, jtilde)
             t3 = self.field(l**2)*(jtilde_prime**2)*Phi_l.partial_yy(j, jtilde)
         
            
             E4_ql = -self.field(48)*atilde
             E6_ql = self.field(864)*btilde
             E4, E6 = self.Eisensteins()
             vii14 = -(t1+t2+t3)/(jprime*Phi_l.partial_x(j, jtilde))
             
             T1 = self.field(l)/self.field(2)*vii14
             T2 = self.field(l)/self.field(4)*((E4)**2/E6-self.field(l)*((E4_ql)**2/E6_ql))
             T3 = self.field(l)/self.field(3)*(E6/E4-self.field(l)*(E6_ql/E4_ql))
             '''coefficient of x^((l-1)/2-1)'''
             p1 = T1+T2+T3
            
             atilde = self.field(l**4)*atilde
             btilde = self.field(l**6)*btilde
             d = int((l-1)/2)
             
             C_k, C_ktilde = self.WeierstrassCoeffecients(d, self.field(self.a), self.field(self.b), atilde, btilde)
             
             coeffs = self.KernelPolyCoefficients(p1,d,C_k,C_ktilde)
             return self.P(coeffs)
             
       
        def ElkiesProcedureV2(self, l, Phi,gcd):
            E_4 = -self.field(self.a)/self.field(3)
            E_6 = -self.field(self.b)/self.field(2)
            Delta = (E_4**3-E_6**2)/self.field(1728)
            j = self.j_inv()
            if gcd.deg()==1:
             print('deg = 1')
             g= -gcd.coefficients[0]
            if gcd.deg() == 2:
             print('2 roots')
             disc = gcd.coefficients[1]**2-self.field(4)*gcd.coefficients[0]*gcd.coefficients[2]
             disc = self.field(modular_sqrt(disc.n,self.p))
             g= (-gcd.coefficients[1]-disc)/(self.field(2)*gcd.coefficients[2])
           # if gcd.deg()==l+1:
               # print('l+1 roots')
              #  g = gcd.findroot()
            Dg = g*Phi.partial_x(g,j)
            Dj = j*Phi.partial_y(g,j)
            s = int(12/math.gcd(l-1,12))
            Delta_l = (g**(int(12/s)))*Delta/(self.field(l)**12)
            if Dj == self.field(0):
                print("Dj==0")
                E4_l = self.E_4/(self.field(l)**2)
                atilde = self.field(-3*(l**4))*E4_l
                jtilde = E4**3/Delta_l
                btilde = self.field(2*(l**6))*self.field(modular_sqrt((((jtilde-self.field(1728))*Delta_l)).n, self.p))
                p1 = self.field(0)
            else:
                E_2 =-self.field(12)*E_6*Dj/(self.field(s)*E_4*Dg)
                gprime = -(self.field(s)/self.field(12))*E_2*g
                jprime = -(E_4**2)*E_6/Delta
                E_0 = E_6/(E_4*E_2)
                Dgprime = gprime*Phi.partial_x(g,j)+g*(gprime*Phi.partial_xx(g,j)+jprime*Phi.partial_xy(g,j))
                Djprime = jprime*Phi.partial_y(g,j)+j*(jprime*Phi.partial_yy(g,j)+gprime*Phi.partial_xy(g,j))
                E0_prime = (((self.field(-s)/self.field(12))*Dgprime)-E_0*Djprime)/Dj
                E4_l = (self.field(1)/self.field(l**2))*(E_4-E_2*((self.field(12)*E0_prime/E_0) +(self.field(6)*(E_4**2/E_6))-(self.field(4)*E_6/E_4))+E_2**2)
                jtilde = E4_l**3/Delta_l
                f = (self.field(l)**s)/g
                fprime = self.field(s)*E_2*f/self.field(12)
                Dgstar = Phi.partial_x(f,jtilde )
                Djstar = Phi.partial_y(f, jtilde)
                jtilde_prime = -fprime*Dgstar/(self.field(l)*Djstar)
                E6_l = -E4_l*jtilde_prime/jtilde
                atilde = self.field(-3*(l**4))*E4_l
                btilde = self.field(-2*(l**6))*E6_l
                p1 = self.field(-l)*E_2
            d = int((l-1)/2)
            Ck,Ck_tilde = self.WeierstrassCoeffecients(d, self.field(self.a), self.field(self.b), atilde, btilde)
            coeffs = self.KernelPolyCoefficients(p1,d,Ck,Ck_tilde)
            return self.P(coeffs)
                
        def WeierstrassCoeffecients(self,d,a,b,atilde, btilde):
            '''compute coefficients of weierstrass p function for E and E-tilde'''
            Ck = [-a/self.field(5), -b/self.field(7)]
            Ck_tilde = [-atilde/self.field(5), -btilde/self.field(7)]
            for k in range(3,d+1):
                coeffk = self.field(3)/(self.field((k-2)*(2*k+3)))
                c_k = coeffk*sum(Ck[j-1]*Ck[k-2-j] for j in range(1,k-1))
                c_ktilde = coeffk*sum(Ck_tilde[j-1]*Ck_tilde[k-2-j] for j in range(1,k-1))
                Ck.append(c_k)
                Ck_tilde.append(c_ktilde)
            return Ck, Ck_tilde
        
   
        def KernelPolyCoefficients(self,p1, d, Ck, Ck_tilde):
            '''Compute Coefficients for factor of division polynomial that 
            vanhishes on eigenspace of the frobenius endomorphism'''
            Ak = [-p1/self.field(2)]
            for j in range(1,d):
                 Ak.append(-(Ck_tilde[j-1]-self.field(2*d+1)*Ck[j-1])/self.field((2*j+1)*(2*j+2)))
            C_series = [self.P.one]
            
            f = self.P(Ak)
            A = self.P([self.field(1)]+Ak)
            s = self.P(Ck)
            C = self.P([self.field.zero]+Ck)

            C_series.append(C)
            j = 2
            g = f
            h = self.P(Ck)
            while j <= d:
                g = g.mul_upto(f,d-j)
                q = self.P([self.field(1)/self.field(math.factorial(j))])*g.shift(j)
                A = A+q
                h = h.mul_upto(s, d-j)
                C_series.append(h.shift(j))
                j = j+1
            Coeffs = [self.field.zero]*(d-1)+[-p1/self.field(2)]+[self.field.one]
            for i in range(1,d+1):
                fl = A.coefficients[i]
                for k in range(1,i+1):
                    c = sum(self.field(math.factorial(d-i+k)/(math.factorial(k-j)*math.factorial(d-i+j)))*C_series[k-j].coefficients[j] for j in range(0, k))
                    fl = fl-c*Coeffs[d-i+k]
                Coeffs[d-i] = fl
            return Coeffs
            
       
        def FrobeniusOrder(self,l,Phi):
            '''compute the order of frobenius r in the case of an atkin prime l so that 
            the lth modular polynomial splits over F_p^r'''
            b = modular_sqrt(l,self.p)
            if b == 0:
                b = -1
            else:
                b = 1
            Q = GenQuotientRing(Phi,self.P)
            k=0
            h = Q.x
            for i in range(1,l+2):
                if (l+1)%i == 0:
                    for c in range(1, i-k+1):
                     h = h**self.p
                    k=i
                   #+ s = h.polynomial-Q.x.polynomial
                   # print(s)
                    g = (h.polynomial-Q.x.polynomial).gcd(Phi)
                    
                   # print(g)
                    if g.deg()==l+1:
                      return i
            
        def Atkin(self, l,Phi,r):
           '''Process an Atkin prime'''
           Fl = IntegersModP(l)
           Pl = PolynomialsOver(Fl)
           d = get_non_square(l)
           irreducible_modulus = Pl([-d,0,1])
           Ql = GenQuotientRing(irreducible_modulus, Pl)
           divs = list(factors(l**2-1))
           divs.sort()
           divs = divs[1:-1]
           G = Ql.one
           while True:
               #find a generator of multiplicative group F_l^2
               a = np.random.randint(0,l)
               g = Ql(Pl([a,1])) #sqrt{d}-a
               t = 1
               for d in divs:
                   if g**((l**2-1)/d)== Ql.one:
                       t = 0
                       break
               if t == 1:
                   G = g
                   break
               
           #compute candidates for rth roots of unity in F_l
           g = G**((l**2-1)/r)
           s = Ql.one
           unities = []
           k=0
           relprimes = get_relprimes(r)
           for i in relprimes:
               s = s*(g**(i-k))
               k=i
               unities.append(s)
           Traces = []
           
           for gamma_r in unities:
               
               g1 = gamma_r.polynomial.coefficients[0].n
               x1_sq = (Fl(self.p)*Fl(g1+1)/Fl(2))
               if x1_sq  == Fl.zero:
                   Traces.append(0)
               x1 = modular_sqrt(x1_sq.n,l)
           
               if x1 != 0:
                   Traces.append(2*x1%l)
                   Traces.append(-2*x1%l)
              
                   
           return list(set(Traces))
   
           
        def AtkinTraces(self,Atkin,number_of_traces):
            L = len(Atkin[0])
            number_of_atkin = L
            i = 0
            traces_count = 0
            S1 = []
            S1_traces = []
            if L == 1:
                S1 = Atkin[1]
                S1_traces = Atkin[0]
                T1 = [(trace, S1[0]) for trace in S1_traces]
                S2 = []
                S2_traces = []
                return T1, []
           
                
            while i <L:
                if traces_count +len(Atkin[0][i]) > number_of_traces/2:
                    S2 = Atkin[1][i:]
                    S2_traces = Atkin[0][i:]
                    break
                traces_count += len(Atkin[0][i])
                S1.append(Atkin[1][i])
                S1_traces.append(Atkin[0][i])
                i = i+1
            
            S1_traces_product = itertools.product(*S1_traces)
            S2_traces_product = itertools.product(*S2_traces)
           
            
            T1 = [CRT(trace, S1) for trace in S1_traces_product]
            T2  = [CRT(trace, S2) for trace in S2_traces_product]
           
            return T1,T2
        
        def BSGS(self, Elkies, Atkin1, Atkin2):
            t3,m3 = CRT(Elkies[0], Elkies[1])
            PointFactory = pointfactory(self.a, self.b, self.p, self.field)
           
            T1, T2 = Atkin1, Atkin2
            if T2 == []:
                m1 = T1[0][1]
               # print(T1)
                P = PointFactory()
              
                for t1 in T1[0][0]:
                    t,_ = CRT([t1,t3], [m1,m3])
                    if abs(t)>2*math.sqrt(self.p):
                        t = t-_
                   # print(P.scal(self.p+1-t))
                    if P.scal(self.p+1-t).z == self.field.zero:
                        return t
                            
            m1, m2 = T1[0][1],T2[0][1]
           # print(m2,m1)
            R1, R2 = [], []
            R1, R2 = sorted(R1, key = abs), sorted(R2, key = abs)
        
            
            for t1,m1 in T1:
             r1 = ((t1-t3)*extendedEuclideanAlgorithm(m2*m3, m1)[0])%m1
             if r1>m1/2:
                 R1.append(r1-m1)
             else:
                 R1.append(r1)
            for t2, m2 in T2:
             r2 = ((t2-t3)*extendedEuclideanAlgorithm(m1*m3, m2)[0])%m2
             R2.append(r2)
             R2.append(r2-m2)
             i =1
            while True:
             print("loop iters"+str(i))
             P = PointFactory()
             #print(P.x, P.y)
             #print("P",P)
             Q = P.scal(self.p+1-t3)
             m3P = P.scal(m3)
             m1m3P = m3P.scal(m1)
             m2m3P = m3P.scal(m2) 
             Q_r1= {}
          
             for r1_i in R1:
                Q_r1_i = Q.sub(m2m3P.scal(r1_i))
                if Q_r1_i.z != self.field.zero:
                 Q_r1[Q_r1_i.x.n] = Q_r1_i.y.n, r1_i
                if Q_r1_i.z == self.field.zero:
                    Q_r1["O"] = r1_i
            
             for r2_i in R2:
                 r2_im1m3P = m1m3P.scal(r2_i)
                 if r2_im1m3P.z != self.field.zero:
                  if r2_im1m3P.x.n in Q_r1:
                    
                     y_r_1_i, r_1_i = Q_r1[r2_im1m3P.x.n]
                     if  y_r_1_i== r2_im1m3P.y.n:
                        # print("match")
                         t = t3+m3*(m2*r_1_i+m1*r2_i)
                         return t
                     else:
                      
                         t = t3+m3*(m2*r_1_i-m1*r2_i)
                         return t
                 if  r2_im1m3P.z == self.field.zero:
                    if "O" in Q_r1:
                        r_1_i = Q_r1["O"]
                        t = t3+m3*(m2*r_1_i+m1*r2_i)
                        return t
             i+=1
           
                       
            
            
            
                
                
        
       
            
        def SEA(self):

            primeslist = primes(self.p)
            x,y = self.a_mod2(), 2
            prod = 2
            E1 = [x]
            E2 = [y]
            
            A1 = []
            A2 = []
            number_of_traces = 0
            counter  =1
            while prod <=  4*math.sqrt(self.p):
                l = Primes[counter]
                Phi_l = CanonicalModularPolynomial(self.p, l)
             
                j = self.j_inv()
                Phi = Phi_l.eval_at_y(j)
                GCD = self.PrimeCheck(l, Phi)
                Fl = IntegersModP(l)
                if GCD.deg()==l+1:
                    counter += 1
                  #  print("skip")
                    continue
                
                if GCD.deg() == 0: 
                  print("Atkin")
                  if (len(set(primesfac(l+1)))<= 2 and l >20) or l>100:
                      counter +=1
                      print(str(l)+' rejected')
                      continue
                  r= self.FrobeniusOrder(l,Phi)
                  #print(r)
                  traces = self.Atkin(l,Phi,r)
                 # print(l, traces)
                  number_of_traces += len(traces)
                  A1.append(traces)
                  A2.append(l)
                  counter +=1
                  prod = l*prod
                
                if GCD.deg()>0:
                 
                 print("Elkies",l)
                 F_l = self.ElkiesProcedureV2(l, Phi_l, GCD)
        
                 Ltorring = GenQuotientRing(F_l, self.P)
                 Ltor = GenFracField(Ltorring)
                 Point = GenBasePoint(self.EC, self.a,self.b)
                 P = Point(Ltor.x, Ltor.one)
                 FrobP = P.Frobenius(self.p)
                 if GCD.deg() == 2:
                    print('deg 2')
                    print()
                    bound = int((l-1)/2)+1
                    j = 1
                    jP = P
                    x_j, y_j = jP.x, jP.y
                    light = 0
                    while j <= bound:   
                        if (x_j-FrobP.x).num.polynomial.gcd(F_l).deg()>0:
                      #   print("Hit",j)
                         light = 1
                         break
                        
                        jP = jP+P
                        j = j+1
                        x_j, y_j = jP.x, jP.y 
                    if light == 0:
                        print('no eigenvalue found')
                    if (y_j-FrobP.y).num.polynomial.gcd(F_l).deg()>0:
                      t= ((Fl(j)+(Fl(self.p)/Fl(j))).n)%l
                    else:
                     t =((Fl(-j)+(Fl(self.p)/Fl(-j))).n)%l
                  
            
                    
                 else:
                  #print('deg 1')
                  w = modular_sqrt(self.p,l)
                  wP = P.multiply(w)
                  if (wP.y-FrobP.y).num.polynomial.gcd(F_l).deg()>0:
                     t =( 2*w%l)
                  else:
                      t = ( -2*w%l)
                 E1.append(t)
                 E2.append(l)
                 counter +=1
                 print(str(l)+" : "+str(t))
                 prod = l*prod
                 
                  
                    
            Elkies = E1,E2
            if A1 == []:
                t,m = CRT(E1,E2)
                if abs(t) > 2*math.sqrt(self.p):
                    return t-m
                return t
            Atkin = A1,A2
         
            Atkin1, Atkin2 = self.AtkinTraces(Atkin, number_of_traces)
            t = self.BSGS(Elkies, Atkin1, Atkin2)
            
            return "#E(F_p) = " + str(self.p+1-t)




      