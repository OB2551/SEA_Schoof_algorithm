from PolynomialsV2 import*
#from FiniteFieldZP import *
from QuotientRingV2 import *
from Euclidean import*
from ClassicalModularPolynomialsV2 import*
from PointPolynomialArithmeticV2 import*
from CanonicalModularPolynomialsV2 import*



from sympy.ntheory import factorint
from FracV2 import*
import scipy.special
import numpy as np
import time
import math
import itertools
import random as rand
import timeit
import gmpy2
Primes = get_primes(1000)
from gmpy2 import*

def pointfactory(a,b,p):
    
 class CurvePoint(object):
     def __init__(self,*args):
        
           if not args:
            while True:
             x = rand.randint(1, p-1)
             m = (pow(x,3,p)+a*x+b)%p

             if m == 0:
                self.x, self.y, self.z = x, 0,1
                break
             y = modular_sqrt(m,p)

             if y!= 0:
                self.x,self.y,self.z=x,y,1
                break
           else:
              self.x,self.y,self.z = args

     def Double(self):
           
            if self.z ==0:
                return self
            if self.y == 0:
                return CurvePoint(0, 1, 0)
            else:
                m = ((3*(self.x**2)+a)*pow((2*self.y), p-2, p))%p
                X = (m**2-self.x-self.x)%p
                Y = (-m*(X-self.x)-self.y)%p
                return CurvePoint(X,Y,1)
            
            
     def Add(self,Q):
         
            if self.z == 0:
                return Q
            if Q.z == 0:
                return self
            if self.x ==Q.x and self.y == Q.y:
                return self.Double()
            if self.x==Q.x:
                return CurvePoint(0, 1, 0)
            else:
                m = ((self.y-Q.y)*pow((self.x-Q.x),p-2,p))%p
                x = m**2-Q.x-self.x
                y = m*(Q.x-x)-Q.y
                return CurvePoint(x%p,y%p,1)
        
     def scal(self,m):
            P = self
            n = abs(m)
            res = CurvePoint(0, 1, 0)
            while n>0:
                 if n%2 == 1:
                     res = res.Add(P)
                 P = P.Double()
             
                 n = n//2
            if m<0 and res.z != 0:
                return CurvePoint(res.x%p, -res.y%p, 1)
            return res
        
     def sub(self,Q):
            if Q.z == 0:
                return self
            R = CurvePoint(Q.x%p, -Q.y%p, 1)
            return self.Add(R) 
        
        
        
 return CurvePoint


class EllipticCurve():
        def __init__(self, p1, a,b):
            #self.field = IntegersModP(p)
            p = mpz(p1)
            self.P = PolynomialsOver(p)
            self.a = mpz(a)%p
            self.b = mpz(b)%p
            self.EC = self.P([self.b,self.a,mpz(0),mpz(1)])
            self.p = mpz(p)
            if self.j_inv()%self.p == 0:
                print("Singular Curve")
                
        def j_inv(self):
            '''compute j-invariant of Elliptic Curve'''
            A = self.a
            B = self.b
            j = 1728*4*(A**3)*pow((4*(A**3) +27*(B**2)), self.p-2, self.p)
            return j%self.p
                
            
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
                r = self.P([pow(2,-1,self.p)])*Phi[m]*((Phi[m+2]*(Phi[m-1]**2))-(Phi[m-2]*(Phi[m+1]**2)))
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
           if ql>l//2:
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
                           return (-2*w)%l
                       else:
                           return (2*w)%l
               
           Q = Q1+Q2
      
           bound = int((l-1)/2)

           j = 1

           jP = P.Frobenius(q)
           B = jP
           while j <= bound:   
               if jP.x == Q.x:
                   break
               j +=1
               jP = jP+B


           if jP.y == Q.y:

               return j
           if jP.y == -Q.y:

               return -j
           else:
               print("fail")
              
        def a_mod2(self):
            q = self.p
            GCD= (self.P.x.powmod(q, self.EC)-self.P.x).gcd(self.EC)
            if GCD.deg() == 0:
                return 1
            return 0
           
    
  
        def Schoof(self):
            primeslist = primes(self.p)

            L = primeslist[-1]

            Phi = self.DivisionPoly(L)
            c1 = []
            c2 = []
            x,y = self.a_mod2(), 2
            c1.append(x)
            c2.append(2)
            for l in primeslist[1:]:
                 x =  self.odd_l_subroutine(Phi, l)
                 c1.append(x)
                 c2.append(l)

            a, Mod = CRT(c1,c2)

            if a>2*math.sqrt(self.p):

                a = a-Mod


            return "#E(F_p) = "+ str(self.p+1-a)
            
        
        def Eisensteins(self):
            E_4 = -48*self.a
            E_6 = 864*self.b
            return E_4%self.p, E_6%self.p
        
        
            
        
        
        def j_inv_prime(self):
            E_4,E_6 = self.Eisensteins()
            jprime = (-E_6*pow(E_4, self.p-2, self.p))*self.j_inv()
            return jprime%self.p
        
        def PrimeCheck(self, l, ModularPolynomial):
            t1 =time.time()
            pol=(self.P.x.powmod(self.p, ModularPolynomial)-self.P.x)
            t2 = time.time()
            #print("splitter time", t2-t1)
            t1 = time.time()
            GCD =pol.gcd(ModularPolynomial)
            t2 = time.time()
            #print("gcd compute time", t2-t1)
            lead = GCD.coefficients[-1]

            GCD = GCD*self.P([pow(lead, -1, self.p)])
            return GCD
        
        def IsogenousCurve(self,l, ModularPolynomial, gcd):
            j = self.j_inv()
            q = self.p
            if gcd.deg()==1:
             jtilde= -gcd.coefficients[0]%q
            if gcd.deg() == 2:
             Delta = gcd.coefficients[1]**2-(4)*gcd.coefficients[0]*gcd.coefficients[2]%q
             disc = modular_sqrt(Delta,q)%q
             jtilde= (-gcd.coefficients[1]-disc)*pow(2*gcd.coefficients[2], -1, q)%q
            jprime = self.j_inv_prime()
            Phi_x = ModularPolynomial.partial_x(j, jtilde)%q
            Phi_y = ModularPolynomial.partial_y(j, jtilde)%q
            jtilde_prime = -(jprime*Phi_x)*pow(l*Phi_y, -1, q)%q
            atilde = -(jtilde_prime**2)*pow(48*jtilde*(jtilde-1728), -1, q)%q
            btilde = -(jtilde_prime**3)*pow(864*(jtilde**2)*(jtilde-(1728)), -1, q)%q
            return j,jprime, jtilde, jtilde_prime, atilde, btilde
            
        def ElkiesProcedure(self, l, Phi_l, gcd):
             q = self.p
             j,jprime, jtilde, jtilde_prime, atilde, btilde = self.IsogenousCurve(l,Phi_l,gcd)
             t1 = jprime**2*Phi_l.partial_xx(j, jtilde)%q
             t2 = 2*l*jprime*jtilde_prime*Phi_l.partial_xy(j, jtilde)%q
             t3 = (l**2)*(jtilde_prime**2)*Phi_l.partial_yy(j, jtilde)%q
         
            
             E4_ql = -48*atilde%q
             E6_ql = 864*btilde%q

             E4, E6 = self.Eisensteins()
             vii14 = -(t1+t2+t3)*pow(jprime*Phi_l.partial_x(j, jtilde), -1, q)%q
             
             T1 = pow(2, -1, q)*l*vii14%q
             T2 = pow(4, -1, q)*l*((E4**2)*pow(E6, -1, q)-l*((E4_ql**2)*pow(E6_ql, -1, q)))%q
             T3 = pow(3, -1, q)*l*(E6*pow(E4, -1, q)-l*(E6_ql*pow(E4_ql, -1, q)))%q
             '''coefficient of x^((l-1)/2-1)'''
             p1 = (T1+T2+T3)%q

            
             atilde = (l**4)*atilde%q
             btilde = (l**6)*btilde%q

             d = int((l-1)/2)
             
             C_k, C_ktilde = self.WeierstrassCoeffecients(d, self.a, self.b, atilde, btilde)
             
             coeffs = self.KernelPolyCoefficients2(p1,d,C_k,C_ktilde)
             return self.P(coeffs)
             
       
        def ElkiesProcedureV2(self, l, Phi,gcd):
            q = self.p
            E_4 = -self.a*pow(3, -1, q)%q
            E_6 = -self.b*pow(2, -1, q)%q
            Delta = (E_4**3-E_6**2)*pow(1728, -1, q)%q
            j = self.j_inv()
            if gcd.deg()==1:
             g= -gcd.coefficients[0]
            if gcd.deg() == 2:
             disc = gcd.coefficients[1]**2-4*gcd.coefficients[0]*gcd.coefficients[2]%q
             disc = (modular_sqrt(disc,q))%q
             g= (-gcd.coefficients[1]-disc)*pow(2*gcd.coefficients[2], -1, q)
            Dg = g*Phi.partial_x(g,j)%q
            Dj = j*Phi.partial_y(g,j)%q
            s = int(12/math.gcd(l-1,12))
            Delta_l = (g**(int(12/s)))*Delta*pow(l**12, -1, q)%q
            if Dj == 0:
                print("Dj==0")
                E4_l = E_4*pow(l**2, -1, q)%q
                atilde = -3*(l**4)*E4_l%q
                jtilde = (E4_l**3)*pow(Delta_l, -1, q)%q
                btilde = (2*(l**6))*(modular_sqrt((((jtilde-(1728))*Delta_l)).n, q))%q
                p1 = 0
            else:
                E_2 =-12*E_6*Dj*pow(s*E_4*Dg, -1, q)%q
                gprime = -s*pow(12, -1, q)*E_2*g%q
                jprime = -(E_4**2)*E_6*pow(Delta, -1, q)%q
                E_0 = E_6*pow(E_4*E_2, -1, q)%q
                
                Dgprime = gprime*Phi.partial_x(g,j)+g*(gprime*Phi.partial_xx(g,j)+jprime*Phi.partial_xy(g,j))%q
                Djprime = jprime*Phi.partial_y(g,j)+j*(jprime*Phi.partial_yy(g,j)+gprime*Phi.partial_xy(g,j))%q
                
                E0_prime = ((-s*pow(12,-1,q))*Dgprime-E_0*Djprime)*pow(Dj,-1,q)%q
                E4_l = pow(l**2, -1, q)*(E_4-E_2*(12*E0_prime*pow(E_0,-1,q) +(6*(E_4**2)*pow(E_6,-1,q))-(4*E_6*pow(E_4, -1,q)))+E_2**2)%q
                jtilde = (E4_l**3)*pow(Delta_l, -1,q)%q
                f = (l**s)*pow(g,-1, q)
                fprime = s*E_2*f*pow(12, -1,q)%q
                Dgstar = Phi.partial_x(f,jtilde )%q
                Djstar = Phi.partial_y(f, jtilde)%q
                jtilde_prime = -fprime*Dgstar*pow(l*Djstar, -1, q)%q
                E6_l = -E4_l*jtilde_prime*pow(jtilde,-1,q)%q
                atilde = -3*(l**4)*E4_l%q
                btilde = -2*(l**6)*E6_l%q
                p1 = -l*E_2%q
            d = int((l-1)/2)
            Ck,Ck_tilde = self.WeierstrassCoeffecients(d, (self.a), (self.b), atilde, btilde)
            coeffs = self.KernelPolyCoefficients2(p1,d,Ck,Ck_tilde)
            return self.P(coeffs)
                
        def WeierstrassCoeffecients(self,d,a,b,atilde, btilde):
            q = self.p
            '''compute coefficients of weierstrass p function for E and E-tilde'''
            Ck = [-a*pow(5,-1, q)%q, -b*pow(7, -1, q)%q]
            Ck_tilde = [-atilde*pow(5,-1, q)%q, -btilde*pow(7, -1, q)%q]
            for k in range(3,d+1):
                coeffk = 3* pow((k-2)*(2*k+3), -1, q)%q
                c_k = coeffk*sum(Ck[j-1]*Ck[k-2-j] for j in range(1,k-1))%q
                c_ktilde = coeffk*sum(Ck_tilde[j-1]*Ck_tilde[k-2-j] for j in range(1,k-1))%q
                Ck.append(c_k)
                Ck_tilde.append(c_ktilde)
            return Ck, Ck_tilde
        


        def KernelPolyCoefficients2(self, p1, d, Ck, Ck_tilde):
            '''Compute Coefficients for factor of division polynomial that
            vanhishes on eigenspace of the frobenius endomorphism'''
            Ak = [-p1 * pow(2, -1, self.p)]
            Ak = Ak + [
                (-(Ck_tilde[j - 1] - (2 * d + 1) * Ck[j - 1]) * pow(((2 * j + 1) * (2 * j + 2)), -1, self.p) % self.p)
                for j in range(1, d)]
            C_series = [self.P.one]
            multiplier = self.P(Ak)
            A = self.P.one + self.P(Ak).shift(1)
            s = self.P(Ck)
            C = self.P([0] + Ck)
            C_series.append(C)
            j = 2
            g = multiplier
            h = s
            mod = self.P.x.shift(d)
            while j <= d:
                modj = self.P.x.shift(d - j)
                g = g * multiplier % modj
                q = self.P([pow((mpz(math.factorial(j))), -1, self.p)]) * g.shift(j)
                A = A + q
                h = h * s % modj
                C_series.append(h.shift(j))
                j = j + 1
            Coeffs = [0] * (d - 1) + [-p1 * pow(2, -1, self.p) % self.p] + [1]
            for i in range(1, d + 1):
                fl = A.coefficients[i]
                for k in range(1, i + 1):
                    vec = [
                        mpz(math.factorial(d - i + k) // (math.factorial(k - j) * math.factorial(d - i + j))) % self.p *
                        C_series[k - j].coefficients[j] for j in range(0, k)]
                    c = sum(vec) % self.p
                    fl = fl - c * Coeffs[d - i + k] % self.p
                Coeffs[d - i] = fl
            return Coeffs

        def FrobeniusOrder(self,l,Phi):
            '''compute the order of frobenius r in the case of an atkin prime l so that 
            the lth modular polynomial splits over F_p^r'''
            b = modular_sqrt(self.p,l)
            if b == 0:
                b = -1
            else:
                b = 1
            k=0
            h = self.P.x
            if (l in [37,61,73]):
               if ((h.powmod(self.p**2, Phi))-self.P.x).gcd(Phi).deg() == l+1:
                   return 2
               else:
                   if b == 1:
                     return (l+1)//2
                   if b == -1:
                       return l+1
            for i in range(1,l+2):

                if (l+1)%i == 0 and ((-1)**((l+1)//i)==b):
                    for c in range(1, i-k+1):
                     h = h.powmod(self.p, Phi)
                    k=i
                    g = (h-self.P.x).gcd(Phi)
                    if g.deg()==l+1:
                      return i
            
        def Atkin(self, l,Phi,r):
           '''Process an Atkin prime'''
           
           Pl = PolynomialsOver(l)
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
               
               g1 = gamma_r.polynomial.coefficients[0]
               x1_sq = (self.p*(g1+1)*pow(2,-1,l))%l
               if x1_sq  == 0:
                   
                   Traces.append(0)
               x1 = modular_sqrt(x1_sq,l)
           
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

            PointFactory = pointfactory(self.a, self.b, self.p)
           
            T1, T2 = Atkin1, Atkin2

            if T2 == []:
                m1 = T1[0][1]
                P = PointFactory()
              
                for t1 in T1[0][0]:
                    t,_ = CRT([t1,t3], [m1,m3])
                    if abs(t)>2*math.sqrt(self.p):
                        t = t-_

                    if P.scal(self.p+1-t).z == 0:
                        return t
                            
            m1, m2 = T1[0][1],T2[0][1]

            R1, R2 = [], []
            R1, R2 = sorted(R1, key = abs), sorted(R2, key = abs)
        
            
            for t1,m1 in T1:
             r1 = ((t1-t3)*pow(m2*m3, -1, m1))%m1
             if r1>m1/2:
                 R1.append(r1-m1)
             else:
                 R1.append(r1)
            for t2, m2 in T2:
             r2 = ((t2-t3)*pow(m1*m3, -1,m2))%m2
             R2.append(r2)
             R2.append(r2-m2)
             i =1
            while i<40:
             #print("loop iters"+str(i))
             P = PointFactory()
             Q = P.scal(self.p+1-t3)
             m3P = P.scal(m3)
             m1m3P = m3P.scal(m1)
             m2m3P = m3P.scal(m2) 
             Q_r1= {}
          
             for r1_i in R1:
                Q_r1_i = Q.sub(m2m3P.scal(r1_i))
                if Q_r1_i.z != 0:
                 Q_r1[Q_r1_i.x] = Q_r1_i.y, r1_i
                if Q_r1_i.z == 0:
                    Q_r1["O"] = r1_i
         
             for r2_i in R2:
                 r2_im1m3P = m1m3P.scal(r2_i)
                 if r2_im1m3P.z != 0:
               
                  if r2_im1m3P.x%self.p in Q_r1:
                     y_r_1_i, r_1_i = Q_r1[r2_im1m3P.x]
                     if  (y_r_1_i- r2_im1m3P.y)%self.p == 0:
                         t = t3+m3*(m2*r_1_i+m1*r2_i)
                         return t
                     else:
                      
                         t = t3+m3*(m2*r_1_i-m1*r2_i)
                         return t
                 if  r2_im1m3P.z%self.p == 0:
                    if "O" in Q_r1:
                        r_1_i = Q_r1["O"]
                        t = t3+m3*(m2*r_1_i+m1*r2_i)
                        return t
             i+=1
            print("Failed at BSGS step")
                       
            
            
            
                
                
        
       
            
        def SEA(self,Canonical = True, filter = True):

            primeslist = primes(self.p)
            x,y = self.a_mod2(), 2
            prod = 2
            E1 = [x]
            E2 = [y]
            
            A1 = []
            A2 = []
            number_of_traces = 0
            counter  =1
            elkiescount = 0
            atkincount = 0
            val =4*math.sqrt(self.p)
            while prod <=  val:
                l = Primes[counter]
                t1 = time.time()
                if Canonical and l<= 197:
                 Phi_l = CanonicalModularPolynomial(self.p, l)
                else:
                    Phi_l = ModularPolynomial(self.p, l)
                t2 = time.time()

                j = self.j_inv()
                t1 = time.time()
                Phi = Phi_l.eval_at_y(j)
                t2 = time.time()

                t1 =time.time()
                GCD = self.PrimeCheck(l, Phi)
                t2 = time.time()

                if GCD.deg()==l+1:
                    counter += 1
                   # print("pathalogical")
                    continue
                
                if GCD.deg() == 0:

                  if filter:
                       if ( l>30 and ( l!=73 or l != 61)):
                        counter +=1
                        print(str(l)+' rejected')
                        continue
                  atkincount += 1
                  s1 = time.time()
                  r= self.FrobeniusOrder(l,Phi)
                  s2 = time.time()

                  traces = self.Atkin(l,Phi,r)

                  number_of_traces += len(traces)
                  A1.append(traces)
                  A2.append(l)
                  counter +=1
                  prod = l*prod
                
                if GCD.deg()>0:

                 elkiescount+=1
                 if Canonical and l<= 197:
                  F_l = self.ElkiesProcedureV2(l, Phi_l, GCD)
                 else:

                  F_l = self.ElkiesProcedure(l, Phi_l, GCD)

                 Ltorring = GenQuotientRing(F_l, self.P)
                 Ltor = GenFracField(Ltorring)
                 Point = GenBasePoint(self.EC, self.a,self.b)
                 P = Point(Ltor.x, Ltor.one)
                 FrobP = P.Frobenius(self.p)
                 if GCD.deg() == 2:

                    
                    bound = (l-1)//2
                    j = 1
                    jP = P
                    x_j, y_j = jP.x, jP.y
                    light = 0
                    while j <= bound:   
                        if (x_j-FrobP.x).num.polynomial.gcd(F_l).deg()>0:

                         light = 1
                         break
                        
                        jP = jP+P
                        j = j+1
                        x_j, y_j = jP.x, jP.y 
                    if light == 0:
                        print('no eigenvalue found, skipping')
                        counter += 1
                        continue
                    if (y_j-FrobP.y).num.polynomial.gcd(F_l).deg()>0:

                      t= (j+self.p*pow(j, l-2, l))%l

                    else:

                        t= (-j+self.p*pow(-j,l-2,l))%l

                  
            
                    
                 else:

                  w = modular_sqrt(self.p,l)
                  wP = P.multiply(w)
                  if (wP.y-FrobP.y).num.polynomial.gcd(F_l).deg()>0:

                     t =( 2*w%l)

                  else:
                      t = ( -2*w%l)

                 E1.append(t)
                 E2.append(l)
                 counter +=1
                 prod = l*prod
                 
                  
                    
            Elkies = E1,E2
            if A1 == []:
                t,m = CRT(E1,E2)
                if abs(t) > 2*math.sqrt(self.p):
                    return "E(F_p)="+str(self.p+1-(t-m))
                return "E(F_p)="+str(self.p+1-t)
            Atkin = A1,A2
         
            Atkin1, Atkin2 = self.AtkinTraces(Atkin, number_of_traces)
            t = self.BSGS(Elkies, Atkin1, Atkin2)
            
            return "#E(F_p) = " + str(self.p+1-t)

            # noinspection PyUnreachableCode



        def BabyStep(self):
                # noinspection PyTypeChecker
                m = rand.randint(math.floor(math.sqrt(math.sqrt(self.p))), math.ceil(math.sqrt(math.sqrt(self.p))+10))
                PointFactory = pointfactory(self.a, self.b, self.p)
                L = 1
                thresh = 4*math.sqrt(self.p)
                while L<thresh:
                    P = PointFactory()
                    jPs = [(P.scal(j),j) for j in range(0,m+1)]
                    jPsx = [point[0].x for point in jPs]
                    jPsy = [(point[0].y, point[1]) for point in jPs]
                    D = dict(zip(jPsx, jPsy))
                    Q = P.scal(self.p+1)
                    twomP = P.scal(2*m)
                    k = 0
                    M = 0
                    while k<self.p+2*math.sqrt(self.p):

                        Z = Q.Add(twomP.scal(k))
                        if Z.x in D.keys():
                            if Z.y == D[Z.x][0]:
                                j = D[Z.x][1]
                            else:
                                j = - D[Z.x][1]
                            M = self.p+1+(2*m*k)-j
                            break
                        k += 1
                    if M == 0:

                        continue
                    primess = primesfac(M)

                    prod = 1
                    for prim in primess:
                        if P.scal(M//prim).z%self.p == 0:
                            M = M//prim

                    L = (M*L)//math.gcd(M,L)
                for n in range(int(self.p+1-2*math.sqrt(self.p)), int(self.p+1+2*math.sqrt(self.p))):
                    if n%L == 0:
                        return "#E(F_p) = " + str(n)

        def findpoint(self, divs):
            Pointfac = pointfactory(self.a, self.b, self.p)

            while True:
                indicator = 0
                P = Pointfac()
                for d in divs[:-1]:
                    if P.scal(d).z == 0:
                        indicator = 1
                        break

                if indicator == 0:
                    if P.scal(divs[-1]).z == 0:
                        return P
                else:
                    continue






