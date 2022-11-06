from EllipticCurve import*
import random

p = 197
A, B = random.randint(0, p), random.randint(0,p)
E = EllipticCurve(p, A, B)
ord1, ord2, ord3 = E.sea(), E.schoof(), E.baby_step()
print("SEA: ", ord1)
print("Schoof: ", ord2)
print("Baby-step giant-step: ", ord3)

p = 2**128-159
A, B = random.randint(0, p), random.randint(0,p)


E = EllipticCurve(p, A, B)
order = E.sea()
print("128-bit example with SEA: ", order)




#Example Output (will be random):
#SEA:  #E(F_p) = 201
#Schoof:  #E(F_p) = 201
#Baby-step giant-step:  #E(F_p) = 201
#128-bit example with SEA:  #E(F_p) = 340282366920938463470648655491441425069
