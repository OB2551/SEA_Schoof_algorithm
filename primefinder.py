import gmpy2
import random as rand
import math
def randomnbitprime(n):
    lower, upper = 2**(n-1), (2**n)-1
    while True:
     m = rand.randint(lower, upper)
     prime = gmpy2.next_prime(m)
     bits = math.floor(math.log(prime,2))+1
     if int(math.floor(math.log(prime,2))+1) == n:
         return prime

def checktype(s):
    alpha = "abcdefghijklmnopqrstuvwxyz%!$#^&*()-+~/'{}="
    for i in s:
        if i in alpha:
            return False
    return True