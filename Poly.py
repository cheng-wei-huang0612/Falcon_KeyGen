# Python program to demonstrate working of extended 

# function for extended Euclidean Algorithm 
def gcdExtended(a, b): 
    # Base Case 
    if a == 0 : 
        return b,0,1
             
    gcd,x1,y1 = gcdExtended(b%a, a) 
     
    # Update x and y using results of recursive 
    # call 
    x = y1 - (b//a) * x1 
    y = x1 
     
    return gcd,x,y 
     
 
def Inverse(to_be_inverted, p):
    tmp = gcdExtended(to_be_inverted,p)
    if not tmp[0] == 1:
        print("Not invertible")

    if tmp[1] < 0:
        return (tmp[1]+p)
    return tmp[1]


# polynomial sampler
import numpy as np
def prng():
    return int(abs((np.random.normal(0, 3, 1))[0]))
    
def Polynomial_Generator(n):
    return [prng() for _ in range(n)]



def Coefficient_reduction(f,q): #Very time consuming 
    for i in range(len(f)):
        f[i] = f[i] % q


# LinearConv.
def Linear_Conv(f,g,p):
    n = len(f)
    result = [0]*(2*n-1)
    for i in range(n):
        for j in range(n):
            result[i+j] += (f[i]*g[j] % p)
    return result

def NegaCyclic_Conv(f,g,p):
    n = len(f)
    result = Linear_Conv(f,g,p)
    for i in range(n-1):
        result[i] -= result[i+n] % p
    return result[:n]
