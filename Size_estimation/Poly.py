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




def modular_exponentiation(base, exp, mod):
    """快速模指數計算 base^exp % mod"""
    result = 1
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        base = (base * base) % mod
        exp //= 2
    return result



# polynomial sampler
import numpy as np
import math
def prng(sigma):
    return round(((np.random.normal(0, sigma, 1))[0]))



def Falcon_poly_generator(n):
    q = 12289
    sigma = math.sqrt(q/(2*n)) * 1.17

    f = []
    for i in range(n):
        f.append(prng(sigma))

    # Check for irreducibility

    # Check the inner product is small enough



    return f


def Coefficient_reduction(f,q): #Very time consuming
    for i in range(len(f)):
        f[i] = f[i] % q

def to_star(f):
    f_new = []

    for i in range(len(f)):
        f_new.append(f[i]*((-1)**i))

    return f_new

# LinearConv.
def Linear_Conv(f,g,p):
    if p != 0:
        
        n = len(f)
        result = [0]*(2*n-1)
        for i in range(n):
            for j in range(n):
                result[i+j] += (f[i]*g[j] % p)
        return result

    else:

        n = len(f)
        result = [0]*(2*n-1)
        for i in range(n):
            for j in range(n):
                result[i+j] += (f[i]*g[j])
        return result



def NegaCyclic_Conv(f,g,p):
    if p != 0:

        n = len(f)
        result = Linear_Conv(f,g,p)
        for i in range(n-1):
            result[i] -= result[i+n] % p
        return result[:n]

    else:

        n = len(f)
        result = Linear_Conv(f,g,p)
        for i in range(n-1):
            result[i] -= result[i+n] 
        return result[:n]

def Pointwise_Multiply(f,g,p):
    if p != 0:
        n = len(f)
        result = [0]*n
        for i in range(n):
            result[i] = (f[i]*g[i]) % p
        return result
    else:
        n = len(f)
        result = [0]*n
        for i in range(n):
            result[i] = (f[i]*g[i])
        return result


# logn = 10
# f = Falcon_poly_generator(n=(1<<logn))
# g = Falcon_poly_generator(n=(1<<logn))
# h = NegaCyclic_Conv(f,g,12289)
# 
# f_ntt = NTT(f, 1<<logn, root = 12282, mod = 12289)
# g_ntt = NTT(g, 1<<logn, root = 12282, mod = 12289)
# 
# h_ntt = Pointwise_Multiply(f_ntt, g_ntt, 12289)
# 
# f_f = INTT(f_ntt, 1<<logn, root = 12282, mod = 12289)
# g_g = INTT(g_ntt, 1<<logn, root = 12282, mod = 12289)
# 
# h_h = INTT(h_ntt, 1<<logn, root = 12282, mod = 12289)
# 
# 
# print(f == f_f)
# print(g == g_g)
# print(h == h_h)
# print(h_h)
# print(h)
