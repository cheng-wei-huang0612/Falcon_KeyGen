def check_arrays(a,b,p):
    if a == b:
        print("exactly equal")
        return True
    print("not exactly equal")

    if len(a) != len(b):
        print("lengths not equal")
        return False


    for i in range(len(a)):
        if (a[i] - b[i]) % p != 0:
            print("element", i, "not equal even in modulus class")
            return False
    print("But equal in modulus class")
    return True


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
    
def number_order(temp, p):
    maximum = max([abs(x) for x in temp])

    degree = 0
    while True:
        if degree*(p>>1) < maximum < (degree+1)*(p>>1):
            return degree
        degree += 1