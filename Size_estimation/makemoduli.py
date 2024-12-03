


# Initialization function: compute the contents of primetest_pp[].
def mk_primetest():
    primetest_pp = []

    max_p = 65536   # ceil(sqrt(2^31))
    max_q = 256     # ceil(sqrt(max_p))
    tab = [True] * max_p
    tab[0] = False
    tab[1] = False
    for i in range(2, max_q):
        if tab[i]:
            for j in range(2*i, max_p, i):
                tab[j] = False
    for i in range(0, max_p):
        if tab[i]:
            primetest_pp.append(i)
    
    return primetest_pp

    



# Function to test if a number is prime.
def is_prime(p,primetest_pp):
    if p < 2:
        return p == 2
    if (p & 1) == 0: # p is even
        return False
    for j in range(0, len(primetest_pp)):
        t = primetest_pp[j]
        if t*t > p:
            return True
        if (p % t) == 0:
            return False
    return True

def mk_moduli(bit_length = 31, num_of_moduli = 241):

    primetest_pp = mk_primetest()


    moduli = []
    p = 2**bit_length + 1
    while len(moduli) < num_of_moduli:
        if is_prime(p,primetest_pp):
            moduli.append(p)
        p -= 2048

    return moduli


