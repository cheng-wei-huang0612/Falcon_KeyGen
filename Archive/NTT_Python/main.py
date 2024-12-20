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



