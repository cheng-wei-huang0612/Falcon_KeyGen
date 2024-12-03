# This python file is used to 
# - Generate random polynomials 
# - Multiply two polynomials with school book multiplication

# Parameters:
logn = 3
n = 1 << logn
p = 33550337




# polynomial sampler
import numpy as np
def prng():
    return int(abs((np.random.normal(0, 3, 1))[0]))

def Polynomial_Generator(n):
    return [prng() for _ in range(n)]


# Convolutions
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



f = Polynomial_Generator(n)
g = Polynomial_Generator(n)
    
print("the first polynomial f is: ",f)
print("the second polynomial g is: ",g)

# Write them into files called "data.h" and "data.c"

correct = NegaCyclic_Conv(f,g,p)
print("the correct result is: ",correct)


with open("data.c", "w") as file:
    file.write("#include <stdint.h>\n")
    file.write("#define N %d\n" % n)
    
    # Write the f array
    file.write("int32_t f[N] = {")
    for i in range(n):
        file.write("%d" % f[i])
        if i < n - 1:
            file.write(", ")
    file.write("};\n")
    
    # Write the g array
    file.write("int32_t g[N] = {")
    for i in range(n):
        file.write("%d" % g[i])
        if i < n - 1:
            file.write(", ")
    file.write("};\n")

    # Write the correct result
    file.write("int32_t correct[N] = {")
    for i in range(n):
        file.write("%d" % correct[i])
        if i < n - 1:
            file.write(", ")


with open("data.h", "w") as file:
    file.write("#define N %d\n" % n)
    file.write("extern int32_t f[N];\n")
    file.write("extern int32_t g[N];\n")
    file.write("extern int32_t correct[N];\n")


