

# We use this file to generate: 
# - data.c
# - data.h

# Key ingredients:

logn_top = 10
P = [33550337, 33540097, 33538049, 33533953]
omega_2048 = [33547300, 33533628, 33510748, 33448611]


# data contain:
# - logn
# - n
# - P[4]
# - omegas[8192]
# - mu_omegas[8192]
# 
# - f[1024],g[1024] to be sampled

# Here begin the computation
def MOD(a, q):
    result = a % q
    if result >= (q>>1):
        result = result - q
    return result


# Bit-reverse index
def bitrev(i):
    rev = 0
    for j in range(10):
        rev = (rev << 1) | (i & 1)
        i >>= 1
    return rev



n = 1 << logn_top
print(bitrev(0))

# polynomial sampler
import numpy as np
import math
def prng(sigma):
    return int(abs((np.random.normal(0, sigma, 1))[0]))

def Falcon_poly_generator(n):
    q = 12289
    sigma = math.sqrt(q/(2*n)) * 1.17

    f = []
    for i in range(n):
        f.append(prng(sigma))

    # Check for irreducibility

    # Check the inner product is small enough
    
    return f

def Mu(b,q):
    mu = 2*int(((b << 31)/(2*q))+0.5)
    return mu

def omegas_gen(omega,n,mod):
    """ 
    returns omegas[n] for each given omega and n
    """
    val = 1
    omegas = [1] * (n)

    for i in range(n):
        omegas[bitrev(i)] = val
        val = MOD((val * omega) , mod)

    return omegas

def mu_omegas_gen(omegas,mod):
    """
    returns mu_omegas[n] for given omegas and modulo mod
    """
    mu_omegas = [Mu(omega,mod) for omega in omegas]
    return mu_omegas


omegas = [0] * 4
mu_omegas = [0] * 4
for i in range(4):
    omegas[i] = omegas_gen(omega_2048[i],n,P[i])
    mu_omegas[i] = mu_omegas_gen(omegas[i],P[i])

f_poly = Falcon_poly_generator(n)
g_poly = Falcon_poly_generator(n)


# Write the header file
with open("data.h", "w") as header_file:
    header_file.write("#ifndef DATA_H\n")
    header_file.write("#define DATA_H\n\n")
    header_file.write("#include <stdint.h>\n\n")
    header_file.write("extern const uint16_t logn_top;\n")
    header_file.write("extern const int32_t P[4];\n")
    header_file.write("extern const int32_t omegas[%d];\n" % (4 * n))
    header_file.write("extern const int32_t mu_omegas[%d];\n" % (4 * n))
    header_file.write("extern const int8_t f[%d];\n"  % n)
    header_file.write("extern const int8_t g[%d];\n\n"% n)
    header_file.write("#endif // DATA_H\n")


# Write the file
with open("data.c", "w") as file:
    file.write("#include \"data.h\"\n\n")
    file.write("const uint16_t logn_top = %d;\n" % logn_top)
    file.write("const uint16_t n = %d;\n" % n)
    file.write("const int32_t P[4] = {%d,%d,%d,%d};\n" % (P[0], P[1], P[2], P[3]))
    
    # Flatten omegas array
    file.write("const int32_t omegas[%d] = {\n" % (4 * n))
    for i in range(4):
        for j in range(n):
            file.write("%d," % omegas[i][j])
            if (j + 1) % 8 == 0:
                file.write("\n\t")
    file.write("};\n")
    
    # Flatten mu_omegas array
    file.write("const int32_t mu_omegas[%d] = {\n" % (4 * n))
    for i in range(4):
        for j in range(n):
            file.write("%d," % mu_omegas[i][j])
            if (j + 1) % 8 == 0:
                file.write("\n\t")
    file.write("};\n")
    
    file.write("const int8_t f[%d] = {\n" % n)
    for i in range(n):
        file.write("%d," % f_poly[i])
        if (i + 1) % 8 == 0:
            file.write("\n\t")
    file.write("};\n")
    
    file.write("const int8_t g[%d] = {\n" % n)
    for i in range(n):
        file.write("%d," % g_poly[i])
        if (i + 1) % 8 == 0:
            file.write("\n\t")
    file.write("};\n")


# write data.py file
with open("data.py", "w") as file:
    file.write("logn_top = %d\n" % logn_top)
    file.write("P = [%d, %d, %d, %d]\n" % (P[0], P[1], P[2], P[3]))
    file.write("omega_2048 = [%d, %d, %d, %d]\n" % (omega_2048[0], omega_2048[1], omega_2048[2], omega_2048[3]))
    file.write("omegas = [\n")
    for i in range(4):
        for j in range(n):
            file.write("%d, " % omegas[i][j])
    file.write("]\n")
    file.write("mu_omegas = [\n")
    for i in range(4):
        for j in range(n):
            file.write("%d, " % mu_omegas[i][j])

    file.write("]\n")
    file.write("f = [")
    for i in range(n):
        file.write("%d, " % f_poly[i])
    file.write("]\n")
    file.write("g = [")
    for i in range(n):
        file.write("%d, " % g_poly[i])
    file.write("]\n")