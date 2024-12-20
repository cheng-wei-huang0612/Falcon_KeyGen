def mod_mul_standard(a, b, p):
    result = (a * b) % p
    if result > (p>>1):
        result -= p
    return result

def mod_add(a, b, p):
    result = (a + b) 
    if result > (p>>1):
        result -= p
    if result < -(p>>1):
        result += p
    return result

def mod_sub(a, b, p):
    result = (a - b) 
    if result > (p>>1):
        result -= p
    if result < -(p>>1):
        result += p
    return result


def Mu(b,q):
    mu = 2*int(((b << 31)/(2*q))+0.5)
    return mu

def Barrett_mul(a, b, p, mu_b):
    # Perform the multiplication
    z = (a * b) % (1 << 32)
    if z >= (1 << 31):
        z -= (1 << 32)  # Convert to signed 32-bit integer

    # Perform the multiplication for mu_b and get the high 32 bits
    temp = (a * (mu_b << 1)) % (1 << 64)
    t_high = (temp >> 32) % (1 << 32)
    if t_high >= (1 << 31):
        t_high -= (1 << 32)  # Convert to signed 32-bit integer
    
    result = z - t_high * p
    # Simulate 32-bit signed integer overflow for the result
    result %= (1 << 32)
    if result >= (1 << 31):
        result -= (1 << 32)  # Convert to signed 32-bit integer

    return result

def Barrett_red(A, N, V, i=25):
    # Step 1: Doubling multiplication, high half (sqdmulh equivalent)
    T0 = (A * V) >> (i + 1)
    
    # Step 2: Signed shift right with rounding (srshr equivalent)
    T1 = (T0 + (1 << (i - 2))) >> (i - 1)
    
    # Step 3: Subtract modulus (mls equivalent)
    z = A - (T1 * N)
    
    # Ensure result is in range -N/2 ≤ z < N/2
    if z >= N // 2:
        z -= N
    elif z < -N // 2:
        z += N
    
    return z


