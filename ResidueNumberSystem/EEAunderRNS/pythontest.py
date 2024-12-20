from math import gcd

def extended_euclid(a, b):
    if b == 0:
        return (1, 0, a)
    x2, y2, g = extended_euclid(b, a % b)
    return (y2, x2 - (a // b)*y2, g)

def crt(x_list, p_list):
    P = 1
    for p in p_list:
        P *= p
    result = 0
    for (x_i, p_i) in zip(x_list, p_list):
        P_i = P // p_i
        inv_P_i = pow(P_i, p_i-2, p_i)
        result += x_i * P_i * inv_P_i
    return result % P

f_big = 13
g_big = 41
p_list = [3, 5, 7]
assert gcd(f_big, g_big) == 1

# 整數域中求解
#x,y,d = extended_euclid(f_big, g_big)
#assert d == 1
#G = x
#F = -y

f_i = [f_big % p for p in p_list]
g_i = [g_big % p for p in p_list]

G_i_list = []
F_i_list = []

for fi, gi, p in zip(f_i, g_i, p_list):
    G_i = G % p
    F_i = F % p
    lhs = (fi*G_i - gi*F_i) % p
    assert lhs == 1, f"Check failed under mod {p}: got {lhs}"
    G_i_list.append(G_i)
    F_i_list.append(F_i)

# 用CRT合併G
G_rebuilt = crt(G_i_list, p_list)
F_rebuilt = crt(F_i_list, p_list)
lhs_full = f_big * G_rebuilt - g_big * F_rebuilt

assert lhs_full == 1
print("Success!")
print(f"G={G_rebuilt}, F={F_rebuilt}, and {f_big}*G - {g_big}*F = {lhs_full}")
