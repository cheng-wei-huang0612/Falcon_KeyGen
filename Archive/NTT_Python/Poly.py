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


def NTT(poly, n, root = 12282, mod = 12289):
    """
    FFT/NTT 
    poly: 多項式係數列表。
    n: 點數（必須為 2 的次冪）。
    root: 原根。
    mod: 質數模數。
    """
    length = n  # 初始處理完整數據長度
    while length > 1:
        # 計算當前層的旋轉因子 wlen（模原根的冪次）
        wlen = modular_exponentiation(root, (mod - 1) // length, mod)

        half = length // 2  # 每組的半長度
        for i in range(0, n, length):  # 每次處理一個完整的組
            w = 1  # 當前區塊的旋轉因子初始化
            for j in range(half):
                # 蝴蝶結構計算
                u = poly[i + j]  # “左分支”數據
                v = poly[i + j + half]  # “右分支”數據

                # 合併結果計算
                poly[i + j] = (u + v * w) % mod  # 左分支更新為和
                poly[i + j + half] = (u - v * w) % mod  # 右分支更新為差

                if poly[i + j + half] < 0:  # 確保模數下非負
                    poly[i + j + half] += mod

                # 更新旋轉因子
                w = (w * wlen) % mod

        # 每層結束後，縮小處理範圍
        length //= 2

    return poly

def INTT(poly, n, root=12282, mod=12289):
    """
    反數論傅立葉變換 (Inverse NTT, INTT)
    poly: 點值表示的多項式。
    n: 點數（必須為 2 的次冪）。
    root: 原根。
    mod: 質數模數。
    """
    # 計算 root 的模反元素
    root_inv = modular_exponentiation(root, mod - 2, mod)
    
    length = 2  # 初始處理的長度，從最小開始（與 NTT 的 length 方向相反）
    while length <= n:
        # 計算當前層的旋轉因子 wlen
        wlen = modular_exponentiation(root_inv, (mod - 1) // length, mod)
        
        half = length // 2  # 每組的半長度
        for i in range(0, n, length):  # 每次處理一個完整的組
            w = 1  # 當前區塊的旋轉因子初始化
            for j in range(half):
                # 蝴蝶結構計算
                u = poly[i + j]  # “左分支”數據
                v = poly[i + j + half]  # “右分支”數據

                # 合併結果計算（注意和 NTT 相反）
                poly[i + j] = (u + v) % mod  # 左分支更新為和
                poly[i + j + half] = ((u - v) * w) % mod  # 右分支更新為差並乘旋轉因子

                if poly[i + j + half] < 0:  # 確保模數下非負
                    poly[i + j + half] += mod

                # 更新旋轉因子
                w = (w * wlen) % mod

        # 每層結束後，擴大處理範圍
        length *= 2

    # 除以 n，乘以模數下的 n 的逆元
    n_inv = modular_exponentiation(n, mod - 2, mod)
    for i in range(n):
        poly[i] = (poly[i] * n_inv) % mod

    return poly


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

