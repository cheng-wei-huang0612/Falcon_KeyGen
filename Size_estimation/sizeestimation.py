from Poly import *
from makemoduli import *



logn = 10
bit_length = 30
num_of_moduli = 1

num_of_used_moduli_fg = [1, 1, 1, 2,  4,  8, 16, 28,  56, 108, 216]
num_of_used_moduli_FG = [1, 2, 4, 8, 12, 24, 44, 84, 160, 320, 216]





fg_max = [4.00, 10.99, 24.07, 50.37, 101.62, 202.22, 400.67, 794.17, 1576.87, 3138.35, 6307.52]
fg_max_std = [0.00, 0.08, 0.25, 0.53, 1.02, 1.87, 3.10, 4.98, 7.49, 12.25, 24.48]
FG_max = [19.61, 39.82, 78.20, 153.65, 303.49, 599.81, 1188.68, 2361.84, 4703.30, 9403.29, 6319.66]
FG_max_std = [0.49, 0.41, 0.73, 1.39, 2.38, 3.87, 6.04, 9.31, 14.77, 27.55, 24.51]


def test_with_stddev(mul_stddev):
    coverage_fg = [i * bit_length for i in num_of_used_moduli_fg]
    coverage_FG = [i * bit_length for i in num_of_used_moduli_FG]

    c_fg = [True for i in range(logn+1)]
    c_FG = [True for i in range(logn+1)]

    for i in range(logn+1):
        c_fg[i] = (coverage_fg[i] > (fg_max[i] + (fg_max_std[i] * mul_stddev)))
        c_FG[i] = (coverage_FG[i] > (FG_max[i] + (FG_max_std[i] * mul_stddev)))

    return c_fg, c_FG

(c_fg, c_FG) = test_with_stddev(5)
print(c_fg)
print(c_FG)
