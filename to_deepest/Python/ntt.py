from modular import *

def CT_4pts(P, t0, t1, t2, t3, omega_1, omega2, omega3, mu_omega_1, mu_omega2, mu_omega3):
    # First layer
    t2 = Barrett_mul(t2, omega_1, P, mu_omega_1)
    t3 = Barrett_mul(t3, omega_1, P, mu_omega_1)

    temp0 = mod_add(t0, t2, P)
    t2 = mod_sub(t0, t2, P)
    t0 = temp0

    temp1 = mod_add(t1, t3, P)
    t3 = mod_sub(t1, t3, P)
    t1 = temp1

    # Second layer
    t1 = Barrett_mul(t1, omega2, P, mu_omega2)
    t3 = Barrett_mul(t3, omega3, P, mu_omega3)

    temp0 = mod_add(t0, t1, P)
    t1 = mod_sub(t0, t1, P)
    t0 = temp0

    temp1 = mod_add(t2, t3, P)
    t3 = mod_sub(t2, t3, P)
    t2 = temp1

    return t0, t1, t2, t3

def GS_4pts(P, t0, t1, t2, t3, omega_1, omega2, omega3, mu_omega_1, mu_omega2, mu_omega3):
    # First layer
    temp0 = mod_add(t0, t1, P)
    t1 = mod_sub(t0, t1, P)
    t0 = temp0

    temp1 = mod_add(t2, t3, P)
    t3 = mod_sub(t2, t3, P)
    t2 = temp1

    t1 = Barrett_mul(t1, omega2, P, mu_omega2)
    t3 = Barrett_mul(t3, omega3, P, mu_omega3)

    # Second layer
    temp0 = mod_add(t0, t2, P)
    t2 = mod_sub(t0, t2, P)
    t0 = temp0

    temp1 = mod_add(t1, t3, P)
    t3 = mod_sub(t1, t3, P)
    t1 = temp1

    t2 = Barrett_mul(t2, omega_1, P, mu_omega_1)
    t3 = Barrett_mul(t3, omega_1, P, mu_omega_1)

    return t0, t1, t2, t3



def forward_10(temp,ptr,P,omegas,mu_omegas):
    n = 1 << 10
    ptrf = ptr
    ptrg = ptrf + n


    for depth in range(0, 10, 2):
        logn_top = 10 - depth
        n = 1 << logn_top
        hhn = n >> 2
        power_index = 1 << depth

        for j in range(0, 1 << depth):
            for i in range(0, hhn):
                temp[ptrf + j*n + i], temp[ptrf + j*n + hhn + i], temp[ptrf + j*n + 2*hhn + i], temp[ptrf + j*n + 3*hhn + i] = CT_4pts(P, temp[ptrf + j*n + i], temp[ptrf + j*n + hhn + i], temp[ptrf + j*n + 2*hhn + i], temp[ptrf + j*n + 3*hhn + i], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1])
                temp[ptrg + j*n + i], temp[ptrg + j*n + hhn + i], temp[ptrg + j*n + 2*hhn + i], temp[ptrg + j*n + 3*hhn + i] = CT_4pts(P, temp[ptrg + j*n + i], temp[ptrg + j*n + hhn + i], temp[ptrg + j*n + 2*hhn + i], temp[ptrg + j*n + 3*hhn + i], omegas[power_index], omegas[power_index<<1], omegas[(power_index<<1)+1], mu_omegas[power_index], mu_omegas[power_index<<1], mu_omegas[(power_index<<1)+1])
            power_index += 1
    
def backward_10(temp,ptr,P,omegas,mu_omegas):
    n = 1 << 10
    ptrf = ptr
    ptrg = ptrf + n


    for depth in range(8, -1, -2):
        logn = 10 - depth
        n = 1 << logn
        hhn = n >> 2
        power_index = 2*(1 << depth)-1

        for j in range(0, 1 << depth):
            for i in range(0, hhn):
                temp[ptrf + j*n + i], temp[ptrf + j*n + hhn + i], temp[ptrf + j*n + 2*hhn + i], temp[ptrf + j*n + 3*hhn + i] = GS_4pts(P, temp[ptrf + j*n + i], temp[ptrf + j*n + hhn + i], temp[ptrf + j*n + 2*hhn + i], temp[ptrf + j*n + 3*hhn + i], -omegas[power_index], -omegas[(power_index<<1)+1], -omegas[power_index<<1], -(mu_omegas[power_index]-2), -(mu_omegas[(power_index<<1)+1]-2), -(mu_omegas[power_index<<1]-2))
                temp[ptrg + j*n + i], temp[ptrg + j*n + hhn + i], temp[ptrg + j*n + 2*hhn + i], temp[ptrg + j*n + 3*hhn + i] = GS_4pts(P, temp[ptrf + j*n + i], temp[ptrf + j*n + hhn + i], temp[ptrf + j*n + 2*hhn + i], temp[ptrf + j*n + 3*hhn + i], -omegas[power_index], -omegas[(power_index<<1)+1], -omegas[power_index<<1], -(mu_omegas[power_index]-2), -(mu_omegas[(power_index<<1)+1]-2), -(mu_omegas[power_index<<1]-2))
            power_index -= 1

    