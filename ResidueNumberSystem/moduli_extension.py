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



def RNS2Integer(rns,moduli):
    n = len(moduli)
    result = 0

    P = 1
    for i in range(n):
        P *= moduli[i]

    v = [0]*n
    for i in range(n):
        v[i] = P//moduli[i]

    q = [0]*n
    for i in range(n):
        q[i] = Inverse(v[i],moduli[i])

    for i in range(n):
        result += rns[i]*v[i]*q[i]


    return result % P

def Integer2RNS(integer,moduli):
    n = len(moduli)
    rns = [0]*n

    for i in range(n):
        rns[i] = integer % moduli[i]

    return rns  


def Multiplier(moduli):

    n = len(moduli)

    P = 1
    for i in range(n):
        P *= moduli[i]

    v = [0]*n
    for i in range(n):
        v[i] = P//moduli[i]

    q = [0]*n
    for i in range(n):
        q[i] = Inverse(v[i],moduli[i])


    M = [v[i]*q[i] % P for i in range(n)]
    return M
    

def Linear_transform(old_rns):
    a1,a2,a3 = old_rns
    P = a1*a2*a3
    a4 = ((6*a1 + 17*a2 + 39*a3) % P) % 41
    a5 = ((15*a1 + 24*a2 + 21*a3) % P) % 43

    return [a1,a2,a3,a4,a5]


def Conversion_routine(old_rns):
    integer = RNS2Integer(old_rns,old_moduli)
    new_rns = Integer2RNS(integer,new_moduli)

    print("Integer: ",integer)
    print("Old RNS: ",old_rns)
    print("New RNS: ",new_rns)
    print("\n")



old_moduli = [13,23,37]
new_moduli = [13,23,37,41,43]


M = Multiplier(old_moduli)

P = 1
for i in range(len(old_moduli)):
    P *= old_moduli[i]



def Conversion_with_Multiplier(old_rns):
    new_rns = old_rns + [0]*(len(new_moduli)-len(old_rns))

    for i in range(len(old_moduli) , len(new_moduli)):
        new_rns[i] = sum([old_rns[j]*M[j] for j in range(len(old_moduli))]) % P
        new_rns[i] = new_rns[i] % new_moduli[i]

    return new_rns


def Conversion_with_Linear_Transformation(old_rns):
    new_rns = Linear_transform(old_rns)
    return new_rns



old_rns = [4,2,3]
new_rns = Conversion_with_Linear_Transformation(old_rns)

print(RNS2Integer(old_rns,old_moduli))
print(RNS2Integer(new_rns,new_moduli))
