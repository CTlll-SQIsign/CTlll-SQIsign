import random
from sage.all import is_prime, GF, sqrt, ZZ, gcd, i, floor, QQ

def Cornacchia_prime(M,d=1):
    if(M==1):
        return(1,0)
    if(M==2):
        if(d==1):
            return(1,1)
        elif(d==2):
            return(0,1)
        else:
            return False
    if(ZZ(M).is_prime()):
        one = GF(M)(-d)
        u = one.sqrt()
        if(u not in GF(M)):
            return False
        a = M
        b = ZZ(u)
        r = M
        while(r**2>=M):
            a = b
            b = r
            r = a%b
        x = r
        y = ZZ(sqrt((M-x**2)//d))
        if(x**2 + d*y**2 == M):
            return(x,y)
        else:
            return(False)

#Extended Cornacchia
def Corvo_imperiale(M):
    bad = 9758835872005600424541432018720696484164386911251445952809873289484557900064682994157437552186699082331048989
    if(gcd(bad,M)!=1):
        return(False)

    good_primes=[2,5,13,17,29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157, 173, 181, 193, 197, 229, 233, 241, 257, 269, 277, 281, 293, 313, 317, 337, 349, 353, 373, 389, 397, 401, 409, 421, 433, 449, 457, 461]
    n = M
    n_val = [0 for i in range (len(good_primes))]
    for k in range(len(good_primes)):
        p = good_primes[k]
        n_val[k] = n.valuation(p)
        n = n // (p**n_val[k])

    if(n %4 == 3):
        return(False)
    if not ((ZZ(n)).is_prime() or n == 1):
        return False
    x,y = Cornacchia_prime(n,1)
    c = ZZ[i](x+i*y)
    for k in range(len(good_primes)):
        if(n_val[k]!= 0):
            x_p,y_p = Cornacchia_prime(good_primes[k],1)
            c_p = ZZ[i](x_p+i*y_p)
            c = c*c_p**n_val[k]
    return(abs(c[0]),abs(c[1]))

#Argument n is the target integer, B the quaternion algebra
def representInteger(n, B, t=1000):
    ib,jb,kb = B.gens()
    p = B.ramified_primes()[0]
    bi = floor((n//p).sqrt())
    if(bi<=0):
        return False
    for it in range(t):
        x = random.randrange(0,bi)
        bj = floor(((n-x*x)//p).sqrt())
        y = random.randrange(0,bj)
        r = Corvo_imperiale(n-p*(x*x+y*y))
        if(r!= False):
            res = r[0]+ib*r[1]+jb*x+kb*y
            assert(res.reduced_norm()==n)
            if(gcd(ZZ(res.reduced_norm()),n*n) == n):
                return(res)
    return(False)
