alltypes = (type(1),int,type(1.0),sage.rings.real_mpfr.RealNumber,sage.rings.rational.Rational, float)
inttypes = (type(1),int)
rattypes = (type(1.0),sage.rings.real_mpfr.RealNumber,sage.rings.rational.Rational,float)

# - - - Declaring constants
pint = 256
logpint = 8
logpint = ceil( log(pint,2) )
rejection_num = 64 #number of samples in rnd
stepnum = 5040

# - - - Some utilities
def binlog(s):
    #Computes ceil for log_2(v) for v<2^256.
    r = ( s > (2 << (2**(logpint-1)) - 1 ) ) << (logpint-1)
    s >>= r
    for i in range(logpint - 2,0,-1):
        shift = ( s > ( (1 << (2**i)) - 1 ) ) << i
        s >>= shift
        r |= shift
    r |= (s >> 1)
    return r+1

def int_sqrt(s):
    """
    Integer square root algorithm without fixed point arithmetic.
    For debug purposes.
    """
    # // Zero yields zero
    # // One yields one
    s = int(s)

    iters = 0
    z = s>0
    s = z*s +(1-z)*1

    # // Initial estimate (must be too high)
    x1 =  2 << ( (binlog(s)//2) ) #2 << ( (log2(s)>>1) + 1 )
    a0, a1 = x1, 0
    for i in range(logpint+1):
        x0 = x1
        x1 = (x0 + s // x0) //2
        a0, a1 = x1, a0
    return min(a0,a1)

# - - - fixed-point arithmetic emulator
class fxp:
    def __init__(self,a, prec=pint, format=True):
        if format:
            if isinstance(a,rattypes):
                a = round( a*2**prec )
            elif isinstance(a,inttypes):
                a = a * 2**(prec)
        self.value = round( a )
        self.prec = prec

    def __eq__(self,other):
        if isinstance(other, inttypes):
            s = self.value / 2**self.prec
            return s == other
        if isinstance(other, rattypes):
            s = self.value / 2**self.prec
            return abs(s - other) < 2**-self.prec
        elif isinstance(other, type(self)):
            assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
            return self.value == other.value

    def __lt__(self,other):
        if isinstance(other, inttypes):
            s = self.value / 2**self.prec
            return s < other
        if isinstance(other, rattypes):
            s = self.value / 2**self.prec
            return s < other
        elif isinstance(other, type(self)):
            assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
            return self.value < other.value
        raise NotImplementedError(f"No output in lt: {type(other)}")

    def __le__(self,other):
        if isinstance(other, inttypes):
            s = self.value / 2**self.prec
            return s <= other
        if isinstance(other, rattypes):
            s = self.value / 2**self.prec
            return s <= other
        elif isinstance(other, type(self)):
            assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
            return self.value <= other.value
        raise NotImplementedError(f"No output in lt: {type(other)}")

    def __gt__(self,other):
        if isinstance(other, inttypes):
            s = self.value / 2**self.prec
            return s > other
        if isinstance(other, rattypes):
            s = self.value / 2**self.prec
            return s > other
        elif isinstance(other, type(self)):
            assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
            return self.value > other.value
        raise NotImplementedError(f"No output in gt: {type(other)}")

    def __gt__(self,other):
        if isinstance(other, inttypes):
            s = self.value / 2**self.prec
            return s >= other
        if isinstance(other, rattypes):
            s = self.value / 2**self.prec
            return s >= other
        elif isinstance(other, type(self)):
            assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
            return self.value >= other.value
        raise NotImplementedError(f"No output in gt: {type(other)}")

    def __neg__(self):
        return fxp( -self.value, self.prec )

    def __abs__(self):
        return fxp( abs(self.value), self.prec )

    def __add__(self,other):
        if isinstance(other, alltypes) and not isinstance(other, type(self)):
                T = self.value*2**(-self.prec)+other
                return fxp( T, self.prec )
        elif isinstance(other, type(self)):
                assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
                T = self.value+other.value
                return fxp( T, self.prec, format=False )
        raise NotImplementedError(f"No output in add: {type(other)}")

    def __radd__(self,other):
        return self.__add__(other)

    def __sub__(self,other):
        if isinstance(other, alltypes) and not isinstance(other, type(self)):
                T = self.value*2**(-self.prec)-other
                return fxp( T, self.prec )
        elif isinstance(other, type(self)):
                assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
                T = self.value-other.value
                return fxp( T, self.prec, format=False )
        raise NotImplementedError(f"No output in add: {type(other)}")

    def __rsub__(self,other):
        if isinstance(other, alltypes) and not isinstance(other, type(self)):
                T = other - self.value*2**(-self.prec)
                return fxp( T, self.prec )
        raise NotImplementedError(f"No output in add: {type(other)}")

    def __mul__(self,other):
        if isinstance(other, alltypes) and not isinstance(other, type(self)):
                T = self.value*other*2**(-self.prec)
                return fxp( T, self.prec )
        elif isinstance(other, type(self)):
                assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec} {self,other}"
                T = self.value*other.value*2**(-2*self.prec)
                return fxp( T, self.prec )
        raise NotImplementedError(f"No output in lt: {type(other)}")

    def __rmul__(self,other):
        return self.__mul__(other)

    def __lshift__(self, other):
        return fxp( self.value * 2**(other), prec=self.prec )

    def __rshift__(self, other):
        return fxp( self.value *2**(-other), prec=self.prec )

    def __and__(self, other):
        assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
        return fxp( self.value & other.value, prec=self.prec )

    def __or__(self, other):
        assert self.prec == other.prec, f"Precision does not align {self.prec} vs {other.prec}"
        return fxp( self.value | other, prec=self.prec )


    def __repr__(self):
        return f"fxp-{self.prec}: {( self.value/2**self.prec ).n()}"

    def n(self):
        return ( self.value/2**self.prec ).n()

    def floor(self):
        return floor( self.value * 2**(-self.prec) )

    def ceil(self):
        return self.floor()+1

    def round(self):
        return round( self.value / 2**self.prec )

# - - - ISqrt code - - -
def RecSqrt(a, n=0):
    """
    Reciprocal square root algorithm. See Algorithm 5 in: Stan Korzilius and Berry Schoenmakers. Divisions and square roots
    with tight error analysis from newton-raphson iteration in secure fixed-point arithmetic. Cryptogr., 7(3):43, 2023
    """
    # a = fxp( a,prec=pint )
    k = (-2*(floor(binlog(a.round())/2)))
    b = a*2**(k) #scale a. Now  b.n() >= 0.5 and b.n() < 2
    # assert b.n() >= 0.5 and b.n() < 2

    beta = round( RR( ( sqrt(2)-1 )/4 )*2**pint ) / 2**pint
    tau = round( 3/sqrt(2)*2**pint ) / 2**pint
    theta = ceil( log(log(tau*2^(-(pint+n)),tau*beta),2) )  #number of loop's iterations

    c = 3/2 + beta - b*( 2**( 2*round((-binlog(floor(b)))/2)) ).n()  #initial guess
    for i in range(theta):  #Newton-Rapson's method
        z1 = round( 2**(pint+n)*c*b ) / 2**(pint+n)
        z2 = 3 - round( 2**(pint+n)*c*z1  ) / 2**(pint+n)
        c = round( 2**(pint+n)*c*z2 ) / 2**(pint+n+1)
    d = round( 2**(-k/2)*c ) * 2**(k)  #scale back to obtain an answer
    return d

def IntSqrt( a, n=0 ):
    """
    Integer square root algorithm. See Algorithm 7 in: Stan Korzilius and Berry Schoenmakers. Divisions and square roots
    with tight error analysis from newton-raphson iteration in secure fixed-point arithmetic. Cryptogr., 7(3):43, 2023
    """
    a = fxp( a,prec=pint+n )
    reca = floor( a*RecSqrt( a, n=n ) )
    if reca**2 > a:
        return reca - 1
    return reca

# - - - CVP code - - -
def uni(l,r):
    N = r-l+1
    log2_N = binlog( N ) #ceil( log N )
    rnd = 0
    for _ in range(rejection_num):
        """
        Here pint random bits should be sampeled while only log2_N kept.
        We do not do so for acceleration in this particular version of code.
        """
        u = randrange(2**log2_N)
        is_candidate = int(u<N)
        rnd = rnd*(1-is_candidate) + u*is_candidate #constant time if condition
    return l+rnd

def enum_old(B, rad2, target):
    """
    Classical enumeration algorithm. Returns vectors within radius rad2^{1/2} around
    target of vectors from lattice(B).
    param b:       basis matrix in Z^{2 \times 2}
    param rad2:    squared radius of enumeration
    param target:  target vector
    """
    b0, b1 = B[0], B[1]
    r00 = b0.dot_product(b0)
    mu10 = b1.dot_product(b0) / r00
    bstar1 = (b1-mu10*b0)
    r11 = bstar1.dot_product(bstar1)

    out = []
    temp = target * B^-1
    tx, ty = temp[0], temp[1]
    yinter = ( rad2 / r11 ).sqrt()
    print(f"log yinter: {RR(log(yinter,2))}")
    for cy in range( ceil(ty - yinter),floor(ty+yinter)+1 ):
        num = (rad2-r11*(ty-cy)^2)
        if num < 0:
            continue
        xinter = ( num / r00 ).sqrt()
        xcenter = tx + mu10*(ty-cy)
        for cx in range( ceil(xcenter-xinter), floor(xcenter+xinter)+1 ):
            # cand = vector([cx,cy])*B
            # cand = (cand-target)
            # if cand.dot_product(cand) <= rad2:
            out.append( vector([cx,cy]) )
    return out

def enum(B, rad2, target, stepnum=stepnum):
    """
    CVP sampling within radius rad2 around target of vectors from lattice(B).
    param b:       basis matrix in Z^{2 \times 2}
    param rad2:    squared radius of enumeration
    param target:  target vector
    """
    scale = 2**4
    sqrtscale = sqrt(scale)

    #Gram-Schmidt computations.
    b0, b1 = B[0], B[1]
    r00 = b0.dot_product(b0)
    mu10 = b1.dot_product(b0) / r00
    bstar1 = (b1-mu10*b0)
    r11 = bstar1.dot_product(bstar1)

    temp = target * B^-1  #get coefficients of target w.r.t. the basis
    tx, ty = temp[0], temp[1]
    yinter = IntSqrt(round((rad2 / r11)*scale)) / (sqrtscale)  #obtain bounds on the second coefficient's interval length
    center2 = round(target.dot_product(bstar1) / r11)  #obtain second coefficient's interval center
    l, r = floor(center2-yinter)+1, ceil(center2+yinter)
    # print(f"log yinter: {RR(log(yinter,2))}")
    out=[]
    for counter in range(stepnum):  #sample stepnum vectors
        x2 = randrange(l,r) #sample second coefficient
        num = max(0, rad2-r11*(x2-ty)^2)
        xinter = ( IntSqrt(round((num / r00 )*scale)) ) / (sqrtscale)  #obtain bounds on the first coefficient's interval length
        center1 = tx + mu10*(ty-x2)  #obtain first coefficient's interval center
        lx, rx = floor(center1-xinter), ceil(center1+xinter)
        x1 = randrange(lx,rx) #uni(lx,rx)  #sample first coefficient
        out.append( (vector([x1,x2])) )
    return out

if __name__ == "__main__":
    save_plots = True #change to False if no graphs dumps needed
    N = next_prime(ceil(2**144.02))
    C = randrange(1,N)
    rad2 = ceil( 1640*N )
    # rad2 = ceil( 5.5*N**1.29 )

    B = matrix(ZZ,[
        [N, 0],
        [C,1],
    ]).LLL()

    # target = vector( [ round( uniform(C**0.5/4,C**0.5/4) ) for j in range(2) ] )
    nrm = ceil( norm(B[0]).n() )
    answer = vector(ZZ, [ round( uniform(-20,21) ) for j in range(2) ] )

    err = vector( [ round( uniform(-nrm*0.14,nrm*0.14) ) for j in range(2) ] )
    shvec = answer*B
    target = shvec + err

    nerr = norm(err).n()
    print( rad2**0.5 )
    print(nerr)
    print(answer)
    norm(B[0]).n(), norm(B[1]).n()

    print( rad2, norm(B[0].n())**2, norm(B[1]-(B[1].dot_product(B[0]))/(B[0].dot_product(B[0]))*B[0]).n()**2 )
    rad_per_se =  rad2
    print(f"Begining exaustive enumeration...")
    out = enum_old(B,rad_per_se,target)

    l = []
    overshoot = []
    rad = rad2**0.5
    for c in out:
        nrm = norm(c*B-target).n()
        l.append(nrm/rad2**0.5)
        overshoot.append(nrm/nerr)
    print( f"Vectors found: {len(l)}, max ||v||/gh: {max(l)}" )
    H1 = histogram( l, bins=80 )
    H1.show( title="Classic enumeration" )
    if save_plots:
        H1.save_image(  "enum_histogram.png", title="Classic enumeration" )

    xmax, ymax = 0, 0
    for o in out:
        xx, yy = abs(o[0]), abs(o[1])
        if xmax<xx:
            xmax=xx
        if ymax<yy:
            ymax = yy

    outs_classic = set([tuple(o) for o in out])
    allvect = len(out)
    # print(len(out), len(outs_classic))

    xmax = xmax if xmax%2 == 1 else xmax+1
    ymax = ymax if ymax%2 == 1 else ymax+1
    M = [[0 for j in range(-ymax,ymax+1)] for i in range(-xmax,xmax+1)]

    for o in out:
        M[o[0]+xmax][o[1]+ymax] += 1

    Pold = matrix_plot(M)
    Pold.show( figsize=12, title="Exaustive enumeration output" )
    if save_plots:
        Pold.save_image( "enum_domain.png", figsize=12, title="Classic enumeration domain" )

    rad_per_se =  rad2 #min( rad2, N**2/4 ).n()
    print( "Begining CVP sampling" )
    out = enum(B,rad_per_se,target)

    l = []
    overshoot = []
    rad = rad2**0.5
    for c in out:
        nrm = norm(c*B-target).n()
        l.append(nrm/rad2**0.5)
        overshoot.append(nrm/nerr)
    # print( len(l), max(l) )
    H1 = histogram( l, bins=80 )
    H1.show( title="CVP sampling" )
    if save_plots:
        H1.save_image(  "prob_enum_histogram.png", title="CVP sampling" )

    xmax, ymax = 0, 0
    for o in out:
        xx, yy = abs(o[0]), abs(o[1])
        if xmax<xx:
            xmax=xx
        if ymax<yy:
            ymax = yy

    outs_new = set([tuple(o) for o in out])
    print( f"{len(outs_new)} out of requested {stepnum} unique vectors found. max ||v||/gh: {max(l)}" )
    print(f"Saturation percentage: {(len(outs_new)/stepnum).n()}")  #why is this one is close to 1 / golden_ratio ?

    xmax = xmax if xmax%2 == 1 else xmax+1
    ymax = ymax if ymax%2 == 1 else ymax+1
    M = [[0 for j in range(-ymax,ymax+1)] for i in range(-xmax,xmax+1)]

    for o in out:
        M[o[0]+xmax][o[1]+ymax] += 1

    Pnew = matrix_plot(M)
    Pnew.show( figsize=12 )
    if save_plots:
        Pnew.save_image( "prob_enum_domain.png", figsize=12, title="CVP sampler domain" )
