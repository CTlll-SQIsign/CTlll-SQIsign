class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

RR = RealField( 172 )
# scale_factor = 180
lagrange_max_steps = 10
lll_max_tours = 18

from fpylll import *
try:
    from multiprocess import Pool  # you might need pip install multiprocess
except ModuleNotFoundError:
    from multiprocessing import Pool

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

RR = RealField( 192 )
RealNumber = RR
lagrange_max_steps = 48
lll_max_tours = 94

from fpylll import *

def bilin( u,v,p=1 ):
    return sum( u[i]*v[i] for i in range(len(u)//2) ) + p*sum( u[i]*v[i] for i in range(len(u)//2,len(u)) )

class GSO_obj: #checked

    def __init__(self,B,p=1):
        """
        GSO_obj constructor.
        param B: list of vectors over QQ - module basis.
        """

        n = len(list(B)) #lattice rank
        m = len(B[0])
        self.n, self.m = n,m
        self.G = [ [None for j in range(i+1)] for i in range(n) ]
        self.p = p

        for i in range(n):
            for j in range(i+1):
                self.G[i][j] =  bilin(B[i],B[j],p) #B[i].dot_product(B[j])

        self.Mu, self.Rr = [[None for j in range(i)] for i in range(n)], [[None for j in range(i+1)] for i in range(n)]

    def __len__(self):
        return len(self.G)

    def __str__(self):
        s=""
        for i in range(self.n):
            for j in range(i+1):
                s+= str(self.G[i][j]) + ", "
            if i!= len(self)-1: s+="\n"
        return(s)

    def __repr__(self):
        s=""
        for i in range(self.n):
            for j in range(i+1):
                print(i,j)
                s+= str(self.G[i][j]) + ", "
            if i!= len(self)-1: s+="\n"
        return(s)


    def update_G(self, B, start):
        """
        Updates G starting at start.
        param B: list of nf_vect - module basis.
        param G: position <n.
        """
        for i in range(start, self.n):
            for j in range(i+1):
                self.G[i][j] = bilin(B[i],B[j],p) #B[i].dot_product(B[j])

    def compute_GSO(self,start=0, end=None):
        """
        Modifies GSO coefficients Mu and Rr at positions start,...,end-1.
        param start: position to start modifying Mu and Rr
        param end: position to end modifying Mu and Rr

        The algorithm is an adaptation of [https://perso.ens-lyon.fr/damien.stehle/downloads/fpLLL_journal.pdf].
        Notice that <au,bv> is a*b*<u,v>.
        """
        n = self.n
        if end is None:
            end=n
        assert end<=n, "Wrong dimensions!"
        assert start<end, "start>=end!"

        for i in range(start,end):
            for j in range(i):
                self.Rr[i][j] = self.G[i][j]
                for k in range(j):
                    self.Rr[i][j] -= self.Mu[j][k]*self.Rr[i][k]
                self.Mu[i][j] = self.Rr[i][j] / self.Rr[j][j]
            self.Rr[i][i]= self.G[i][i]
            for j in range(i):
                self.Rr[i][i]-=self.Mu[i][j]*self.Rr[i][j]

    def size_reduce(self,B, start=0, end=0, debug=0):
        """
          This function returns size reduced B_start,...,B_{end-1} where B is the basis of lattice.
          param B: basis to be updated
          param start: position to start at.
          params end: position to end before.
          Updates U, self and B s.t. U*B entries starting at start and ending at end-1 are size reduced and self is valid from 0 to end-1.
        """
        n = self.n
        for i in range(start,end):
            for j in range(i-1,-1,-1):
                self.compute_GSO( start=0, end=i+1 )
                delta_nf = round( self.Mu[i][j] ) #we will make transformation b_i <- b_i-delta_nf*b_j

                gij, gjj = self.G[i][j], self.G[j][j]
                for kappa in range(i):
                     self.G[i][kappa] -=  delta_nf*self.G[j][kappa] if kappa<=j else delta_nf*self.G[kappa][j]  #<b_i-delta_nf*b_j,b_kappa> = G_i_kappa - delta_nf*G_j_kappa
                for kappa in range(i+1,n):
                    self.G[kappa][i] -=  delta_nf*self.G[kappa][j] #<b_kappa,b_i-delta_nf*b_j> = G_kappa_i - delta_nf*G_kappa_j
                self.G[i][i] -= ( 2*delta_nf*gij-delta_nf*delta_nf*gjj ) #<b_i-delta_nf*b_j,b_i-delta_nf*b_j> = G_i_i -2*delta_nf*G_i_j + delta_nf^2*G_j_j
                B[i] -= delta_nf *B[j] #update the basis
                self.compute_GSO( start=i, end=i+1 )
                self.check_integrity( B,0,self.n )
        return B

    def update_after_svp_oracle_rank_2(self,U_,i):
        """
        After SVP oracle in lll is done, updates self.
        param U_: transformation matrix (list of nf_vect) returned by SVP oracle
        param i: position at which the SVP oracle has been called
        """
        n=self.n
        u0,u1 = U_
        g0, g1, g2 = self.G[i][i], self.G[i+1][i+1], self.G[i+1][i]

        """
        For the both rows (zone I) we have that:
            <b'[i],b'[ell]> = U_[0][0]*<b[i],b[ell]> + U_[0][1]*<b[i+1],b[ell]>
            <b'[i+1],b'[ell]> = U_[1][0]*<b[i],b[ell]> + U_[1][1]*<b[i+1],b[ell]>
        """
        for ell in range(i):
            tmp = U_[0][0]*self.G[i][ell] + U_[0][1]*self.G[i+1][ell]
            self.G[i+1][ell] = U_[1][0]*self.G[i][ell] + U_[1][1]*self.G[i+1][ell]
            self.G[i][ell] = tmp

        """
        For the both columns (zone III) we have that:
            <b'[kappa],b'[i]>   = U_[0][0]*<b[kappa],b[i]> + U_[0][1]*<b[kappa],b[i+1]>
            <b'[kappa],b'[i+1]> = U_[1][0]*<b[kappa],b[i]> + U_[1][1]*<b[kappa],b[i+1]>
        """
        for kappa in range(i+2,n):
            tmp =                U_[0][0]*self.G[kappa][i] + U_[0][1]*self.G[kappa][i+1]
            self.G[kappa][i+1] = U_[1][0]*self.G[kappa][i] + U_[1][1]*self.G[kappa][i+1]
            self.G[kappa][i] = tmp

        """
        the rest is dealt with with similar formulas
        self.G[i][i] =      u0[0]*u0[0]*g0 + u0[0]*u0[1]*g2+u0[1]*u0[0]*g2+u0[1]*u0[1]*g1
        self.G[i+1][i+1] =  u1[0]*u1[0]*g0 + u1[0]*u1[1]*g2+u1[1]*u1[0]*g2+u1[1]*u1[1]*g1
        self.G[i+1][i] =    u1[0]*u0[0]*g0 + u1[0]*u0[1]*g2+u1[1]*u0[0]*g2+u1[1]*u0[1]*g1
        """
        self.G[i][i] =      u0[0]*u0[0]*g0 + (u0[0]*u0[1]+u0[1]*u0[0])*g2+u0[1]*u0[1]*g1
        self.G[i+1][i+1] =  u1[0]*u1[0]*g0 + (u1[0]*u1[1]+u1[1]*u1[0])*g2+u1[1]*u1[1]*g1
        self.G[i+1][i] =    u1[0]*u0[0]*g0 + (u1[0]*u0[1]+u1[1]*u0[0])*g2+u1[1]*u0[1]*g1

    def check_consistency(self, B, start=0, end=1):
        """
        Checks:
        1) If Mu factor of B is equal to self.Mu
        2) If Rr factor of B is equal to self.Rr
        3) If Mu[i][j] == Rr[i][j] / Rr[j][j]
        param B: up-to-date basis of module (consists of nf_vect)
        param start: position to start at.
        params end: position to end before.
        """
        n = len(self.G)
        m = n
        self.compute_GSO(start=0,end=len(self))
        G = GSO_obj(B,p=self.p)
        G.compute_GSO(start=0, end=end)

        Mu, Rr = G.Mu, G.Rr
        Mu_, Rr_ = self.Mu, self.Rr
        for i in range(end):
            for j in range(i):
                if abs( Mu[i][j] - Mu_[i][j] ) > 10**-3:
                    print(f"{bcolors.WARNING}Non consistent GSO at Mu_{i,j}:{bcolors.ENDC} {abs( Mu[i][j] - Mu_[i][j] ).n() }")

        for i in range(start, end):
            for j in range(i+1):
                if Rr[i][j] is None:
                    continue
                else:
                    if abs( Rr[i][j] - Rr_[i][j] ) > 10**-3:
                        print(f"{bcolors.WARNING}Non consistent GSO at Rr_{i,j}:{bcolors.ENDC} {abs( Rr[i][j] - Rr_[i][j] ).n() }")

        for i in range(start, end):
            for j in range(i):
                tmp0 = Rr[i][j] / Rr[j][j]
                tmp =  Mu[i][j]
                if abs(tmp-tmp0)>10**-3:
                    print(f"{bcolors.WARNING}Warning:{bcolors.ENDC} R/R={abs(ln(tmp))} at {i,j}")

    def check_integrity( self,B,start=0,end=0 ):
        """
        Checks if Gram matrix of self is equal to the Gram matrix of B
        param B: up-to-date basis of module (consists of nf_vect)
        param start: position to start at.
        params end: position to end before.
        """
        n = len(self.G)
        H = GSO_obj(B,p=self.p)
        for i in range(start,max(end,n)):
            for j in range(len(self.G[i])):
                diff = self.G[i][j] - H.G[i][j]
                if  abs(diff) > 10**-3:
                    print(f"{bcolors.WARNING}Troubles at {i,j},{bcolors.ENDC} delta={diff.n(50)}")

    def is_lll_reduced( self,delta,eta ):
        for i in range(1,self.n):
            for j in range(i):
                if abs( self.Mu[i][j] ) > eta:
                    print("Bad size!")
                    return False
        for i in range(self.n-1):
            if delta*self.Rr[i][i] > self.Mu[i+1][i]**2 * self.Rr[i][i] + self.Rr[i+1][i+1]:
                print(f"Bad lovacz at {i,j}, {delta*self.Rr[i][i].n() , self.Mu[i+1][i]**2 * self.Rr[i][i] + self.Rr[i+1][i+1].n()}")
                return False
        return True

def lagrange_reduction_gram(G):
    """
    Finds a transformation U such that U*B is a Lagrange-reduced basis of a 2-dim lattice.
    param G: Gram matrix B*B.transpose().
    """
    mu = G[1][0] / G[0][0]
    #size reduction
    U = matrix.identity(2)

    for cntr in range(lagrange_max_steps):
        mu = G[1][0] / G[0][0]
        #size reduction
        muround = round( mu )
        G[1][1] -= ( 2*muround*G[1][0]-muround**2 * G[0][0] )
        G[1][0] -= muround*G[0][0]
        mu -= muround
        U[1] -= muround*(U[0])

        G[0][0], G[1][1] = G[1][1], G[0][0] #swap
        U[0], U[1] = U[1], U[0]
    mu = G[1][0] / G[0][0]
    condition = (G[1][1]) < G[0][0]
    U[0], U[1] = (U[1], U[0]) if condition else (U[0], U[1])  #alias for the constant time condition*U[0] + (1-condition)*U[1] trick
    return(U)

def bitsize(n):
    # given an integer or rational n, returns the bitsize of integer(s) it consists of
    if isinstance(n,(type(ZZ),type(0))):
        return( len(bin(n))-2 )
    if isinstance( n, type(1/2) ):
        return( max( bitsize(n.numerator()) , bitsize(n.denominator()) ) )

def maxbitsize(G):
    #given GSO object, returns the max bitsize of integers stored within
    mx = 0
    for gg in G.G:
        for g in gg:
            s = bitsize(g)
            if s > mx:
                mx = s
    # print(f"G max: {mx}")

    mu = 0
    for gg in G.Mu:
        for g in gg:
            s = bitsize(g)
            if s > mx:
                mu = s
    # print(f"Mu max: {mu}")

    rr = 0
    for gg in G.Mu:
        for g in gg:
            s = bitsize(g)
            if s > rr:
                rr = s
    # print(f"RR max: {rr}")
    return max([mx,mu,rr])

def bitsize_bnd(B,p):
    #given matrix B and p>0 \in ZZ returns the estimated bitsize required to run BKZ-2
    n = B.nrows()
    C = B*diagonal_matrix( [1]*(n//2) + [sqrt(p)]*(n//2) )
    maxnrm = max( norm(c).n() for c in C )**2
    return n + ceil(3*(n-1)*log( maxnrm,2 )/2)  + log(n,2)/2 + 1

def lll_const( B, p=1, logger=None ):
    """
    Applies constant time LLL reduction to B.
    param B: a square nonsingular matrix over QQ.
    """
    G = GSO_obj( B, p=p )
    G.compute_GSO(0,B.nrows())
    n = len(list(B))
    maxbit = 0
    for t in range(lll_max_tours):
        curbit = maxbitsize(G)
        if curbit > maxbit:
            maxbit = curbit
        # assert bitsize_bnd(B,p) >= curbit, f"Bug: {bitsize_bnd(B,p)} < { curbit }"
        # a0, a1 =  bitsize_bnd(B,p) , maxbitsize(G)
        for i in range(n-1):
            B = G.size_reduce(B,i+1,i+2)

            M =  [ [G.Rr[i][i]], [G.Mu[i+1][i]*G.Rr[i][i],G.Mu[i+1][i]**2*G.Rr[i][i] + G.Rr[i+1][i+1]] ]
            detM = G.Rr[i][i]*(G.Mu[i+1][i]**2*G.Rr[i][i] + G.Rr[i+1][i+1]) - (G.Mu[i+1][i]*G.Rr[i][i])**2 #debuf for heur
            assert detM >= 0, f"AAA: {detM}"

            U = lagrange_reduction_gram(M)
            B[i:i+2]=U*B[i:i+2]
            G.update_after_svp_oracle_rank_2(U,i)
            G.compute_GSO(start=i, end=i+2)

            B = G.size_reduce(B,i,i+2)
            # G.check_consistency(B,0,i+2)
            # G.check_integrity(B,0,i+2)
        if not logger is None:
            # logger.append((G.Rr[i][i] / detM)**0.5)
            logger.append(maxbit)
    # print(f"Is LLL red: {G.is_lll_reduced(0.99,0.501)}")
    return B, maxbit

def test_subroutines():
    #Correctness of Cholesky
    n=10 #dimension 10

    M = matrix([[randrange(-102,103)*(i+1) for j in range(n)] for i in range(n)])

    G = GSO_obj( M )
    G.compute_GSO()

    # print( G )
    Gfpylll = GSO.Mat( IntegerMatrix.from_matrix(M) )
    Gfpylll.update_gso()

    for i in range(n-1):
        assert abs( G.Rr[i][i].n() - Gfpylll.get_r(i,i) ) < 10**-6, f"Cholesky inrorrect for Rr[{i}][{i}]"
        for j in range(i):
            assert abs( G.Mu[i][j].n() - Gfpylll.get_mu(i,j) ) < 10**-6, f"Cholesky inrorrect for Mu[{i}][{j}]"
    print(f"Cholesky correct!")

    #Size reducedness
    n=10 #dimension 10
    M = matrix([[randrange(-102,103)*(i+1) for j in range(n)] for i in range(n)])
    G = GSO_obj( M )
    G.compute_GSO()
    G.size_reduce(M,start=0,end=n)
    for i in range(n):
        for j in range(i):
            assert( abs(G.Mu[i][j]) <0.501 ), f"Size reduction fail at Mu[{i,j}]"
    print("Size reduction test successful!")


def tournum( B,p ):
    """
    Returns a number of tours predicted using the Sandpile Model.
    param B: 4 x 4 basis of a lattice.
    """
    G = GSO_obj( B,p=p )
    G.compute_GSO()

    r2 = [ G.Rr[i][i].n() for i in range(len(G.Rr)) ]
    r = vector( [ (log(rr,2)/2).n() for rr in r2 ] )
    print(r)
    scale = sum(r)/4
    r = vector( [ rr-scale for rr in r ] )

    nu2 = log( 4/3, 2 )/2
    hkzprof = vector( [ 2*nu2, nu2, -nu2, -2*nu2 ] )

    nrm = norm(hkzprof-r)
    # print(f"nrmdeb: {nrm.n()}")
    try:
        return ceil( RR( 2*log(nrm) / log(8/7) ) )
    except:
        return 0
