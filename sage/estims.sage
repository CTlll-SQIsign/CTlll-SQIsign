load("clll.sage")

lagrange_max_steps = 15
lll_max_tours = 15

logger = []
lllred = []
ntests = 200
maxreq, maxbit = 0, 0

p=127
log_bound_on_coeffs = 320.  #how large the coefficients of matrices can be
maxbit = 0
maxcurbit = 0
for tests in range(ntests):  #we ran ntests tests and store maximal bitsize required in maxbit
    B = IntegerMatrix(4,4)
    B.randomize("qary", bits=log_bound_on_coeffs, k=2)
    if tests==0:
        tn = tournum(B,p)
    bt  = bitsize_bnd(matrix(B),p)
    if maxbit < bt:
        maxbit = bt
    elif tests%200==0:
        print(f"{tests} out of {ntests} done")

    C = matrix(B)
    C, curbit = lll_const(C, p=p, logger=None)
    print(f"tournum: {lll_max_tours} vs required {tn} | bit: {curbit} vs required {bt}")
    if curbit>maxcurbit:
        maxcurbit=curbit
        assert curbit <= bt
    nrmmy = [ norm(cc).n() for cc in C ]
    G = GSO.Mat( B, float_type="dd")
    lll = LLL.Reduction( G )
    lll()
    B = matrix( lll.M.B )
    nrmfpylll = [ norm(bb).n() for bb in B ]
    # print(f"fpylll: {nrmfpylll}")
    # print(f"My: {nrmmy}")

    G = GSO_obj( C, p=p )
    G.compute_GSO(0,C.nrows())
    isred = G.is_lll_reduced(0.95,0.55)
    print(f"Approx_fact vs Euclid: {min(nrmmy) / min(nrmfpylll)}. Is LLL reduced: {isred}")
    logger.append( min(nrmmy) / min(nrmfpylll) )
    assert isred, f"NOT RED!"

print( f"Required precision: {maxbit} vs observed precision: {maxcurbit} | required tour number: {tn}" )
