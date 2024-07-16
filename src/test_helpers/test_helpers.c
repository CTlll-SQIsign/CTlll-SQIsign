#include "test_helpers.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void nct_ibz_xgcd_non_ct(nct_ibz_t *gcd, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *n, const nct_ibz_t *d){
    nct_ibz_t r0, r1, u0,u1, v0,v1, q,r;
    int swaped = 0;
    nct_ibz_init(&r0);
    nct_ibz_init(&r1);
    nct_ibz_init(&u0);
    nct_ibz_init(&u1);
    nct_ibz_init(&v0);
    nct_ibz_init(&v1);
    nct_ibz_init(&q);
    nct_ibz_init(&r);
    nct_ibz_abs(&r0,n);
    nct_ibz_abs(&r1,d);
    if(nct_ibz_cmp(&r1,&r0)>0){
        nct_ibz_swap(&r1,&r0);
        swaped = 1;
    }
    nct_ibz_set(&u0,1);
    nct_ibz_set(&u1,0);
    nct_ibz_set(&v0,0);
    nct_ibz_set(&v1,1);
    while(!nct_ibz_is_zero(&r1)){
        nct_ibz_div_floor(&q,&r0,&r0,&r1);
        nct_ibz_swap(&r0,&r1);
        nct_ibz_mul(&r,&u1,&q);
        nct_ibz_sub(&u0,&u0,&r);
        nct_ibz_swap(&u0,&u1);
        nct_ibz_mul(&r,&v1,&q);
        nct_ibz_sub(&v0,&v0,&r);
        nct_ibz_swap(&v0,&v1);
    }
    nct_ibz_copy(gcd,&r0);
    if(swaped){
        nct_ibz_swap(&u0,&v0);
    }
    nct_ibz_set(&r0,0);
    if(nct_ibz_cmp(n,&r0)<0){
        nct_ibz_neg(&u0,&u0);
    }
    if(nct_ibz_cmp(d,&r0)<0){
        nct_ibz_neg(&v0,&v0);
    }
    nct_ibz_copy(u,&u0);
    nct_ibz_copy(v,&v0);
    nct_ibz_finalize(&r0);
    nct_ibz_finalize(&r1);
    nct_ibz_finalize(&u0);
    nct_ibz_finalize(&u1);
    nct_ibz_finalize(&v0);
    nct_ibz_finalize(&v1);
    nct_ibz_finalize(&q);
    nct_ibz_finalize(&r);
}

void nct_ibz_div_towards_zero(nct_ibz_t *q, nct_ibz_t *r, const nct_ibz_t *n, const nct_ibz_t *d){
    nct_ibz_t prod, qx, rx, abs_n, abs_d,zero;
    nct_ibz_init(&prod);
    nct_ibz_init(&abs_n);
    nct_ibz_init(&abs_d);
    nct_ibz_init(&qx);
    nct_ibz_init(&rx);
    nct_ibz_init(&zero);
    nct_ibz_set(&zero,0);
    nct_ibz_abs(&abs_n,n);
    nct_ibz_abs(&abs_d,d);
    nct_ibz_div_floor(&qx,&rx,&abs_n,&abs_d);
    nct_ibz_mul(&prod,n,d);
    int sign_q = (nct_ibz_cmp(&prod,&zero)>0);
    int sign_r = (nct_ibz_cmp(n,&zero)>0);
    nct_ibz_neg(&prod,&qx);
    nct_ibz_conditional_assign(q,&qx,&prod,sign_q);
    nct_ibz_neg(&prod,&rx);
    nct_ibz_conditional_assign(r,&rx,&prod,sign_r);
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&abs_n);
    nct_ibz_finalize(&abs_d);
    nct_ibz_finalize(&qx);
    nct_ibz_finalize(&rx);
    nct_ibz_finalize(&zero);
}

int nct_ibz_vec_4_scalar_div(nct_ibz_vec_4_t *quot, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec){
    int res = 1;
    nct_ibz_t r;
    nct_ibz_init(&r);
    for(int i = 0; i < 4; i++){
        nct_ibz_div_towards_zero(&((*quot)[i]),&r,&((*vec)[i]),scalar);
        res = res && nct_ibz_is_zero(&r);
    }
    nct_ibz_finalize(&r);
    return(res);
}

int nct_ibz_mat_4x4_scalar_div(nct_ibz_mat_4x4_t *quot, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat){
    int res = 1;
    nct_ibz_t r;
    nct_ibz_init(&r);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_div_towards_zero(&((*quot)[i][j]),&r,&((*mat)[i][j]),scalar);
            res = res && nct_ibz_is_zero(&r);
        }
    }
    nct_ibz_finalize(&r);
    return(res);
}

int nct_quat_test_xgcd_verify(const nct_ibz_t *gcd, const nct_ibz_t *u, const nct_ibz_t *v, const nct_ibz_t *d, const nct_ibz_t *n){
    int res = 0;
    nct_ibz_t q,r, prod, tmp;
    nct_ibz_init(&q);
    nct_ibz_init(&r);
    nct_ibz_init(&prod);
    nct_ibz_init(&tmp);
    nct_ibz_mul(&prod,u,n);
    nct_ibz_mul(&tmp,v,d);
    nct_ibz_add(&tmp,&tmp,&prod);
    res = res || (0!=nct_ibz_cmp(&tmp,gcd));
    nct_ibz_div_floor(&q,&r,d,gcd);
    res = res || !(nct_ibz_is_zero(&r));
    nct_ibz_div_floor(&q,&r,n,gcd);
    res = res || !(nct_ibz_is_zero(&r));
    nct_ibz_finalize(&q);
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&tmp);
    return(res);
}





//copied from SQIsign

//Small helper for integers
void nct_ibz_mod_not_zero(nct_ibz_t *res, const nct_ibz_t *x, const nct_ibz_t *mod){
    nct_ibz_t m,t;
    nct_ibz_init(&m);
    nct_ibz_init(&t);
    nct_ibz_mod(&m,x,mod);
    nct_ibz_set(&t, nct_ibz_is_zero(&m));
    nct_ibz_mul(&t,&t,mod);
    nct_ibz_add(res,&m,&t);
    nct_ibz_finalize(&m);
    nct_ibz_finalize(&t);
}

//centered and rather positive then negative
void nct_ibz_centered_mod(nct_ibz_t *remainder, const nct_ibz_t *a, const nct_ibz_t *mod){
    //assert(nct_ibz_cmp(mod,&nct_ibz_const_zero)>0);
    nct_ibz_t tmp, d, t;
    nct_ibz_init(&tmp);
    nct_ibz_init(&d);
    nct_ibz_init(&t);
    nct_ibz_set(&t,2);
    nct_ibz_div_floor(&d,&tmp,mod,&t);
    nct_ibz_mod_not_zero(&tmp,a,mod);
    nct_ibz_set(&t,nct_ibz_cmp(&tmp,&d)>0);
    nct_ibz_mul(&t,&t,mod);
    nct_ibz_sub(remainder,&tmp,&t);
    nct_ibz_finalize(&tmp);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&t);
}

// if c, res = x, else res = y
void nct_ibz_conditional_assign(nct_ibz_t *res, const nct_ibz_t *x, const nct_ibz_t *y, int c){
    nct_ibz_t s, t, r;
    nct_ibz_init(&r);
    nct_ibz_init(&s);
    nct_ibz_init(&t);
    nct_ibz_set(&s,c!=0);
    nct_ibz_set(&t,1);
    nct_ibz_sub(&t,&t,&s);
    nct_ibz_mul(&r,&s,x);
    nct_ibz_mul(res,&t,y);
    nct_ibz_add(res,&r,res);
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&s);
    nct_ibz_finalize(&t);
}

void nct_ibz_vec_4_linear_combination(nct_ibz_vec_4_t *lc, const nct_ibz_t *coeff_a, const nct_ibz_vec_4_t  *vec_a, const nct_ibz_t *coeff_b, const nct_ibz_vec_4_t *vec_b){
    nct_ibz_t prod;
    nct_ibz_vec_4_t sums;
    nct_ibz_vec_4_init(&sums);
    nct_ibz_init(&prod);
    for (int i = 0; i <4; i++){
        nct_ibz_mul(&(sums[i]),coeff_a,&((*vec_a)[i]));
        nct_ibz_mul(&prod,coeff_b,&((*vec_b)[i]));
        nct_ibz_add(&(sums[i]),&(sums[i]),&prod);
    }
    for (int i = 0; i <4; i++){
        nct_ibz_copy(&((*lc)[i]),&(sums[i]));
    }
    nct_ibz_finalize(&prod);
    nct_ibz_vec_4_finalize(&sums);
}

//understand how to put the non-zero on the 1st one
void nct_ibz_xgcd_with_u_not_0(nct_ibz_t *d, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *x, const nct_ibz_t *y){
    nct_ibz_t q, x1,y1, t, old_x, old_y, constants; 
    int c;
    nct_ibz_init(&constants);
    nct_ibz_init(&q);
    nct_ibz_init(&t);
    nct_ibz_init(&x1);
    nct_ibz_init(&y1);
    nct_ibz_init(&old_x);
    nct_ibz_init(&old_y);
    nct_ibz_copy(&x1,x);
    nct_ibz_copy(&y1,y);
    nct_ibz_copy(&old_x,x);
    nct_ibz_copy(&old_y,y);
    //prepare inputs:
    //set x to 1 if both 0
    c = nct_ibz_is_zero(&x1) & nct_ibz_is_zero(&y1);
    nct_ibz_set(&t, c);
    nct_ibz_add(&x1,&x1,&t);
    //xgcd
    nct_ibz_xgcd_non_ct(d,u,v,&x1,&y1);
    //make sure u!=0 (v can be 0 if needed)
    nct_ibz_copy(&x1,&old_x);
    c = nct_ibz_is_zero(&y1);
    nct_ibz_set(&t, c);
    nct_ibz_add(&y1,&y1,&t);
    nct_ibz_div_towards_zero(&q,&t,&x1,&y1);

    c = nct_ibz_is_zero(u);
    nct_ibz_set(&t, c);
    nct_ibz_add(u,u,&t);//if u is 0, set u to 1 else leave it u -> u+c1
    nct_ibz_set(&constants,1);
    nct_ibz_sub(&t,&constants,&q);
    nct_ibz_conditional_assign(&t,v,&t,nct_ibz_is_zero(&old_x));//if x is zero, set t to v (v unchanged)
    nct_ibz_copy(&q,v);
    nct_ibz_conditional_assign(v,&t,&q,c);//prepare u-1 put it into u , else not: 1-q

    //try to match sign rules
    nct_ibz_set(&constants,0);
    c = (nct_ibz_cmp(d,&constants))<0;
    nct_ibz_set(&t,c);
    nct_ibz_set(&constants,2);
    nct_ibz_mul(&t,&t,&constants);
    nct_ibz_set(&constants,1);
    nct_ibz_sub(&t,&constants,&t);
    nct_ibz_mul(u,u,&t);
    nct_ibz_mul(d,d,&t);
    nct_ibz_mul(v,v,&t);

    nct_ibz_finalize(&x1);
    nct_ibz_finalize(&y1);
    nct_ibz_finalize(&old_x);
    nct_ibz_finalize(&old_y);
    nct_ibz_finalize(&q);
    nct_ibz_finalize(&t);
    nct_ibz_finalize(&constants);
}

void nct_ibz_rounded_div(nct_ibz_t *q, const nct_ibz_t *a, const nct_ibz_t *b){
    nct_ibz_t r,sign_q, abs_b, zero;
    nct_ibz_init(&r);
    nct_ibz_init(&sign_q);
    nct_ibz_init(&abs_b);
    nct_ibz_init(&zero);
    nct_ibz_set(&zero,0);

    //assumed to round towards 0
    nct_ibz_abs(&abs_b,b);
    // q is of same sign as a*b (and 0 if a is 0)
    nct_ibz_mul(&sign_q,a,b);
    nct_ibz_div_towards_zero(q,&r,a,b);
    nct_ibz_abs(&r,&r);
    nct_ibz_add(&r,&r,&r);
    nct_ibz_set(&sign_q,(1-2*(nct_ibz_cmp(&sign_q,&zero)<0))*(nct_ibz_cmp(&r,&abs_b)>0));
    nct_ibz_add(q,q,&sign_q);
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&sign_q);
    nct_ibz_finalize(&abs_b);
    nct_ibz_finalize(&zero);
}

// HNF functions
int nct_ibz_mat_4x4_is_hnf(const nct_ibz_mat_4x4_t *mat){
    int res = 1;
    int found = 0;
    int ind = 0;
    nct_ibz_t zero;
    nct_ibz_init(&zero);
    // upper triangular
    for (int i = 0; i < 4; i++){
        // upper triangular
        for (int j = 0; j < i; j++){
            res = res && nct_ibz_is_zero(&((*mat)[i][j]));
        }
        // find first non 0 element of line
        found = 0;
        for (int j = i; j < 4; j++){
            if(found){
                // all values are positive, and first non-0 is the largest of that line
                res = res && (nct_ibz_cmp(&((*mat)[i][j]),&zero)>=0);
                res = res && (nct_ibz_cmp(&((*mat)[i][ind]),&((*mat)[i][j]))>0);
            } else {
                if(!nct_ibz_is_zero(&((*mat)[i][j]))){
                    found = 1;
                    ind = j;
                    // mustbe non-negative
                    res = res && (nct_ibz_cmp(&((*mat)[i][j]),&zero)>0);
                }
            }
        } 
    }
    // check that first nom-zero elements ndex per column is strictly increasing
    int linestart = -1;
    int i = 0;
    for(int j = 0; j<4; j++){
        while((i < 4) &&(nct_ibz_is_zero(&((*mat)[i][j])))){
            i = i+1;
        } if (i != 4) {
            res = res && (linestart < i);
        }
        i = 0;
    }
    nct_ibz_finalize(&zero);
    return res;
}

//centered mod
void nct_ibz_vec_4_linear_combination_mod(nct_ibz_vec_4_t *lc, const nct_ibz_t *coeff_a, const nct_ibz_vec_4_t  *vec_a, const nct_ibz_t *coeff_b, const nct_ibz_vec_4_t *vec_b, const nct_ibz_t *mod){
    nct_ibz_t prod, m;
    nct_ibz_vec_4_t sums;
    nct_ibz_vec_4_init(&sums);
    nct_ibz_init(&prod);
    nct_ibz_init(&m);
    nct_ibz_copy(&m,mod);
    for (int i = 0; i <4; i++){
        nct_ibz_mul(&(sums[i]),coeff_a,&((*vec_a)[i]));
        nct_ibz_mul(&prod,coeff_b,&((*vec_b)[i]));
        nct_ibz_add(&(sums[i]),&(sums[i]),&prod);
        nct_ibz_centered_mod(&(sums[i]),&(sums[i]),&m);
    }
    for (int i = 0; i <4; i++){
        nct_ibz_copy(&((*lc)[i]),&(sums[i]));
    }
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&m);
    nct_ibz_vec_4_finalize(&sums);
}

void nct_ibz_vec_4_copy_mod(nct_ibz_vec_4_t *res, const nct_ibz_vec_4_t *vec, const nct_ibz_t *mod){
    nct_ibz_t m;
    nct_ibz_init(&m);
    nct_ibz_copy(&m,mod);
    for(int i = 0; i < 4; i++){
        nct_ibz_centered_mod(&((*res)[i]),&((*vec)[i]),&m);
    }
    nct_ibz_finalize(&m);
}

// no need to center this, and not 0
void nct_ibz_vec_4_scalar_mul_mod(nct_ibz_vec_4_t *prod, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec, const nct_ibz_t *mod){
    nct_ibz_t m, s;
    nct_ibz_init(&m);
    nct_ibz_init(&s);
    nct_ibz_copy(&s,scalar);
    nct_ibz_copy(&m,mod);
    for(int i = 0; i < 4; i++){
        nct_ibz_mul(&((*prod)[i]),&((*vec)[i]),&s);
        nct_ibz_mod(&((*prod)[i]),&((*prod)[i]),&m);
    }
    nct_ibz_finalize(&m);
    nct_ibz_finalize(&s);
}

//Algorithm used is the one at number 2.4.8 in Henri Cohen's "A Course in Computational Algebraic Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993
// assumes nct_ibz_xgcd outputs u,v which are small in absolute value (as described in the book)
void nct_ibz_mat_4xn_hnf_mod_core(nct_ibz_mat_4x4_t *hnf, int generator_number, const nct_ibz_vec_4_t (*generators)[generator_number], const nct_ibz_t *mod){
    int i,j;
    //assert(generator_number>3);
    int n = generator_number;
    int k = n;
    nct_ibz_t b, u, v, d, q, m, coeff_1, coeff_2, r, one;
    nct_ibz_vec_4_t c;
    nct_ibz_vec_4_t a[generator_number];
    nct_ibz_vec_4_t w[4];
    nct_ibz_init(&one);
    nct_ibz_init(&b);
    nct_ibz_init(&d);
    nct_ibz_init(&u);
    nct_ibz_init(&v);
    nct_ibz_init(&r);
    nct_ibz_init(&m);
    nct_ibz_init(&q);
    nct_ibz_init(&coeff_1);
    nct_ibz_init(&coeff_2);
    nct_ibz_vec_4_init(&c);
    for (int h = 0; h < n; h++){
        if(h<4)
            nct_ibz_vec_4_init(&(w[h]));
        nct_ibz_vec_4_init(&(a[h]));
        nct_ibz_copy(&(a[h][0]), &((*generators)[h][0]));
        nct_ibz_copy(&(a[h][1]), &((*generators)[h][1]));
        nct_ibz_copy(&(a[h][2]), &((*generators)[h][2]));
        nct_ibz_copy(&(a[h][3]), &((*generators)[h][3]));
    }
    nct_ibz_set(&one,1);
    //nct_ibz_set(&d,0);
    //assert(nct_ibz_cmp(mod,&d)>0);
    nct_ibz_copy(&m,mod);
    for (i =3; i >-1;i--){
        k = k-1;
        for (j = k-1; j > -1; j--){
            // assumtion that nct_ibz_xgcd outputs u,v which are small in absolute value is needed here
            // also, needs u non 0, but v can be 0 if needed
            // also requires xgcd to output aki, 1, 0 on input aki,0 (and d=1,u=1,v=0 if aki = 0): (means nothing changes if aji = 0)
            nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&(a[k][i]),&(a[j][i]));
            nct_ibz_vec_4_linear_combination(&c,&u,&(a[k]),&v,&(a[j]));
            nct_ibz_div_towards_zero(&coeff_1,&r, &(a[k][i]),&d);//modified line from nct_ibz_div which might cause serious trouble
            //if coeff1 is 0, set it to 1 to ensure not to loose info on j column
            nct_ibz_set(&r,(nct_ibz_is_zero(&coeff_1) & (nct_ibz_is_zero(&(a[j][i]))))!=0);
            nct_ibz_add(&coeff_1,&r, &coeff_1);
            //compute k coefficient in lc for j
            nct_ibz_div_towards_zero(&coeff_2,&r, &(a[j][i]),&d);//modified line from nct_ibz_div which might cause serious trouble
            nct_ibz_neg(&coeff_2, &coeff_2);
            //ensure this does nothing if a[j][i]=0
            nct_ibz_vec_4_linear_combination_mod(&(a[j]),&coeff_1,&(a[j]),&coeff_2,&(a[k]),&m);//do lin comb mod m
            nct_ibz_vec_4_copy_mod(&(a[k]),&c,&m); // mod m in copy
        }
        nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&(a[k][i]),&m);
        nct_ibz_vec_4_scalar_mul_mod(&(w[i]),&u,&(a[k]),&m);//mod m in scalar mult
        //needs modulo to be in [1,m] for wii
        nct_ibz_mod_not_zero(&(w[i][i]),&(w[i][i]),&m);
        for(int h=i+1;h<4; h++){
            nct_ibz_div_floor(&q,&r,&(w[h][i]),&(w[i][i]));
            nct_ibz_neg(&q,&q);
            nct_ibz_vec_4_linear_combination(&(w[h]),&one, &(w[h]),&q,&(w[i]));
        }
        nct_ibz_div_towards_zero(&m,&r,&m,&d);//modified line from nct_ibz_div which might cause serious trouble
        //assert(nct_ibz_is_zero(&r));
        //needs modulo to be in [1,m] for aki
        nct_ibz_mod_not_zero(&(a[k][i]),&(a[k][i]),&m);
    }
    for (j = 0; j < 4; j++) {
        for(i = 0; i < 4; i++){
            nct_ibz_copy(&((*hnf)[i][j]),&(w[j][i]));
        }
    }

    nct_ibz_finalize(&b);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&u);
    nct_ibz_finalize(&v);
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&q);
    nct_ibz_finalize(&coeff_1);
    nct_ibz_finalize(&coeff_2);
    nct_ibz_finalize(&m);
    nct_ibz_finalize(&one);
    nct_ibz_vec_4_finalize(&c);
    for (int h = 0; h < n; h++){
        if(h<4)
            nct_ibz_vec_4_finalize(&(w[h]));
        nct_ibz_vec_4_finalize(&(a[h]));
    }
}


void nct_quat_alg_equal_denom(nct_quat_alg_elem_t *res_a, nct_quat_alg_elem_t *res_b, const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b){
  nct_ibz_t gcd, r;
  nct_ibz_init(&gcd);
  nct_ibz_init(&r);
  nct_ibz_gcd(&gcd, &(a->denom), &(b->denom));
  //temporarily set res_a.denom to a.denom/gcd, and res_b.denom to b.denom/gcd
  nct_ibz_div_towards_zero(&(res_a->denom), &r, &(a->denom), &gcd);
  nct_ibz_div_towards_zero(&(res_b->denom), &r, &(b->denom), &gcd);
  for (int i = 0; i<4;i++){
    //multiply coordiates by reduced denominators from the other element
    nct_ibz_mul(&(res_a->coord[i]), &(a->coord[i]), &(res_b->denom));
    nct_ibz_mul(&(res_b->coord[i]), &(b->coord[i]), &(res_a->denom));
  }
  // multiply both reduced denominators
  nct_ibz_mul(&(res_a->denom), &(res_a->denom), &(res_b->denom));
  // multiply them by the gcd to get the new common denominator
  nct_ibz_mul(&(res_b->denom), &(res_a->denom), &gcd);
  nct_ibz_mul(&(res_a->denom), &(res_a->denom), &gcd);
  nct_ibz_finalize(&gcd);
  nct_ibz_finalize(&r);
}

//Public Functions

void nct_quat_alg_sub(nct_quat_alg_elem_t *res, const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b){
  nct_quat_alg_elem_t res_a, res_b;
  nct_quat_alg_elem_init(&res_a);
  nct_quat_alg_elem_init(&res_b);
  // put both on the same denominator
  nct_quat_alg_equal_denom(&res_a,&res_b,a,b);
  //then substract
  nct_ibz_copy(&res->denom, &res_a.denom);
  nct_ibz_vec_4_sub(&res->coord,&res_a.coord,&res_b.coord);
  nct_quat_alg_elem_finalize(&res_a);
  nct_quat_alg_elem_finalize(&res_b);
}


void nct_quat_alg_normalize(nct_quat_alg_elem_t *x){
  nct_ibz_t gcd,sign, r,zero;
  nct_ibz_init(&gcd);
  nct_ibz_init(&sign);
  nct_ibz_init(&r);
  nct_ibz_init(&zero);
  nct_ibz_set(&zero,0);
  nct_ibz_content(&gcd,&(x->coord));
  nct_ibz_gcd(&gcd,&gcd,&(x->denom));
  nct_ibz_div_towards_zero(&(x->denom),&r,&(x->denom),&gcd);
  nct_ibz_vec_4_scalar_div(&(x->coord),&gcd,&(x->coord));
  nct_ibz_set(&sign, 2*(0>nct_ibz_cmp(&zero,&(x->denom)))-1);
  nct_ibz_vec_4_scalar_mul(&(x->coord),&sign,&(x->coord));
  nct_ibz_mul(&(x->denom),&sign,&(x->denom));
  nct_ibz_finalize(&gcd);
  nct_ibz_finalize(&sign);
  nct_ibz_finalize(&r);
  nct_ibz_finalize(&zero);
}

int nct_quat_alg_elem_equal(const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b){
  nct_quat_alg_elem_t diff;
  nct_quat_alg_elem_init(&diff);
  nct_quat_alg_sub(&diff,a,b);
  int res = nct_quat_alg_elem_is_zero(&diff);
  nct_quat_alg_elem_finalize(&diff);
  return(res);
}

int nct_quat_alg_elem_is_zero(const nct_quat_alg_elem_t *x){
  int res = nct_ibz_vec_4_is_zero(&(x->coord));
  return(res);
}

// helper functions for lattices
void nct_quat_alg_elem_copy_ibz(nct_quat_alg_elem_t *elem, const nct_ibz_t *denom, const nct_ibz_t *coord0,const nct_ibz_t *coord1,const nct_ibz_t *coord2,const nct_ibz_t *coord3){
    nct_ibz_copy(&(elem->coord[0]), coord0);
    nct_ibz_copy(&(elem->coord[1]), coord1);
    nct_ibz_copy(&(elem->coord[2]), coord2);
    nct_ibz_copy(&(elem->coord[3]), coord3);

    nct_ibz_copy(&(elem->denom),denom);
}

void nct_quat_alg_elem_set(nct_quat_alg_elem_t *elem, int64_t denom, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3){
    nct_ibz_set(&(elem->coord[0]), coord0);
    nct_ibz_set(&(elem->coord[1]), coord1);
    nct_ibz_set(&(elem->coord[2]), coord2);
    nct_ibz_set(&(elem->coord[3]), coord3);

    nct_ibz_set(&(elem->denom),denom);
}
//return 1 if equal, 0 else.
int nct_quat_lattice_equal(const nct_quat_lattice_t *lat1, const nct_quat_lattice_t *lat2){
    int equal;
    nct_ibz_t abs_denom1, abs_denom2;
    nct_ibz_mat_4x4_t m1,m2;
    nct_ibz_init(&abs_denom1);
    nct_ibz_init(&abs_denom2);
    nct_ibz_mat_4x4_init(&m1);
    nct_ibz_mat_4x4_init(&m2);
    // test if both are in HNF as needed
    assert(nct_ibz_mat_4x4_is_hnf(&(lat1->basis)));
    assert(nct_ibz_mat_4x4_is_hnf(&(lat2->basis)));
    // get absolute values of denominators
    nct_ibz_abs(&abs_denom1,&(lat1->denom));
    nct_ibz_abs(&abs_denom2,&(lat2->denom));
    // cross-multiply by denomiators to get both basis on same denominators
    nct_ibz_mat_4x4_triangular_scalar_mul(&m1,&abs_denom2,&(lat1->basis));
    nct_ibz_mat_4x4_triangular_scalar_mul(&m2,&abs_denom1,&(lat2->basis));
    // both are still HNF, so simply test for equality
    equal = !nct_ibz_mat_4x4_equal(&m1,&m2);
    nct_ibz_finalize(&abs_denom1);
    nct_ibz_finalize(&abs_denom2);
    nct_ibz_mat_4x4_finalize(&m1);
    nct_ibz_mat_4x4_finalize(&m2);
    return(equal);
}

//sublattice test
int nct_quat_lattice_inclusion(const nct_quat_lattice_t *sublat, const nct_quat_lattice_t *overlat){
    int res;
    nct_quat_lattice_t sum;
    nct_quat_lattice_init(&sum);
    nct_quat_lattice_add(&sum,overlat,sublat);
    res = nct_quat_lattice_equal(&sum,overlat);
    nct_quat_lattice_finalize(&sum);
    return(res);
}

void nct_quat_lattice_reduce_denom(nct_quat_lattice_t *reduced, const nct_quat_lattice_t *lat){
    nct_ibz_t gcd;
    nct_ibz_init(&gcd);
    nct_ibz_mat_4x4_gcd(&gcd,&(lat->basis));
    nct_ibz_gcd(&gcd,&gcd,&(lat->denom));
    nct_ibz_mat_4x4_scalar_div(&(reduced->basis),&gcd,&(lat->basis));
    nct_ibz_div_towards_zero(&(reduced->denom),&gcd,&(lat->denom),&gcd);//modified from nct_ibz_div, might fail
    nct_ibz_abs(&(reduced->denom),&(reduced->denom));
    nct_ibz_finalize(&gcd);
}

void nct_quat_lattice_add(nct_quat_lattice_t *res, const nct_quat_lattice_t *lat1, const nct_quat_lattice_t *lat2){
    nct_ibz_vec_4_t generators[8];
    nct_ibz_mat_4x4_t tmp;
    nct_ibz_t det1, det2, detprod,dethnf;
    nct_ibz_init(&det1);
    nct_ibz_init(&det2);
    nct_ibz_init(&detprod);
    for(int i = 0; i <8; i++)
        nct_ibz_vec_4_init(&(generators[i]));
    nct_ibz_mat_4x4_init(&tmp);
    nct_ibz_mat_4x4_scalar_mul(&tmp,&(lat1->denom),&(lat2->basis));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_copy(&(generators[j][i]),&(tmp[i][j]));
        }
    }
    nct_ibz_mat_4x4_inv_with_det_as_denom(NULL,&det1,&tmp);
    nct_ibz_mat_4x4_scalar_mul(&tmp,&(lat2->denom),&(lat1->basis));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_copy(&(generators[4+j][i]),&(tmp[i][j]));
        }
    }
    nct_ibz_mat_4x4_inv_with_det_as_denom(NULL,&det2,&tmp);
    //assert(!nct_ibz_is_zero(&det1));
    //assert(!nct_ibz_is_zero(&det2));
    nct_ibz_gcd(&detprod,&det1,&det2);
    nct_ibz_mat_4xn_hnf_mod_core(&(res->basis),8,&generators,&detprod);
    nct_ibz_mul(&(res->denom),&(lat1->denom), &(lat2->denom));
    nct_quat_lattice_reduce_denom(res,res);
    nct_ibz_mat_4x4_finalize(&tmp);
    nct_ibz_finalize(&det1);
    nct_ibz_finalize(&det2);
    nct_ibz_finalize(&detprod);
    for(int i = 0; i <8; i++)
        nct_ibz_vec_4_finalize(&(generators[i]));
}


// lattice assumed of full rank and under HNF, none of both is tested so far
int nct_quat_lattice_contains_without_alg(nct_ibz_vec_4_t *coord, const nct_quat_lattice_t *lat, const nct_quat_alg_elem_t *x){
    int res = 1;
    nct_ibz_vec_4_t work_coord, work_x, column;
    nct_ibz_mat_4x4_t work_mat;
    nct_quat_alg_elem_t test;
    nct_ibz_t r, prod, one;
    nct_ibz_init(&r);
    nct_ibz_init(&prod);
    nct_ibz_init(&one);
    nct_quat_alg_elem_init(&test);
    nct_ibz_vec_4_init(&work_coord);
    nct_ibz_mat_4x4_init(&work_mat);
    nct_ibz_vec_4_init(&work_x);
    nct_ibz_vec_4_init(&column);
    // test if rank 4 lattice under HNF
//#ifndef NDEBUG
    //assert(nct_ibz_mat_4x4_is_hnf(&(lat->basis)));
    //for(int i = 0; i < 4; i++){
    //    assert(!nct_ibz_is_zero(&(lat->basis[i][i])));
    //}
//#endif
    nct_ibz_set(&one,1);
    for(int i = 0; i < 4;i++){
        // put on same denominator, 1st part
        nct_ibz_mul(&(work_x[i]), &(x->coord[i]),&(lat->denom));
        for(int j = i; j < 4;j++){
            // put on same denominator, 1st part
            nct_ibz_mul(&(work_mat[i][j]), &(lat->basis[i][j]),&(x->denom));
        }
    }
    for(int i = 0; i < 4;i++){
        nct_ibz_div_towards_zero(&(work_coord[3-i]), &r,&(work_x[3-i]), &(work_mat[3-i][3-i]));
        for (int j = 0; j < 4;j++){
            nct_ibz_mul(&prod,&(work_coord[3-i]),&(work_mat[j][3-i]));
            nct_ibz_sub(&(work_x[j]),&(work_x[j]),&prod);
        }
    }
    //final test: see if product gives expected result
    nct_ibz_mat_4x4_triangular_eval(&(test.coord),&(lat->basis),&work_coord);
    nct_ibz_copy(&(test.denom),&(lat->denom));
    res = nct_quat_alg_elem_equal(&test,x);
    //copy result
    if(coord != NULL){
        //output result iff res=1
        nct_ibz_set(&prod,res);
        nct_ibz_vec_4_scalar_mul(&work_coord,&prod,&work_coord);
        for(int i = 0; i < 4;i++){
            nct_ibz_copy(&((*coord)[i]),&(work_coord[i]));
        }
    }
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&one);
    nct_ibz_mat_4x4_finalize(&work_mat);
    nct_quat_alg_elem_finalize(&test);
    nct_ibz_vec_4_finalize(&work_coord);
    nct_ibz_vec_4_finalize(&work_x);
    nct_ibz_vec_4_finalize(&column);
    return(res);
}

void nct_quat_lattice_index(nct_ibz_t *index, const nct_quat_lattice_t *sublat, const nct_quat_lattice_t *overlat) {
    nct_ibz_t tmp,det;
    nct_ibz_init(&tmp);
    nct_ibz_init(&det);
    //assert(nct_ibz_mat_4x4_is_hnf(&overlat->basis));
    //assert(nct_ibz_mat_4x4_is_hnf(&sublat->basis));
    
    // det = det(sublat->basis)
    nct_ibz_mat_4x4_inv_with_det_as_denom(NULL,&det,&sublat->basis);
    // tmp = (overlat->denom)⁴
    nct_ibz_mul(&tmp, &overlat->denom, &overlat->denom);
    nct_ibz_mul(&tmp, &tmp, &tmp);
    // index = (overlat->denom)⁴ · det(sublat->basis)
    nct_ibz_mul(index, &det, &tmp);
    // tmp = (sublat->denom)⁴
    nct_ibz_mul(&tmp, &sublat->denom, &sublat->denom);
    nct_ibz_mul(&tmp, &tmp, &tmp);
    // det = det(overlat->basis)
    nct_ibz_mat_4x4_inv_with_det_as_denom(NULL,&det,&overlat->basis);
    // tmp = (sublat->denom)⁴ · det(overlat->basis)
    nct_ibz_mul(&tmp, &tmp, &det);
    // index = index / tmp
    nct_ibz_div_towards_zero(index, &tmp, index, &tmp);
    //assert(nct_ibz_is_zero(&tmp));
    // index = |index|
    nct_ibz_abs(index, index);
    
    nct_ibz_finalize(&tmp);
    nct_ibz_finalize(&det);
}

void nct_quat_lattice_hnf(nct_quat_lattice_t *lat){
    nct_ibz_t mod;
    nct_ibz_vec_4_t generators[4];
    nct_ibz_init(&mod);
    nct_ibz_mat_4x4_inv_with_det_as_denom(NULL,&mod,&(lat->basis));
    nct_ibz_abs(&mod,&mod);
    for(int i = 0; i <4; i++)
        nct_ibz_vec_4_init(&(generators[i]));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_copy(&(generators[j][i]),&(lat->basis[i][j]));
        }
    }
    nct_ibz_mat_4xn_hnf_mod_core(&(lat->basis),4, &generators,&mod);
    nct_quat_lattice_reduce_denom(lat,lat);
    nct_ibz_finalize(&mod);
    for(int i = 0; i <4; i++)
        nct_ibz_vec_4_finalize(&(generators[i]));
}


int nct_quat_alg_norm(nct_ibz_t *norm,const nct_quat_alg_elem_t *elem, const nct_quat_alg_t *alg){
    nct_ibz_t prod;
    nct_ibz_init(&prod);
    nct_ibz_mul(&prod,&(elem->coord[3]),&(elem->coord[3]));
    nct_ibz_mul(norm,&(elem->coord[2]),&(elem->coord[2]));
    nct_ibz_add(norm,norm,&prod);
    nct_ibz_mul(norm,norm,&(alg->p));
    nct_ibz_mul(&prod,&(elem->coord[1]),&(elem->coord[1]));
    nct_ibz_add(norm,norm,&prod);
    nct_ibz_mul(&prod,&(elem->coord[0]),&(elem->coord[0]));
    nct_ibz_add(norm,norm,&prod);
    nct_ibz_mul(&prod,&(elem->denom),&(elem->denom));
    nct_ibz_div_floor(norm,&prod,norm,&prod);
    int res = nct_ibz_is_zero(&prod);
    nct_ibz_finalize(&prod);
    return(res);
}

int nct_quat_lattice_ideal_norm(nct_ibz_t *norm, const nct_quat_lattice_t *lat, const nct_quat_alg_t *alg){
    nct_ibz_t new;
    int res = 1;
    nct_quat_alg_elem_t elem;
    nct_ibz_init(&new);
    nct_quat_alg_elem_init(&elem);
#ifndef NDEBUG
    for(int i = 1; i < 4; i++){
        for(int j = 0; j < i; j++){
            assert(nct_ibz_is_zero(&(lat->basis[i][j])));
        }
    }
#endif
    nct_ibz_set(&new,1);
    for(int i = 0; i < 4; i++){
        nct_ibz_mul(&new,&new,&(lat->basis[i][i]));
    }
    nct_ibz_set(norm,4);
    nct_ibz_mul(norm,&new,norm);
    for(int i = 0; i < 4; i++){
        nct_ibz_div_floor(norm,&new,norm,&(lat->denom));
        res = res & nct_ibz_is_zero(&new);
    }
#ifndef NDEBUG
    nct_ibz_copy(&new,norm);
#endif
    nct_ibz_sqrt_floor(norm,norm);
#ifndef NDEBUG
    nct_ibz_t test;
    nct_ibz_init(&test);
    nct_ibz_mul(&test,norm,norm);
    res = res & (0==nct_ibz_cmp(&test,&new));
    nct_ibz_copy(&(elem.denom),&(lat->denom));
    for(int i = 1; i < 4; i++){
        for(int j = 0; j < i; j++){
            nct_ibz_copy(&(elem.coord[j]),&(lat->basis[i][j]));
        }
        nct_quat_alg_norm(&new,&elem,alg);
        nct_ibz_div_floor(&new,&test,&new,norm);
        res = res & nct_ibz_is_zero(&test);
    }
    nct_ibz_finalize(&test);
    if(!res) printf("Failure in norm computation\n");
#endif
    nct_quat_alg_elem_finalize(&elem);
    nct_ibz_finalize(&new);
    return(res);
}

void nct_quat_lattice_rfactor(nct_ibq_t *rfactor, nct_ibq_t *b0norm, const nct_ibz_mat_4x4_t *mat, const nct_quat_alg_t *alg){
    nct_ibz_t den, num;
    nct_ibq_t tmp;
    nct_ibz_mat_4x4_t g;
    nct_ibz_mat_4x4_t basis;
    nct_ibq_mat_4x4_t rr;
    nct_ibq_mat_4x4_t mu;
    nct_ibz_mat_4x4_init(&g);
    nct_ibz_mat_4x4_init(&basis);
    nct_ibq_mat_4x4_init(&rr);
    nct_ibq_mat_4x4_init(&mu);
    nct_ibq_init(&tmp);
    nct_ibz_init(&den);
    nct_ibz_init(&num);
    // transpose to be able to use algebra elemnts (which are then rows) more easily as vectors
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_copy(&(basis[i][j]),&((*mat)[j][i]));
        }
    }
    // G is the gram matrix: Symmetric, so only compute half
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < i+1; j++){
            nct_ibz_vec_4_bilinear(&(g[i][j]),&(basis[i]),&(basis[j]),alg);
        }
    }
    //copied from GSO_compute in sage version
    nct_ibz_set(&den,1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < i; j++){
            nct_ibq_set(&(rr[i][j]),&(g[i][j]),&den);
            for(int k = 0; k < j+1; k++){
                nct_ibq_mul(&tmp,&(mu[j][k]),&(rr[i][k]));
                nct_ibq_sub(&(rr[i][j]),&(rr[i][j]),&tmp);
            }
            nct_ibq_div(&(mu[i][j]),&(rr[i][j]),&(rr[j][j]));
        }
        nct_ibq_set(&(rr[i][i]),&(g[i][i]),&den);
        for(int j = 0; j < i+1; j++){
                nct_ibq_mul(&tmp,&(mu[i][j]),&(rr[i][j]));
                nct_ibq_sub(&(rr[i][i]),&(rr[i][i]),&tmp);
        }
    }

    nct_ibq_set(rfactor,&nct_ibz_const_one,&nct_ibz_const_one);
    for(int i = 0; i < 4; i++){
        nct_ibq_mul(rfactor,rfactor,&(rr[i][i]));
    }

    nct_ibq_copy(b0norm, &(rr[0][0]));
    nct_ibq_mat_4x4_finalize(&mu);
    nct_ibq_mat_4x4_finalize(&rr);
    nct_ibz_mat_4x4_finalize(&basis);
    nct_ibz_mat_4x4_finalize(&g);
    nct_ibz_finalize(&den);
    nct_ibz_finalize(&num);
    nct_ibq_finalize(&tmp);
}

int nct_quat_test_length(const nct_quat_lattice_t *red, const nct_quat_alg_t *alg, int print_flag){
    int res = 0;
    nct_ibz_t num, denom;
    nct_ibq_t b0norm, rfactor, tmp;

    nct_ibz_init(&num);
    nct_ibz_init(&denom);

    nct_ibq_init(&b0norm);
    nct_ibq_init(&rfactor);
    nct_ibq_init(&tmp);

    nct_quat_lattice_rfactor(&rfactor, &b0norm, &(red->basis),alg);

    //finish computing the RHS: 1 * (4/3)^12
    nct_ibz_set(&num, 16777216);
    nct_ibz_set(&denom, 531441);
    nct_ibq_set(&tmp, &num, &denom);
    nct_ibq_mul(&rfactor, &rfactor, &tmp);

    //compute the LSH ||b0||^8. NOTE: nct_quat_lattice_rfactor returns ||b0||^2
    nct_ibq_mul(&b0norm, &b0norm, &b0norm);
    nct_ibq_mul(&b0norm, &b0norm, &b0norm);

    //test lemma 6
    res = res || !(nct_ibq_cmp(&b0norm, &rfactor) < 0);  // negative

    if(print_flag) {
        printf("Lemma6\n");
        nct_ibz_printf("rhs = %Qd\n", &rfactor);
        nct_ibz_printf("lhs = %Qd\n", &b0norm);
    }

    nct_ibz_finalize(&num);
    nct_ibz_finalize(&denom);

    nct_ibq_finalize(&tmp);
    nct_ibq_finalize(&rfactor);
    nct_ibq_finalize(&b0norm);
    return res;
}
