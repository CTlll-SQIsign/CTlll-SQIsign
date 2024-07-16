#include "ct_intbig.h"
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>


void ibz_finalize(ibz_t *x)
{}

void ibq_finalize(ibq_t *x)
{
    ibz_finalize(&((*x)[0]));
    ibz_finalize(&((*x)[1]));
}

void ibz_init(ibz_t *x){
    for (int i = 0; i < LIMBNUM; i++){
        (*x)[i] = (uint64_t)0;
    }

}

void ibq_init(ibq_t *x){
    ibz_init(&((*x)[0]));
    ibz_init(&((*x)[1]));
    ibz_set(&((*x)[1]),1);
}


void ibz_conditional_swap(ibz_t *a, ibz_t *b, int cond){
    mpn_cnd_swap(cond,*a,*b,LIMBNUM);
}

void ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b)
{
    mpn_add_n(*sum, *a, *b,LIMBNUM);
}

void ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b)
{
    mpn_sub_n(*diff, *a, *b,LIMBNUM);
}

void ibz_neg(ibz_t *n, const ibz_t *x){
    ibz_t zero;
    ibz_init(&zero);
    ibz_set(&zero,0);
    ibz_sub(n,&zero,x);
    ibz_finalize(&zero);
}

void ibz_abs(ibz_t *a, const ibz_t *x){
    ibz_t zero;
    ibz_init(&zero);
    ibz_set(&zero,0);
    int t_neg = (0>ibz_cmp(x,&zero));
    ibz_neg(&zero,x);
    ibz_copy(a,x);
    ibz_conditional_swap(a,&zero,t_neg);
    ibz_finalize(&zero);
}


void ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b)
{
    mp_limb_t res[2*LIMBNUM];
    mp_limb_t tp [MUL_SCRATCH];
    for(int i = 0; i < 2*LIMBNUM; i++) res[i] = 0;
    mpn_sec_mul (res, *a, LIMBNUM, *b, LIMBNUM, tp);
    for(int i = 0; i < LIMBNUM; i++) (*prod)[i] = res[i];
}


int ibz_limbnum(const ibz_t *x){
    uint64_t ff = (int64_t)-1;
    uint64_t zz = 0;
    int c = LIMBNUM;
    int t = 1;
    for (int i=0; i<LIMBNUM;i++){
        t = t & ((zz==(*x)[LIMBNUM-1-i])|(ff==(*x)[LIMBNUM-1-i]));
        c = c - t;
    }
    return(c);
}

//bad cases seem to happen when n negative. So try to work with absolute value of n, if n negative, use -q-1 and d-r: n=qd+r -n =-qd-r = -q(d+1)+d-r
void ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d){
    mp_limb_t tp [DIV_SCRATCH];
    ibz_t zero, dt, dtd, dq, dr;
    ibz_init(&zero);
    ibz_init(&dt);
    ibz_init(&dtd);
    ibz_init(&dq);
    ibz_init(&dr);
    ibz_set(&zero,0);
    int signn = 1-2*(ibz_cmp(n,&zero)<0);
    int signd = 1-2*(ibz_cmp(d,&zero)<0);
    ibz_abs(&dt,n);
    ibz_abs(&dtd,d);
    int dlen = ibz_limbnum(&dtd);
    int res = mpn_sec_div_qr (dq, dt,LIMBNUM, dtd, dlen, tp);
    //need to sanitize outputs: check most significant limb to q, set highest ones of r to 0 and q back to 0
    assert(res==0 | (res = (int64_t)-1));
    for (int i=0;i<LIMBNUM;i++){
        int t = (i<dlen);
        dt[i] = dt[i]*(t);
    }
    ibz_copy(r,&dt);
    for (int i=0;i<LIMBNUM;i++){
        int t = (i<LIMBNUM-dlen);
        (dq)[i] = (dq)[i]*(t);
    }
    ibz_set(&dt,1-ibz_is_zero(r));
    ibz_add(&dt,&dq,&dt);
    ibz_neg(&dt,&dt);
    ibz_conditional_swap(&dq,&dt,signn<0);
    ibz_abs(&dt,d);
    ibz_sub(&dt,&dt,r);
    ibz_conditional_swap(r,&dt,(signn<0)*(1-ibz_is_zero(r)));

    ibz_set(&dt,signd);
    ibz_mul(q,&dt,&dq);

    ibz_finalize(&zero);
    ibz_finalize(&dt);
    ibz_finalize(&dtd);
    ibz_finalize(&dq);
    ibz_finalize(&dr);
}

//only valid if remainder is zero
void ibz_internal_div_no_remainder(ibz_t *q, const ibz_t *n, const ibz_t *d){
    ibz_t zero, dt, dtd, dq;
    mp_limb_t tp [DIV_SCRATCH];
    ibz_init(&zero);
    ibz_init(&dt);
    ibz_init(&dtd);
    ibz_init(&dq);
    ibz_set(&zero,0);
    int signn = 1-2*(ibz_cmp(n,&zero)<0);
    int signd = 1-2*(ibz_cmp(d,&zero)<0);
    ibz_abs(&dt,n);
    ibz_abs(&dtd,d);
    int dlen = ibz_limbnum(&dtd);
    int res = mpn_sec_div_qr (dq, dt,LIMBNUM, dtd, dlen, tp);
    assert(res==0 | (res = (int64_t)-1));
    for (int i=0;i<LIMBNUM;i++){
        int t = (i<dlen);
        dt[i] = dt[i]*(t);
    }
    assert(ibz_is_zero(&dt));
    for (int i=0;i<LIMBNUM;i++){
        int t = (i<LIMBNUM-dlen);
        (dq)[i] = (dq)[i]*(t);
    }
    ibz_neg(&dt,&dq);
    ibz_conditional_swap(&dt,&dq,signd*signn<0);
    ibz_copy(q,&dq);
    ibz_finalize(&zero);
    ibz_finalize(&dt);
    ibz_finalize(&dtd);
    ibz_finalize(&dq);
}

int ibz_cmp(const ibz_t *a, const ibz_t *b)
{
    //try to compute diff, and compare result to remaining
    ibz_t diff;
    ibz_init(&diff);
    ibz_sub(&diff,a,b);
    int neg = (diff[LIMBNUM-1]>>(sizeof(mp_limb_t)-1))&1;
    int zero = ibz_is_zero(&diff);
    int ret = (1-zero)*(1-2*neg);
    ibz_finalize(&diff);
    return ret;
}

int ibz_is_zero(const ibz_t *x)
{
    int res = mpn_zero_p(*x,LIMBNUM);
    return(res);
}

void ibz_set(ibz_t *i, int64_t x)
{
    (*i)[0]=  (int64_t)x;
    int t = (x>=0);
    t = t-1;
    for(int j = 1; j < LIMBNUM; j++){
        (*i)[j]=  (uint64_t)t;
    }
}

//only for positive numbers
int ibz_set_from_str(ibz_t *i, const char *str, int base)
{
    //strsize must be at most such that the largest number in base 10 fits in the LIMBNUM-1 limbs
    int len = strlen(str);
    if(len>STRSIZE){
        return 0;
    }
    char *newstr = malloc(len);
    newstr[len-1] = 0;
    for (int j = 0; j < len; j++){
        newstr[j]=str[j]-48;
    }
    ibz_set(i,0);
    int res = mpn_set_str(*i, newstr,len,base);
    int validres = (res+1<LIMBNUM);
    free(newstr);
    return(validres);
}

void ibz_copy(ibz_t *target, const ibz_t *value)
{
    mpn_copyi(*target,*value,LIMBNUM);
}

// r = x mod m and 0 <=r<abs(m)
void ibz_mod(ibz_t *r, const ibz_t *x, const ibz_t *m){
    ibz_t a;
    ibz_init(&a);
    ibz_abs(&a,m);
    ibz_div_floor(&a,r,x,&a);
    ibz_finalize(&a);
}

void ibz_bitwise_and(ibz_t *res, const ibz_t *a, const ibz_t *b) {
    for (int i = 0; i < LIMBNUM; ++i) {
        (*res)[i] = (*a)[i] & (*b)[i];
    }
}

//leaks t
void ibz_internal_truncate(ibz_t *x, int t){
    if(t==0){
        ibz_set(x,0);
        return;
    }
    if(t<=0 | t>=64*LIMBNUM){
        return;
    }
    int sign = ((*x)[LIMBNUM-1]>>63) &1;
    uint64_t l = 0xffffffffffffffff;
    uint64_t m = 0xffffffffffffffff;
    l = l*sign;
    int limb_n = (int)(t/64);
    int small = t%64;
    for(int i = LIMBNUM-1; i>limb_n;i--){
        (*x)[i] = l;
    }
    l = l <<(small);
    m = m >> (64-small);
    (*x)[(limb_n)] = l + (m&(*x)[limb_n]);
    for(int i =limb_n-1; i>-1;i--){
        (*x)[i] = (uint64_t)0;
    }
}

void ibz_internal_one_shift_right(ibz_t *x){
    for(int i = 0; i <LIMBNUM-1;i++){
        (*x)[i] = ((*x)[i]>>1) + (((*x)[i+1]&1)<<63);
    }
    (*x)[LIMBNUM-1] = ((*x)[LIMBNUM-1]>>1) + ((*x)[LIMBNUM-1]&(((uint64_t)1)<<63));
}

void ibz_internal_divsteps(ibz_t *f, ibz_t *g){
    assert(GCD_CONST_1 <= GCD_CONST_2);
    assert(0 <= GCD_CONST_1);
    int t = GCD_CONST_2;
    int delta = 1;
    int test;
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_internal_truncate(f,t);
    ibz_internal_truncate(g,t);
    for(int i = 0; i<GCD_CONST_1;i++){
        ibz_internal_truncate(f,t);
        test = ((delta>0) & ((*g)[0]&1));
        delta = delta*(1-2*test);
        ibz_neg(&tmp,f);
        mpn_cnd_swap(test,(*g),tmp,LIMBNUM);
        mpn_cnd_swap(test,(*f),tmp,LIMBNUM);
        delta = 1+delta;
        ibz_add(&tmp,f,g);
        ibz_conditional_swap(&tmp,g,(*g[0])&1);
        ibz_internal_one_shift_right(g);
        t = t-1;
        ibz_internal_truncate(g,t);
    }
    ibz_finalize(&tmp);
}

void ibz_remove_factor_two(ibz_t *factor, ibz_t *a, ibz_t *b){
    int ctr = 0;
    int not_found = 1;
    int not_found2 = 1;
    uint64_t j = 1;
    uint64_t target = 0;
    uint64_t part_fact = 0;
    for(int i = 0; i < LIMBNUM; i++){
        target = ((uint64_t)not_found)*((*a)[i] | (*b)[i]) + target*((uint64_t)(1-not_found));
        not_found = not_found & (target==0);
        ctr = ctr + not_found;
    }
    for(int i = 0; i < 64; i++){
        part_fact = not_found2*j + (1-not_found2)*part_fact;
        not_found2 = not_found2 & (0==(target&j));
        j = j <<1;
    }
    for(int i = 0; i < LIMBNUM; i++){
        (*factor)[i] = part_fact*(ctr==i);
    }
    assert(!ibz_is_zero(factor));
    ibz_internal_div_no_remainder(a,a,factor);
    ibz_internal_div_no_remainder(b,b,factor);
}


void ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b)
{
    ibz_t at, bt, res;
    ibz_init(&at);
    ibz_init(&bt);
    ibz_init(&res);
    ibz_abs(&at,a);
    ibz_abs(&bt,b);
    ibz_set(&res,1);
    ibz_remove_factor_two(&res,&at,&bt);
    int test = !(at[0]&1);
    mpn_cnd_swap(test,at,bt,LIMBNUM);
    assert(at[0]&1);
    ibz_internal_divsteps(&at, &bt);
    ibz_abs(&at,&at);
    ibz_mul(gcd, &at,&res);
    ibz_finalize(&res);
    ibz_finalize(&at);
    ibz_finalize(&bt);
}

void ibz_swap(ibz_t *a, ibz_t *b){
    mpn_cnd_swap(1,*a,*b,LIMBNUM);
}

//helper function: reduce fraction
void ibq_red(ibq_t *red){
    ibz_t gcd, r;
    ibz_init(&gcd);
    ibz_init(&r);
    ibz_gcd(&gcd,&((*red)[0]),&((*red)[1]));
    ibz_div_floor(&((*red)[0]),&r,&((*red)[0]),&gcd);
    assert(ibz_is_zero(&r));
    ibz_div_floor(&((*red)[1]),&r,&((*red)[1]),&gcd);
    assert(ibz_is_zero(&r));
    ibz_finalize(&gcd);
    ibz_finalize(&r);
}


void ibq_add(ibq_t *sum, const ibq_t *a, const ibq_t *b)
{
    ibz_t y;
    ibq_t x;
    ibz_init(&y);
    ibq_init(&x);
    ibz_mul(&(x[0]),&((*a)[0]),&((*b)[1]));
    ibz_mul(&(x[1]),&((*a)[1]),&((*b)[1]));
    ibz_mul(&y,&((*a)[1]),&((*b)[0]));
    ibz_add(&(x[0]),&(x[0]),&y);
    ibq_copy(sum,&x);
    ibq_red(sum);
    ibz_finalize(&y);
    ibq_finalize(&x);
}


void ibq_sub(ibq_t *diff, const ibq_t *a, const ibq_t *b)
{
    ibq_t x;
    ibq_init(&x);
    ibq_neg(&x,b);
    ibq_add(diff,a,&x);
    ibq_finalize(&x);
}


void ibq_neg(ibq_t *neg, const ibq_t *x){
    ibz_t zero;
    ibz_init(&zero);
    ibz_sub(&((*neg)[0]),&zero,&((*x)[0]));
    ibz_copy(&((*neg)[1]),&((*x)[1]));
    ibz_finalize(&zero);
}

void ibq_abs(ibq_t *abs, const ibq_t *x){
    ibz_t sign,zero;
    ibz_init(&sign);
    ibz_init(&zero);
    ibz_set(&zero,0);
    ibz_mul(&sign,&((*x)[0]),&((*x)[1]));
    int s = ibz_cmp(&sign,&zero)>0;
    s = 2*s-1;
    ibz_set(&sign,s);
    ibz_mul(&((*abs)[0]),&((*x)[0]),&sign);
    ibz_copy(&((*abs)[1]),&((*x)[1]));
    ibz_finalize(&sign);
    ibz_finalize(&zero);
}

void ibq_mul(ibq_t *prod, const ibq_t *a, const ibq_t *b)
{
    ibz_mul(&((*prod)[0]),&((*a)[0]),&((*b)[0]));
    ibz_mul(&((*prod)[1]),&((*a)[1]),&((*b)[1]));
    ibq_red(prod);
}

void ibq_div(ibq_t *quot, const ibq_t *a, const ibq_t *b){
    ibq_t q;
    ibq_init(&q);
    ibz_mul(&(q[0]),&((*a)[0]),&((*b)[1]));
    ibz_mul(&(q[1]),&((*a)[1]),&((*b)[0]));
    ibq_red(&q);
    ibq_copy(quot,&q);
    ibq_finalize(&q);
}

int ibq_cmp(const ibq_t *a, const ibq_t *b)
{
    //careful about sign, which reducing does not manage
    ibq_t diff;
    ibz_t d,zero;
    ibq_init(&diff);
    ibz_init(&zero);
    ibz_init(&d);
    ibz_set(&zero,0);
    ibq_sub(&diff,a,b);
    ibq_denom(&d,&diff);
    int signd = ibz_cmp(&d,&zero);
    ibq_num(&d,&diff);
    int signn = ibz_cmp(&d,&zero);
    int res = signd*signn;
    ibq_finalize(&diff);
    ibz_finalize(&zero);
    ibz_finalize(&d);
    return(res);
}

int ibq_is_zero(const ibq_t *x)
{
    return ibz_is_zero(&((*x)[0]));
}

int ibq_set(ibq_t *q, const ibz_t *a, const ibz_t *b)
{
    ibz_t one, bi;
    ibz_init(&one);
    ibz_init(&bi);
    one[0]=1;
    ibz_copy(&bi,b);
    int t_swap = ibz_is_zero(b);
    ibz_conditional_swap(&bi,&one,t_swap);
    ibz_copy(&((*q)[0]),a);
    ibz_copy(&((*q)[1]),&bi);
    ibq_red(q);
    return(!ibz_is_zero(b));
}

void ibq_copy(ibq_t *target, const ibq_t *value)
{
    ibz_copy(&((*target)[0]),&((*value)[0]));
    ibz_copy(&((*target)[1]),&((*value)[1]));
    ibq_red(target);
}

void ibq_swap(ibq_t *a, ibq_t *b){
    ibz_swap(&((*a)[0]),&((*b)[0]));
    ibz_swap(&((*a)[1]),&((*b)[1]));
}

void ibq_denom(ibz_t *d, const ibq_t *x)
{
    ibz_copy(d,&((*x)[1]));
}

void ibq_num(ibz_t *n, const ibq_t *x)
{
    ibz_copy(n,&((*x)[0]));
}

int ibq_to_ibz(ibz_t *z, const ibq_t *q) {
    ibz_t n, d, quo, rem;
    int res;

    ibz_init(&n);
    ibz_init(&d);
    ibz_init(&quo);
    ibz_init(&rem);

    ibq_num(&n, q);
    ibq_denom(&d, q);

    ibz_div_floor(&quo, &rem, &n, &d);
    ibz_copy(z, &quo);
    res = ibz_is_zero(&rem);

    ibz_finalize(&n);
    ibz_finalize(&d);
    ibz_finalize(&quo);
    ibz_finalize(&rem);
    return res;
}


void ibz_mat_4x4_transpose(ibz_mat_4x4_t *transposed, const ibz_mat_4x4_t *mat){
    ibz_mat_4x4_t work;
    ibz_mat_4x4_init(&work);
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            ibz_copy(&(work[i][j]),&((*mat)[j][i]));
        }
    }
    ibz_mat_4x4_copy(transposed,&work);
    ibz_mat_4x4_finalize(&work);
}

void ibz_mat_4x4_copy(ibz_mat_4x4_t *copy, const ibz_mat_4x4_t *copied){
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            ibz_copy(&((*copy)[i][j]),&((*copied)[i][j]));
        }
    }
}

void ibz_mat_4x4_init(ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_init(&(*mat)[i][j]);
        }
    }
}
void ibz_mat_2x2_init(ibz_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            ibz_init(&(*x)[i][j]);
        }
    }
}
void ibq_mat_2x2_init(ibq_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            ibq_init(&(*x)[i][j]);
        }
    }
}
void ibq_mat_4x4_init(ibq_mat_4x4_t *x){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibq_init(&(*x)[i][j]);
        }
    }
}
void ibz_vec_4_init(ibz_vec_4_t *x){
    for(int j = 0; j < 4; j++){
        ibz_init(&(*x)[j]);
    }
}

void ibz_mat_4x4_finalize(ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_finalize(&(*mat)[i][j]);
        }
    }
}
void ibz_mat_2x2_finalize(ibz_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            ibz_finalize(&(*x)[i][j]);
        }
    }
}
void ibq_mat_2x2_finalize(ibq_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            ibq_finalize(&(*x)[i][j]);
        }
    }
}
void ibq_mat_4x4_finalize(ibq_mat_4x4_t *x){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibq_finalize(&(*x)[i][j]);
        }
    }
}
void ibz_vec_4_finalize(ibz_vec_4_t *x){
    for(int j = 0; j < 4; j++){
        ibz_finalize(&(*x)[j]);
    }
}


//helper

//result char must have size PRINTSIZE
void ibz_to_str(char *result, const ibz_t *z){
    char format[3] = "%s";
    format[2]=0;
    char out[PRINTSIZE] = "abc\n";
    ibz_t cp;
    ibz_init(&cp);
    ibz_copy(&cp,z);
    int len = mpn_get_str(out,10,cp,LIMBNUM);
    for (int i=0; i < len; i++){
        out[i] = out[i]+48;
    }
    out[len]=0;
    ibz_finalize(&cp);
    sprintf(result,format, out);
}

//needs %s for integer in format
void ibz_printf_1(char *format, const ibz_t *z){
    char out[PRINTSIZE] = "abc\n";
    ibz_t cp;
    ibz_init(&cp);
    ibz_copy(&cp,z);
    int len = mpn_get_str(out,10,cp,LIMBNUM);
    //printf("len %d\n",len);
    for (int i=0; i < len; i++){
        out[i] = out[i]+48;
    }
    out[len]=0;
    ibz_finalize(&cp);
    printf(format, out);
}

//needs 2 %s for num, denom in format
void ibq_printf_1(char *format, const ibq_t *q){
    char out[PRINTSIZE] = "abc\n";
    char out_d[PRINTSIZE] = "abc\n";
    ibz_t cp;
    ibz_init(&cp);
    ibz_copy(&cp,&((*q)[0]));
    int len = mpn_get_str(out,10,cp,LIMBNUM);
    ibz_copy(&cp,&((*q)[1]));
    int len_d = mpn_get_str(out_d,10,cp,LIMBNUM);
    //printf("len %d\n",len);
    for (int i=0; i < len; i++){
        out[i] = out[i]+48;
    }
    for (int i=0; i < len_d; i++){
        out_d[i] = out_d[i]+48;
    }
    out[len]=0;
    out_d[len_d]=0;
    ibz_finalize(&cp);
    printf(format, out, out_d);
}

void ibz_print_limbs(const ibz_t *z){
    for(int i = 0; i<LIMBNUM; i++){
        printf("%"PRIx64", ",(*z)[i]);
    }
    printf("\n");
}