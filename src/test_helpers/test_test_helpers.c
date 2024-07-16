
#include "test_test_helpers.h"

int nct_quat_test_xgcd_non_ct(){
    int res = 0;
    nct_ibz_t u,v,n,d,gcd;
    nct_ibz_init(&u);
    nct_ibz_init(&v);
    nct_ibz_init(&n);
    nct_ibz_init(&d);
    nct_ibz_init(&gcd);
    nct_ibz_set(&n,100);
    nct_ibz_set(&d,80);
    nct_ibz_xgcd_non_ct(&gcd,&u,&v,&n,&d);
    res = res || nct_quat_test_xgcd_verify(&gcd,&u,&v,&d,&n);
    nct_ibz_set(&n,-120);
    nct_ibz_set(&d,80);
    nct_ibz_xgcd_non_ct(&gcd,&u,&v,&n,&d);
    res = res || nct_quat_test_xgcd_verify(&gcd,&u,&v,&d,&n);
    if(res){
        printf("Test xgcd failed\n");
    }
    nct_ibz_finalize(&u);
    nct_ibz_finalize(&v);
    nct_ibz_finalize(&n);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&gcd);
    return(res);
}


//int nct_ibz_vec_4_scalar_div(nct_ibz_vec_4_t *quot, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec);
int nct_quat_test_dim4_nct_ibz_vec_4_scalar_div(){
    int res = 0;
    int s;
    nct_ibz_t scalar;
    nct_ibz_vec_4_t quot,vec;
    nct_ibz_vec_4_init(&vec);
    nct_ibz_vec_4_init(&quot);
    nct_ibz_init(&scalar);

    s = 5;
    nct_ibz_set(&scalar,s);
    for(int i = 0; i < 4; i++){
        nct_ibz_set(&(vec[i]),(i)*s);
    }
    res = res || !nct_ibz_vec_4_scalar_div(&quot,&scalar,&vec);
    for(int i = 0; i < 4; i++){
        res = res || (0!=nct_quat_test_helper_nct_ibz_equal_i(&(quot[i]),i));
    }

    res = res || !nct_ibz_vec_4_scalar_div(&vec,&scalar,&vec);
    for(int i = 0; i < 4; i++){
        res = res || (0!=nct_quat_test_helper_nct_ibz_equal_i(&(vec[i]),i));
    }

    if (res != 0){
        printf("Quaternion unit test dim4_nct_ibz_vec_4_scalar_div failed\n");
    }
    nct_ibz_vec_4_finalize(&vec);
    nct_ibz_vec_4_finalize(&quot);
    nct_ibz_finalize(&scalar);
    return(res);
}


//int nct_ibz_mat_4x4_scalar_div(nct_ibz_mat_4x4_t *quot, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat);
int nct_quat_test_dim4_nct_ibz_mat_4x4_scalar_div(){
    int res = 0;
    int s;
    nct_ibz_t scalar;
    nct_ibz_mat_4x4_t quot,mat,cmp;
    nct_ibz_mat_4x4_init(&mat);
    nct_ibz_mat_4x4_init(&cmp);
    nct_ibz_mat_4x4_init(&quot);
    nct_ibz_init(&scalar);

    s = 5;
    nct_ibz_set(&scalar,s);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_set(&(mat[i][j]),(i+j)*s);
            nct_ibz_set(&(cmp[i][j]),(i+j));
        }
    }
    nct_ibz_mat_4x4_scalar_div(&quot,&scalar,&mat);
    res = res || (nct_ibz_mat_4x4_equal(&quot,&cmp));

    nct_ibz_mat_4x4_scalar_div(&mat,&scalar,&mat);
    res = res || (nct_ibz_mat_4x4_equal(&mat,&cmp));

    if (res != 0){
        printf("Quaternion unit test dim4_nct_ibz_mat_4x4_scalar_div failed\n");
    }
    nct_ibz_mat_4x4_finalize(&mat);
    nct_ibz_mat_4x4_finalize(&cmp);
    nct_ibz_mat_4x4_finalize(&quot);
    nct_ibz_finalize(&scalar);
    return(res);
}

//void nct_ibz_mod_not_zero(nct_ibz_t *res, const nct_ibz_t *x, const nct_ibz_t *mod);
int nct_quat_test_integer_mod_not_zero(){
    int res = 0;
    nct_ibz_t m, a,r, d, one;
    int s,x;
    nct_ibz_init(&m);
    nct_ibz_init(&a);
    nct_ibz_init(&r);
    nct_ibz_init(&d);
    nct_ibz_init(&one);
    nct_ibz_set(&one,1);
    s = 5;
    nct_ibz_set(&m,s);
    for(x = -20; x < 20; x++){
        nct_ibz_set(&a,x);
        nct_ibz_mod_not_zero(&r,&a,&m);
        nct_ibz_sub(&d,&r,&a);
        res = res || !nct_ibz_divides(&d,&m);
        res = res || !(nct_ibz_cmp(&r,&m)<=0);
        res = res || !(nct_ibz_cmp(&one,&r)<=0);
        nct_ibz_mod_not_zero(&a,&a,&m);
        res = res || !(nct_ibz_cmp(&a,&r)==0);
        res = res || !(nct_ibz_cmp(&r,&m)<=0);
        res = res || !(nct_ibz_cmp(&one,&r)<=0);
    }


    if (res != 0){
        printf("Quaternion unit test integer_nct_ibz_mod_not_zero failed\n");
    }
    nct_ibz_finalize(&m);
    nct_ibz_finalize(&a);
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&one);
    return(res);
}

//void nct_ibz_centered_mod(nct_ibz_t *remainder, const nct_ibz_t *a, const nct_ibz_t *mod);
int nct_quat_test_integer_nct_ibz_centered_mod(){
    int res = 0;
    nct_ibz_t m, a,r, d, h, two;
    int s,x;
    nct_ibz_init(&m);
    nct_ibz_init(&a);
    nct_ibz_init(&r);
    nct_ibz_init(&d);
    nct_ibz_init(&h);
    nct_ibz_init(&two);
    nct_ibz_set(&two,2);
    s = 5;
    nct_ibz_set(&m,s);
    for(x = -20; x < 20; x++){
        nct_ibz_set(&a,x);
        nct_ibz_centered_mod(&r,&a,&m);
        nct_ibz_sub(&d,&r,&a);
        res = res || !nct_ibz_divides(&d,&m);
        nct_ibz_mul(&h,&r,&two);
        res = res || !(nct_ibz_cmp(&h,&m)<=0);
        nct_ibz_neg(&h,&h);
        res = res || !(nct_ibz_cmp(&h,&m)<0);
        nct_ibz_centered_mod(&a,&a,&m);
        res = res || !(nct_ibz_cmp(&a,&r)==0);
    }

    
    if (res != 0){
        printf("Quaternion unit test integer_nct_ibz_centered_mod failed\n");
    }
    nct_ibz_finalize(&m);
    nct_ibz_finalize(&a);
    nct_ibz_finalize(&r);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&h);
    nct_ibz_finalize(&two);
    return(res);
}

//void nct_ibz_conditional_assign(nct_ibz_t *res, const nct_ibz_t *x, const nct_ibz_t *y, int c);
int nct_quat_test_integer_nct_ibz_conditional_assign(){
    int res = 0;
    nct_ibz_t x,y,r;
    nct_ibz_init(&x);
    nct_ibz_init(&y);
    nct_ibz_init(&r);
    nct_ibz_set(&x,5);
    nct_ibz_set(&y,-9);
    nct_ibz_conditional_assign(&r,&x,&y,1);
    res = res || !(nct_ibz_cmp(&x,&r)==0);
    nct_ibz_conditional_assign(&r,&x,&y,0);
    res = res || !(nct_ibz_cmp(&y,&r)==0);
    nct_ibz_conditional_assign(&x,&x,&y,0);
    res = res || !(nct_ibz_cmp(&y,&x)==0);
    nct_ibz_set(&x,-5);
    nct_ibz_set(&y,-0);
    nct_ibz_conditional_assign(&y,&x,&y,1);
    res = res || !(nct_ibz_cmp(&y,&x)==0);


    if (res != 0){
        printf("Quaternion unit test integer_nct_ibz_conditional_assign failed\n");
    }
    nct_ibz_finalize(&x);
    nct_ibz_finalize(&y);
    nct_ibz_finalize(&r);
    return(res);
}

//void nct_ibz_xgcd_with_u_not_0(nct_ibz_t *d, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *x, const nct_ibz_t *y);
int nct_quat_test_integer_nct_ibz_xgcd_with_u_not_0(){
    int res = 0;
    nct_ibz_t x,y,u,v,d,s,p, zero, one;
    nct_ibz_init(&x);
    nct_ibz_init(&y);
    nct_ibz_init(&u);
    nct_ibz_init(&v);
    nct_ibz_init(&d);
    nct_ibz_init(&s);
    nct_ibz_init(&p);
    nct_ibz_init(&zero);
    nct_ibz_init(&one);
    nct_ibz_set(&zero,0);
    nct_ibz_set(&one,1);

    nct_ibz_set(&x, 75);
    nct_ibz_set(&y, 50);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,25));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&x);
    nct_ibz_neg(&p,&p);
    nct_ibz_rounded_div(&p,&p,&d);
    nct_ibz_abs(&s,&y);
    nct_ibz_rounded_div(&s,&s,&d);
    nct_ibz_set(&x, 1-2*(nct_ibz_cmp(&x,&zero)<0));
    nct_ibz_mul(&u,&u,&x);
    nct_ibz_set(&y, 1-2*(nct_ibz_cmp(&y,&zero)<0));
    nct_ibz_mul(&v,&v,&y);
    res = res || !(nct_ibz_cmp(&v,&p)>0);
    res = res || !(nct_ibz_cmp(&v,&zero)<=0);
    res = res || !(nct_ibz_cmp(&u,&s)<=0);
    res = res || !(nct_ibz_cmp(&u,&one)>=0);

    nct_ibz_set(&x, -75);
    nct_ibz_set(&y, 50);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,25));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&x);
    nct_ibz_neg(&p,&p);
    nct_ibz_rounded_div(&p,&p,&d);
    nct_ibz_abs(&s,&y);
    nct_ibz_rounded_div(&s,&s,&d);
    nct_ibz_set(&x, 1-2*(nct_ibz_cmp(&x,&zero)<0));
    nct_ibz_mul(&u,&u,&x);
    nct_ibz_set(&y, 1-2*(nct_ibz_cmp(&y,&zero)<0));
    nct_ibz_mul(&v,&v,&y);
    res = res || !(nct_ibz_cmp(&v,&p)>0);
    res = res || !(nct_ibz_cmp(&v,&zero)<=0);
    res = res || !(nct_ibz_cmp(&u,&s)<=0);
    res = res || !(nct_ibz_cmp(&u,&one)>=0);

    nct_ibz_set(&x, -75);
    nct_ibz_set(&y, -50);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,25));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&x);
    nct_ibz_neg(&p,&p);
    nct_ibz_rounded_div(&p,&p,&d);
    nct_ibz_abs(&s,&y);
    nct_ibz_rounded_div(&s,&s,&d);
    nct_ibz_set(&x, 1-2*(nct_ibz_cmp(&x,&zero)<0));
    nct_ibz_mul(&u,&u,&x);
    nct_ibz_set(&y, 1-2*(nct_ibz_cmp(&y,&zero)<0));
    nct_ibz_mul(&v,&v,&y);
    res = res || !(nct_ibz_cmp(&v,&p)>0);
    res = res || !(nct_ibz_cmp(&v,&zero)<=0);
    res = res || !(nct_ibz_cmp(&u,&s)<=0);
    res = res || !(nct_ibz_cmp(&u,&one)>=0);


    nct_ibz_set(&x, 75);
    nct_ibz_set(&y, -50);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,25));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&x);
    nct_ibz_neg(&p,&p);
    nct_ibz_rounded_div(&p,&p,&d);
    nct_ibz_abs(&s,&y);
    nct_ibz_rounded_div(&s,&s,&d);
    nct_ibz_set(&x, 1-2*(nct_ibz_cmp(&x,&zero)<0));
    nct_ibz_mul(&u,&u,&x);
    nct_ibz_set(&y, 1-2*(nct_ibz_cmp(&y,&zero)<0));
    nct_ibz_mul(&v,&v,&y);
    res = res || !(nct_ibz_cmp(&v,&p)>0);
    res = res || !(nct_ibz_cmp(&v,&zero)<=0);
    res = res || !(nct_ibz_cmp(&u,&s)<=0);
    res = res || !(nct_ibz_cmp(&u,&one)>=0);

    nct_ibz_set(&x, 50);
    nct_ibz_set(&y, 50);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,50));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&x);
    nct_ibz_neg(&p,&p);
    nct_ibz_rounded_div(&p,&p,&d);
    nct_ibz_abs(&s,&y);
    nct_ibz_rounded_div(&s,&s,&d);
    nct_ibz_set(&x, 1-2*(nct_ibz_cmp(&x,&zero)<0));
    nct_ibz_mul(&u,&u,&x);
    nct_ibz_set(&y, 1-2*(nct_ibz_cmp(&y,&zero)<0));
    nct_ibz_mul(&v,&v,&y);
    res = res || !(nct_ibz_cmp(&v,&p)>0);
    res = res || !(nct_ibz_cmp(&v,&zero)<=0);
    res = res || !(nct_ibz_cmp(&u,&s)<=0);
    res = res || !(nct_ibz_cmp(&u,&one)>=0);

    nct_ibz_set(&x, 0);
    nct_ibz_set(&y, -50);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,50));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&v);
    nct_ibz_abs(&s,&u);
    res = res || !(nct_ibz_cmp(&s,&one)<=1);
    res = res || !(nct_ibz_cmp(&p,&one)<=1);
    res = res || !(nct_ibz_cmp(&u,&zero)!=0);

    nct_ibz_set(&x, -50);
    nct_ibz_set(&y, 0);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,50));
    nct_ibz_mul(&s,&u,&x);
    nct_ibz_mul(&p,&v,&y);
    nct_ibz_add(&s,&p,&s);
    res = res || !(nct_ibz_cmp(&s,&d)==0);
    res = res || !(nct_ibz_cmp(&d,&zero)>0);
    nct_ibz_abs(&p,&v);
    nct_ibz_abs(&s,&u);
    res = res || !(nct_ibz_cmp(&s,&one)<=1);
    res = res || !(nct_ibz_cmp(&p,&one)<=1);
    res = res || !(nct_ibz_cmp(&u,&zero)!=0);
    
    nct_ibz_set(&x, 0);
    nct_ibz_set(&y, 0);
    nct_ibz_xgcd_with_u_not_0(&d,&u,&v,&x,&y);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&d,1));
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&u,1));
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&v,0));

    if (res != 0){
        printf("Quaternion unit test integer_nct_ibz_xgcd_with_u_not_0 failed\n");
    }
    nct_ibz_finalize(&x);
    nct_ibz_finalize(&y);
    nct_ibz_finalize(&u);
    nct_ibz_finalize(&v);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&s);
    nct_ibz_finalize(&p);
    nct_ibz_finalize(&zero);
    nct_ibz_finalize(&one);
    return(res);
}


//void nct_ibz_rounded_div(nct_ibz_t *q, const nct_ibz_t *a, const nct_ibz_t *b);
int nct_quat_test_integer_nct_ibz_rounded_div(){
    int res = 0;
    nct_ibz_t q, a, b;
    nct_ibz_init(&a);
    nct_ibz_init(&b);
    nct_ibz_init(&q);

    // basic tests
    nct_ibz_set(&a,15);
    nct_ibz_set(&b,3);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,5));
    nct_ibz_set(&a,16);
    nct_ibz_set(&b,3);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,5));
    nct_ibz_set(&a,17);
    nct_ibz_set(&b,3);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,6));
    nct_ibz_set(&a,37);
    nct_ibz_set(&b,5);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,7));
    // test sign combination
    nct_ibz_set(&a,149);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,12));
    nct_ibz_set(&a,149);
    nct_ibz_set(&b,-12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-12));
    nct_ibz_set(&a,-149);
    nct_ibz_set(&b,-12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,12));
    nct_ibz_set(&a,-149);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-12));
    nct_ibz_set(&a,151);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,13));
    nct_ibz_set(&a,-151);
    nct_ibz_set(&b,-12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,13));
    nct_ibz_set(&a,151);
    nct_ibz_set(&b,-12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-13));
    nct_ibz_set(&a,-151);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-13));
    //divisibles with sign
    nct_ibz_set(&a,144);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,12));
    nct_ibz_set(&a,-144);
    nct_ibz_set(&b,-12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,12));
    nct_ibz_set(&a,144);
    nct_ibz_set(&b,-12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-12));
    nct_ibz_set(&a,-144);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-12));
    // tests close to 0
    nct_ibz_set(&a,-12);
    nct_ibz_set(&b,-25);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,0));
    nct_ibz_set(&a,12);
    nct_ibz_set(&b,25);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,0));
    nct_ibz_set(&a,-12);
    nct_ibz_set(&b,25);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,0));
    nct_ibz_set(&a,12);
    nct_ibz_set(&b,-25);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,0));
    nct_ibz_set(&a,-12);
    nct_ibz_set(&b,-23);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,1));
    nct_ibz_set(&a,12);
    nct_ibz_set(&b,23);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,1));
    nct_ibz_set(&a,-12);
    nct_ibz_set(&b,23);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-1));
    nct_ibz_set(&a,12);
    nct_ibz_set(&b,-23);
    nct_ibz_rounded_div(&q,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&q,-1));
    // test output equal input
    nct_ibz_set(&a,-151);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&a,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&a,-13));
    nct_ibz_set(&a,-151);
    nct_ibz_set(&b,12);
    nct_ibz_rounded_div(&b,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&b,-13));
    // test for cmp not returning 1 or -1 or 0
    nct_ibz_set(&a,4292606433816540);
    nct_ibz_set(&b,864673106105940);
    nct_ibz_rounded_div(&b,&a,&b);
    res = res || !(0==nct_quat_test_helper_nct_ibz_equal_i(&b,5));

    if (res != 0){
        printf("Quaternion unit test integer_nct_ibz_rounded_div failed\n");
    }
    nct_ibz_finalize(&a);
    nct_ibz_finalize(&b);
    nct_ibz_finalize(&q);
    return(res);
}

//void nct_ibz_vec_4_linear_combination(nct_ibz_vec_4_t *lc, const nct_ibz_t *coeff_a, const nct_ibz_vec_4_t  *vec_a, const nct_ibz_t *coeff_b, const nct_ibz_vec_4_t *vec_b){
int nct_quat_test_dim4_nct_ibz_vec_4_linear_combination(){
    int res = 0;
    nct_ibz_vec_4_t a, b, lc, cmp;
    nct_ibz_t ca, cb;
    nct_ibz_init(&ca);
    nct_ibz_init(&cb);
    nct_ibz_vec_4_init(&a);
    nct_ibz_vec_4_init(&b);
    nct_ibz_vec_4_init(&lc);
    nct_ibz_vec_4_init(&cmp);
    nct_ibz_vec_4_set(&a,1,2,3,4);
    nct_ibz_vec_4_set(&b,-2,1,3,-3);
    nct_ibz_set(&ca,2);
    nct_ibz_set(&cb,-1);
    nct_ibz_vec_4_set(&cmp,4,3,3,11);
    nct_ibz_vec_4_linear_combination(&lc,&ca,&a,&cb,&b);
    for(int i = 0; i < 4; i++){
        res = res || nct_ibz_cmp(&(lc[i]),&(cmp[i]));
    }
    nct_ibz_vec_4_set(&cmp,1,2,3,4);
    nct_ibz_vec_4_linear_combination(&a,&ca,&a,&cb,&a);
    for(int i = 0; i < 4; i++){
        res = res || nct_ibz_cmp(&(a[i]),&(cmp[i]));
    }
    if (res != 0){
        printf("Quaternion unit test dim4_nct_ibz_vec_4_linear_combination failed\n");
    }
    nct_ibz_finalize(&ca);
    nct_ibz_finalize(&cb);
    nct_ibz_vec_4_finalize(&a);
    nct_ibz_vec_4_finalize(&b);
    nct_ibz_vec_4_finalize(&lc);
    nct_ibz_vec_4_finalize(&cmp);
    return(res);
}


//int nct_ibz_mat_4x4_is_hnf(const nct_ibz_mat_4x4_t *mat);
int nct_quat_test_dim4_is_hnf(){
    int res = 0;
    nct_ibz_mat_4x4_t mat;
    nct_ibz_mat_4x4_init(&mat);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            nct_ibz_set(&(mat[i][j]),0);
        }
    }
    res = res || (!nct_ibz_mat_4x4_is_hnf(&mat));
    nct_ibz_set(&(mat[0][0]),7);
    nct_ibz_set(&(mat[0][1]),6);
    nct_ibz_set(&(mat[0][2]),5);
    nct_ibz_set(&(mat[0][3]),4);
    nct_ibz_set(&(mat[1][1]),6);
    nct_ibz_set(&(mat[1][2]),5);
    nct_ibz_set(&(mat[1][3]),4);
    nct_ibz_set(&(mat[2][2]),5);
    nct_ibz_set(&(mat[2][3]),4);
    nct_ibz_set(&(mat[3][3]),4);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&mat));

    nct_ibz_set(&(mat[0][0]),7);
    nct_ibz_set(&(mat[0][1]),0);
    nct_ibz_set(&(mat[0][2]),5);
    nct_ibz_set(&(mat[0][3]),4);
    nct_ibz_set(&(mat[1][1]),0);
    nct_ibz_set(&(mat[1][2]),0);
    nct_ibz_set(&(mat[1][3]),0);
    nct_ibz_set(&(mat[2][2]),5);
    nct_ibz_set(&(mat[2][3]),4);
    nct_ibz_set(&(mat[3][3]),4);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&mat));

    // negative tests
    nct_ibz_set(&(mat[0][0]),7);
    nct_ibz_set(&(mat[0][1]),0);
    nct_ibz_set(&(mat[0][2]),5);
    nct_ibz_set(&(mat[0][3]),4);
    nct_ibz_set(&(mat[1][1]),1);
    nct_ibz_set(&(mat[1][2]),5);
    nct_ibz_set(&(mat[1][3]),9);
    nct_ibz_set(&(mat[2][2]),5);
    nct_ibz_set(&(mat[2][3]),4);
    nct_ibz_set(&(mat[3][3]),4);
    res = res || (nct_ibz_mat_4x4_is_hnf(&mat));

    nct_ibz_set(&(mat[0][0]),7);
    nct_ibz_set(&(mat[0][1]),0);
    nct_ibz_set(&(mat[0][2]),5);
    nct_ibz_set(&(mat[0][3]),4);
    nct_ibz_set(&(mat[1][1]),1);
    nct_ibz_set(&(mat[1][2]),-5);
    nct_ibz_set(&(mat[1][3]),1);
    nct_ibz_set(&(mat[2][2]),5);
    nct_ibz_set(&(mat[2][3]),4);
    nct_ibz_set(&(mat[3][3]),4);
    res = res || (nct_ibz_mat_4x4_is_hnf(&mat));


    nct_ibz_set(&(mat[0][0]),7);
    nct_ibz_set(&(mat[0][1]),0);
    nct_ibz_set(&(mat[0][2]),5);
    nct_ibz_set(&(mat[0][3]),4);
    nct_ibz_set(&(mat[1][0]),2);
    nct_ibz_set(&(mat[1][1]),3);
    nct_ibz_set(&(mat[1][2]),1);
    nct_ibz_set(&(mat[1][3]),1);
    nct_ibz_set(&(mat[2][2]),5);
    nct_ibz_set(&(mat[2][3]),4);
    nct_ibz_set(&(mat[3][3]),4);
    res = res || (nct_ibz_mat_4x4_is_hnf(&mat));


    nct_ibz_set(&(mat[0][0]),7);
    nct_ibz_set(&(mat[0][1]),0);
    nct_ibz_set(&(mat[0][2]),5);
    nct_ibz_set(&(mat[0][3]),4);
    nct_ibz_set(&(mat[1][0]),2);
    nct_ibz_set(&(mat[1][1]),3);
    nct_ibz_set(&(mat[1][2]),-1);
    nct_ibz_set(&(mat[1][3]),7);
    nct_ibz_set(&(mat[2][2]),0);
    nct_ibz_set(&(mat[2][3]),0);
    nct_ibz_set(&(mat[3][3]),4);
    res = res || (nct_ibz_mat_4x4_is_hnf(&mat));

    if (res != 0){
        printf("Quaternion unit test dim4_nct_ibz_mat_4x4_is_hnf failed\n");
    }
    nct_ibz_mat_4x4_finalize(&mat);
    return(res);
}

//void nct_ibz_mat_4xn_hnf_mod_core(nct_ibz_mat_4x4_t *hnf, int generator_number, const nct_ibz_vec_4_t (*generators)[generator_number], const nct_ibz_t *mod);
int nct_quat_test_dim4_nct_ibz_mat_4xn_hnf_mod_core(){
    int res = 0;
    nct_ibz_t det;
    nct_ibz_vec_4_t generators[8];
    nct_ibz_mat_4x4_t hnf,cmp;
    nct_ibz_init(&det);
    nct_ibz_mat_4x4_init(&hnf);
    nct_ibz_mat_4x4_init(&cmp);
    for(int i = 0; i <8; i++)
        nct_ibz_vec_4_init(&(generators[i]));

    nct_ibz_set(&(generators[2][0]),2);
    nct_ibz_set(&(generators[3][1]),3);
    nct_ibz_set(&(generators[4][0]),4);
    nct_ibz_set(&(generators[2][3]),5);
    nct_ibz_set(&(generators[7][3]),6);
    nct_ibz_set(&(generators[7][1]),7);
    nct_ibz_set(&(generators[3][1]),8);
    nct_ibz_set(&(generators[1][1]),9);
    nct_ibz_set(&(generators[6][0]),10);
    nct_ibz_set(&(generators[5][0]),11);
    nct_ibz_set(&(generators[0][0]),12);
    nct_ibz_set(&(generators[5][2]),1);
    nct_ibz_set(&(generators[0][2]),2);
    nct_ibz_set(&det,4);
    nct_ibz_mat_4xn_hnf_mod_core(&hnf,8,&generators,&det);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&hnf));

    // test equality of result to a known hnf
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 8; j++){
            nct_ibz_set(&(generators[j][i]),0);
        }
    }
    nct_ibz_set(&(generators[0][0]),4);
    nct_ibz_set(&(generators[2][0]),3);
    nct_ibz_set(&(generators[4][0]),1);
    nct_ibz_set(&(generators[7][0]),-1);
    nct_ibz_set(&(generators[1][1]),5);
    nct_ibz_set(&(generators[5][1]),-2);
    nct_ibz_set(&(generators[2][2]),3);
    nct_ibz_set(&(generators[6][2]),1);
    nct_ibz_set(&(generators[5][2]),1);
    nct_ibz_set(&(generators[3][3]),7);
    nct_ibz_set(&(generators[7][3]),-3);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            nct_ibz_set(&(cmp[i][j]),0);
        }
        nct_ibz_set(&(cmp[i][i]),1);
    }
    nct_ibz_set(&det,1);
    nct_ibz_mat_4xn_hnf_mod_core(&hnf,8,&generators,&det);
    res = res || (nct_ibz_mat_4x4_equal(&cmp,&hnf));

    // test known hnf encountered in
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 8; j++){
            nct_ibz_set(&(generators[j][i]),0);
        }
    }
    nct_ibz_set(&(generators[4][0]),438);
    nct_ibz_set(&(generators[4][1]),400);
    nct_ibz_set(&(generators[4][2]),156);
    nct_ibz_set(&(generators[4][3]),-2);
    nct_ibz_set(&(generators[5][0]),-400);
    nct_ibz_set(&(generators[5][1]),438);
    nct_ibz_set(&(generators[5][2]),2);
    nct_ibz_set(&(generators[5][3]),156);
    nct_ibz_set(&(generators[6][0]),-28826);
    nct_ibz_set(&(generators[6][1]),-148);
    nct_ibz_set(&(generators[6][2]),220);
    nct_ibz_set(&(generators[6][3]),-122);
    nct_ibz_set(&(generators[7][0]),586);
    nct_ibz_set(&(generators[7][1]),-28426);
    nct_ibz_set(&(generators[7][2]),278);
    nct_ibz_set(&(generators[7][3]),218);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            nct_ibz_set(&(cmp[i][j]),0);
        }
    }
    nct_ibz_set(&(cmp[0][0]),2321156);
    nct_ibz_set(&(cmp[1][1]),2321156);
    nct_ibz_set(&(cmp[0][2]),620252);
    nct_ibz_set(&(cmp[1][2]),365058);
    nct_ibz_set(&(cmp[2][2]),2);
    nct_ibz_set(&(cmp[0][3]),1956098);
    nct_ibz_set(&(cmp[1][3]),620252);
    nct_ibz_set(&(cmp[3][3]),2);
    nct_ibz_set(&det,21551060705344);
    nct_ibz_mat_4xn_hnf_mod_core(&hnf,8,&generators,&det);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&hnf));
    res = res || (nct_ibz_mat_4x4_equal(&cmp,&hnf));
    if (res != 0){
        printf("Quaternion unit test dim4_nct_ibz_mat_4x8_hnf_mod_core failed\n");
    }
    nct_ibz_mat_4x4_finalize(&hnf);
    nct_ibz_mat_4x4_finalize(&cmp);
    nct_ibz_finalize(&det);
    for(int i = 0; i <8; i++)
        nct_ibz_vec_4_finalize(&(generators[i]));
    return(res);
}







//static inline void nct_quat_alg_init_set_ui(nct_quat_alg_t *alg, unsigned int p);
int nct_quat_test_init_set_ui(){
    int res = 0;
    int p = 5;
    nct_quat_alg_t alg;
    nct_ibz_mat_4x4_t cmp;
    nct_ibz_mat_4x4_init(&cmp);
    nct_quat_alg_init_set_ui(&alg, p);
    res = res || (0 != nct_quat_test_helper_nct_ibz_equal_i(&(alg.p),p));
    nct_ibz_mat_4x4_identity(&cmp);
    nct_ibz_set(&(cmp[2][2]),p);
    nct_ibz_set(&(cmp[3][3]),p);
    res = res || nct_ibz_mat_4x4_equal(&cmp,&(alg.gram));
    if (res != 0){
        printf("Quaternion unit test alg_init_set_ui failed\n");
    }
    nct_quat_alg_finalize(&alg);
    nct_ibz_mat_4x4_finalize(&cmp);
    return(res);
}


//void nct_quat_alg_equal_denom(nct_quat_alg_elem_t *res_a, nct_quat_alg_elem_t *res_b, const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b);
int nct_quat_test_alg_equal_denom(){
    int res = 0;
    nct_quat_alg_elem_t a, b, res_a, res_b, cmp_a, cmp_b;
    nct_quat_alg_elem_init(&a);
    nct_quat_alg_elem_init(&b);
    nct_quat_alg_elem_init(&res_a);
    nct_quat_alg_elem_init(&res_b);
    nct_quat_alg_elem_init(&cmp_a);
    nct_quat_alg_elem_init(&cmp_b);

    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),9);
    nct_ibz_set(&(b.coord[0]),-6);
    nct_ibz_set(&(b.coord[1]),2);
    nct_ibz_set(&(b.coord[2]),67);
    nct_ibz_set(&(b.coord[3]),-19);
    nct_ibz_set(&(b.denom),3);
    nct_ibz_set(&(cmp_a.coord[0]),-12);
    nct_ibz_set(&(cmp_a.coord[1]),0);
    nct_ibz_set(&(cmp_a.coord[2]),-7);
    nct_ibz_set(&(cmp_a.coord[3]),19);
    nct_ibz_set(&(cmp_a.denom),9);
    nct_ibz_set(&(cmp_b.coord[0]),-18);
    nct_ibz_set(&(cmp_b.coord[1]),6);
    nct_ibz_set(&(cmp_b.coord[2]),201);
    nct_ibz_set(&(cmp_b.coord[3]),-57);
    nct_ibz_set(&(cmp_b.denom),9);
    nct_quat_alg_equal_denom(&res_a,&res_b,&a,&b);
    res = res || nct_ibz_cmp(&(res_a.denom),&(cmp_a.denom));
    res = res || nct_ibz_cmp(&(res_b.denom),&(cmp_b.denom));
    res = res || nct_ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || nct_ibz_cmp(&(res_a.coord[i]),&(cmp_a.coord[i]));
        res = res || nct_ibz_cmp(&(res_b.coord[i]),&(cmp_b.coord[i]));
    }


    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),9);
    nct_ibz_set(&(b.coord[0]),-6);
    nct_ibz_set(&(b.coord[1]),2);
    nct_ibz_set(&(b.coord[2]),67);
    nct_ibz_set(&(b.coord[3]),-19);
    nct_ibz_set(&(b.denom),6);
    nct_ibz_set(&(cmp_a.coord[0]),-24);
    nct_ibz_set(&(cmp_a.coord[1]),0);
    nct_ibz_set(&(cmp_a.coord[2]),-14);
    nct_ibz_set(&(cmp_a.coord[3]),38);
    nct_ibz_set(&(cmp_a.denom),18);
    nct_ibz_set(&(cmp_b.coord[0]),-18);
    nct_ibz_set(&(cmp_b.coord[1]),6);
    nct_ibz_set(&(cmp_b.coord[2]),201);
    nct_ibz_set(&(cmp_b.coord[3]),-57);
    nct_ibz_set(&(cmp_b.denom),18);
    nct_quat_alg_equal_denom(&res_a,&res_b,&a,&b);
    res = res || nct_ibz_cmp(&(res_a.denom),&(cmp_a.denom));
    res = res || nct_ibz_cmp(&(res_b.denom),&(cmp_b.denom));
    res = res || nct_ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || nct_ibz_cmp(&(res_a.coord[i]),&(cmp_a.coord[i]));
        res = res || nct_ibz_cmp(&(res_b.coord[i]),&(cmp_b.coord[i]));
    }

    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),6);
    nct_ibz_set(&(b.coord[0]),-6);
    nct_ibz_set(&(b.coord[1]),2);
    nct_ibz_set(&(b.coord[2]),67);
    nct_ibz_set(&(b.coord[3]),-19);
    nct_ibz_set(&(b.denom),6);
    nct_ibz_set(&(cmp_a.coord[0]),-12);
    nct_ibz_set(&(cmp_a.coord[1]),0);
    nct_ibz_set(&(cmp_a.coord[2]),-7);
    nct_ibz_set(&(cmp_a.coord[3]),19);
    nct_ibz_set(&(cmp_a.denom),6);
    nct_ibz_set(&(cmp_b.coord[0]),-6);
    nct_ibz_set(&(cmp_b.coord[1]),2);
    nct_ibz_set(&(cmp_b.coord[2]),67);
    nct_ibz_set(&(cmp_b.coord[3]),-19);
    nct_ibz_set(&(cmp_b.denom),6);
    nct_quat_alg_equal_denom(&res_a,&res_b,&a,&b);
    res = res || nct_ibz_cmp(&(res_a.denom),&(cmp_a.denom));
    res = res || nct_ibz_cmp(&(res_b.denom),&(cmp_b.denom));
    res = res || nct_ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || nct_ibz_cmp(&(res_a.coord[i]),&(cmp_a.coord[i]));
        res = res || nct_ibz_cmp(&(res_b.coord[i]),&(cmp_b.coord[i]));
    }

    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),6);
    nct_ibz_set(&(cmp_a.coord[0]),-12);
    nct_ibz_set(&(cmp_a.coord[1]),0);
    nct_ibz_set(&(cmp_a.coord[2]),-7);
    nct_ibz_set(&(cmp_a.coord[3]),19);
    nct_ibz_set(&(cmp_a.denom),6);
    nct_ibz_set(&(cmp_b.coord[0]),-12);
    nct_ibz_set(&(cmp_b.coord[1]),0);
    nct_ibz_set(&(cmp_b.coord[2]),-7);
    nct_ibz_set(&(cmp_b.coord[3]),19);
    nct_ibz_set(&(cmp_b.denom),6);
    nct_quat_alg_equal_denom(&a,&b,&a,&a);
    res = res || nct_ibz_cmp(&(a.denom),&(a.denom));
    res = res || nct_ibz_cmp(&(b.denom),&(b.denom));
    res = res || nct_ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || nct_ibz_cmp(&(a.coord[i]),&(a.coord[i]));
        res = res || nct_ibz_cmp(&(b.coord[i]),&(b.coord[i]));
    }
    if (res != 0){
        printf("Quaternion unit test alg_equal_denom failed\n");
    }
    nct_quat_alg_elem_finalize(&a);
    nct_quat_alg_elem_finalize(&b);
    nct_quat_alg_elem_finalize(&res_a);
    nct_quat_alg_elem_finalize(&res_b);
    nct_quat_alg_elem_finalize(&cmp_a);
    nct_quat_alg_elem_finalize(&cmp_b);
    return(res);
}


//Tests of public functions

//void nct_quat_alg_sub(nct_quat_alg_elem_t *res, const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b);
int nct_quat_test_alg_sub(){
    int res = 0;
    nct_quat_alg_elem_t a, b, c, cmp;
    nct_quat_alg_elem_init(&a);
    nct_quat_alg_elem_init(&b);
    nct_quat_alg_elem_init(&c);
    nct_quat_alg_elem_init(&cmp);

    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),9);
    nct_ibz_set(&(b.coord[0]),-6);
    nct_ibz_set(&(b.coord[1]),2);
    nct_ibz_set(&(b.coord[2]),7);
    nct_ibz_set(&(b.coord[3]),-19);
    nct_ibz_set(&(b.denom),3);
    nct_ibz_set(&(cmp.coord[0]),-12-3*(-6));
    nct_ibz_set(&(cmp.coord[1]),-3*2);
    nct_ibz_set(&(cmp.coord[2]),-7-3*7);
    nct_ibz_set(&(cmp.coord[3]),19-3*(-19));
    nct_ibz_set(&(cmp.denom),9);
    nct_quat_alg_sub(&c, &a, &b);
    res = res || nct_ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){
        res = res || nct_ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }

    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),9);
    nct_ibz_set(&(b.coord[0]),-6);
    nct_ibz_set(&(b.coord[1]),2);
    nct_ibz_set(&(b.coord[2]),7);
    nct_ibz_set(&(b.coord[3]),-19);
    nct_ibz_set(&(b.denom),6);
    nct_ibz_set(&(cmp.coord[0]),-2*12-3*(-6));
    nct_ibz_set(&(cmp.coord[1]),-3*2);
    nct_ibz_set(&(cmp.coord[2]),-2*7-3*7);
    nct_ibz_set(&(cmp.coord[3]),2*19-3*(-19));
    nct_ibz_set(&(cmp.denom),18);
    nct_quat_alg_sub(&a, &a, &b);
    res = res || nct_ibz_cmp(&(a.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || nct_ibz_cmp(&(a.coord[i]),&(cmp.coord[i]));
    }

    nct_ibz_set(&(a.coord[0]),-12);
    nct_ibz_set(&(a.coord[1]),0);
    nct_ibz_set(&(a.coord[2]),-7);
    nct_ibz_set(&(a.coord[3]),19);
    nct_ibz_set(&(a.denom),9);
    nct_ibz_set(&(cmp.coord[0]),0);
    nct_ibz_set(&(cmp.coord[1]),0);
    nct_ibz_set(&(cmp.coord[2]),0);
    nct_ibz_set(&(cmp.coord[3]),0);
    nct_ibz_set(&(cmp.denom),9);
    nct_quat_alg_sub(&a, &a, &a);
    res = res || nct_ibz_cmp(&(a.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || nct_ibz_cmp(&(a.coord[i]),&(cmp.coord[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_sub failed\n");
    }
    nct_quat_alg_elem_finalize(&a);
    nct_quat_alg_elem_finalize(&b);
    nct_quat_alg_elem_finalize(&c);
    nct_quat_alg_elem_finalize(&cmp);
    return(res);
}

//void nct_quat_alg_normalize(nct_quat_alg_elem_t *x);
int nct_quat_test_alg_normalize(){
    int res = 0;
    nct_quat_alg_elem_t x, cmp;
    nct_ibz_t gcd;
    nct_quat_alg_elem_init(&x);
    nct_quat_alg_elem_init(&cmp);
    nct_ibz_init(&gcd);

    // sign change
    nct_ibz_set(&(x.coord[0]), -125);
    nct_ibz_set(&(x.coord[1]), 2);
    nct_ibz_set(&(x.coord[2]), 0);
    nct_ibz_set(&(x.coord[3]), -30);
    nct_ibz_set(&(x.denom), -25);
    nct_ibz_set(&(cmp.coord[0]), 125);
    nct_ibz_set(&(cmp.coord[1]), -2);
    nct_ibz_set(&(cmp.coord[2]), 0);
    nct_ibz_set(&(cmp.coord[3]), 30);
    nct_ibz_set(&(cmp.denom), 25);
    nct_quat_alg_normalize(&x);
    res = res || nct_ibz_cmp(&(x.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || nct_ibz_cmp(&(x.coord[i]),&(cmp.coord[i]));
    }
    //divide by gcd
    nct_ibz_set(&(x.coord[0]), -36);
    nct_ibz_set(&(x.coord[1]), 18);
    nct_ibz_set(&(x.coord[2]), 0);
    nct_ibz_set(&(x.coord[3]), -300);
    nct_ibz_set(&(x.denom), 48);
    nct_ibz_set(&(cmp.coord[0]), -6);
    nct_ibz_set(&(cmp.coord[1]), 3);
    nct_ibz_set(&(cmp.coord[2]), 0);
    nct_ibz_set(&(cmp.coord[3]), -50);
    nct_ibz_set(&(cmp.denom), 8);
    nct_quat_alg_normalize(&x);
    res = res || nct_ibz_cmp(&(x.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || nct_ibz_cmp(&(x.coord[i]),&(cmp.coord[i]));
    }
    //divide by gcd
    nct_ibz_set(&(x.coord[0]), -36);
    nct_ibz_set(&(x.coord[1]), 18);
    nct_ibz_set(&(x.coord[2]), 0);
    nct_ibz_set(&(x.coord[3]), -300);
    nct_ibz_set(&(x.denom), -6);
    nct_ibz_set(&(cmp.coord[0]), 6);
    nct_ibz_set(&(cmp.coord[1]), -3);
    nct_ibz_set(&(cmp.coord[2]), 0);
    nct_ibz_set(&(cmp.coord[3]), 50);
    nct_ibz_set(&(cmp.denom), 1);
    nct_quat_alg_normalize(&x);
    res = res || nct_ibz_cmp(&(x.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || nct_ibz_cmp(&(x.coord[i]),&(cmp.coord[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_normalize failed\n");
    }
    nct_quat_alg_elem_finalize(&x);
    nct_quat_alg_elem_finalize(&cmp);
    nct_ibz_finalize(&gcd);
    return(res);
}


//int nct_quat_alg_elem_equal(const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b);
int nct_quat_test_alg_elem_equal(){
    int res = 0;
    nct_quat_alg_elem_t a,b;
    nct_quat_alg_elem_init(&a);
    nct_quat_alg_elem_init(&b);
    nct_ibz_vec_4_set(&(a.coord),1,-3,-2,2);
    nct_ibz_set(&(a.denom),5);
    nct_ibz_vec_4_set(&(b.coord),3,-9,-6,6);
    nct_ibz_set(&(b.denom),15);
    res = res || !nct_quat_alg_elem_equal(&a,&b);
    res = res || !nct_quat_alg_elem_equal(&a,&a);
    res = res || !nct_quat_alg_elem_equal(&b,&b);
    nct_ibz_vec_4_set(&(a.coord),1,-3,-2,2);
    nct_ibz_set(&(a.denom),5);
    nct_ibz_vec_4_set(&(b.coord),3,-9,-6,3);
    nct_ibz_set(&(b.denom),15);
    res = res || nct_quat_alg_elem_equal(&a,&b);
    nct_ibz_vec_4_set(&(a.coord),5,-15,-10,10);
    nct_ibz_set(&(a.denom),25);
    nct_ibz_vec_4_set(&(b.coord),3,-9,-6,6);
    nct_ibz_set(&(b.denom),15);
    res = res || !nct_quat_alg_elem_equal(&a,&b);
    nct_ibz_vec_4_set(&(a.coord),5,-15,-10,10);
    nct_ibz_set(&(a.denom),25);
    nct_ibz_vec_4_set(&(b.coord),0,-9,-6,6);
    nct_ibz_set(&(b.denom),5);
    res = res || nct_quat_alg_elem_equal(&a,&b);
    if (res != 0){
        printf("Quaternion unit test alg_elem_equal failed\n");
    }
    nct_quat_alg_elem_finalize(&a);
    nct_quat_alg_elem_finalize(&b);
    return(res);
}

//int nct_quat_alg_elem_is_zero(const nct_quat_alg_elem_t *x);
int nct_quat_test_alg_elem_is_zero(){
    int res = 0;
    nct_quat_alg_elem_t x;
    nct_quat_alg_elem_init(&x);
    nct_ibz_set(&(x.denom),1);
    nct_ibz_set(&(x.coord[0]),0);
    nct_ibz_set(&(x.coord[1]),0);
    nct_ibz_set(&(x.coord[2]),0);
    nct_ibz_set(&(x.coord[3]),0);
    res = res | (1-nct_quat_alg_elem_is_zero(&x));
    nct_ibz_set(&(x.denom),56865);
    res = res | (1-nct_quat_alg_elem_is_zero(&x));
    nct_ibz_set(&(x.denom),0);
    // maybe failure should be accepted here, but according to doc, this is still 0
    res = res | (1-nct_quat_alg_elem_is_zero(&x));
    nct_ibz_set(&(x.coord[3]),1);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.denom),56865);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[3]),-1);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[2]),1);
    nct_ibz_set(&(x.coord[3]),0);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[2]),-20);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[1]),1);
    nct_ibz_set(&(x.coord[2]),0);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[1]),-50000);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[0]),1);
    nct_ibz_set(&(x.coord[1]),0);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[0]),-90000);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[0]),0);
    nct_ibz_set(&(x.coord[1]),-500);
    nct_ibz_set(&(x.coord[2]),20);
    nct_ibz_set(&(x.coord[3]),0);
    res = res | nct_quat_alg_elem_is_zero(&x);
    nct_ibz_set(&(x.coord[0]),19);
    nct_ibz_set(&(x.coord[1]),-500);
    nct_ibz_set(&(x.coord[2]),20);
    nct_ibz_set(&(x.coord[3]),-2);
    res = res | nct_quat_alg_elem_is_zero(&x);
    if (res != 0){
        printf("Quaternion unit test alg_elem_is_zero failed\n");
    }
    nct_quat_alg_elem_finalize(&x);
    return(res);
}

//void nct_quat_alg_elem_copy_ibz(nct_quat_alg_elem_t *elem, const nct_ibz_t *denom, const nct_ibz_t *coord0,const nct_ibz_t *coord1,const nct_ibz_t *coord2,const nct_ibz_t *coord3){
int nct_quat_test_alg_elem_copy_ibz(){
    int res = 0;
    nct_ibz_t a,b,c,d,q;
    nct_quat_alg_elem_t elem;
    nct_quat_alg_elem_init(&elem);
    nct_ibz_init(&a);
    nct_ibz_init(&b);
    nct_ibz_init(&c);
    nct_ibz_init(&d);
    nct_ibz_init(&q);
    nct_ibz_set(&a,1);
    nct_ibz_set(&b,2);
    nct_ibz_set(&c,3);
    nct_ibz_set(&d,4);
    nct_ibz_set(&q,5);
    nct_quat_alg_elem_copy_ibz(&elem,&q,&a,&b,&c,&d);
    res = res || nct_ibz_cmp(&(elem.coord[0]),&a);
    res = res || nct_ibz_cmp(&(elem.coord[1]),&b);
    res = res || nct_ibz_cmp(&(elem.coord[2]),&c);
    res = res || nct_ibz_cmp(&(elem.coord[3]),&d);
    res = res || nct_ibz_cmp(&(elem.denom),&q);

    if (res != 0){
        printf("Quaternion unit test alg_elem_copy_ibz failed\n");
    }
    nct_ibz_finalize(&a);
    nct_ibz_finalize(&b);
    nct_ibz_finalize(&c);
    nct_ibz_finalize(&d);
    nct_ibz_finalize(&q);
    nct_quat_alg_elem_finalize(&elem);
    return(res);
}


//void nct_quat_alg_elem_set(nct_quat_alg_elem_t *elem, int64_t denom, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3)
int nct_quat_test_alg_elem_set(){
    int res = 0;
    nct_quat_alg_elem_t elem;
    nct_quat_alg_elem_init(&elem);
    nct_quat_alg_elem_set(&elem,5,1,2,3,4);
    res = res || (0 != nct_quat_test_helper_nct_ibz_equal_i(&(elem.coord[0]),1));
    res = res || (0 != nct_quat_test_helper_nct_ibz_equal_i(&(elem.coord[1]),2));
    res = res || (0 != nct_quat_test_helper_nct_ibz_equal_i(&(elem.coord[2]),3));
    res = res || (0 != nct_quat_test_helper_nct_ibz_equal_i(&(elem.coord[3]),4));
    res = res || (0 != nct_quat_test_helper_nct_ibz_equal_i(&(elem.denom),5));

    if (res != 0){
        printf("Quaternion unit test alg_elem_set failed\n");
    }
    nct_quat_alg_elem_finalize(&elem);
    return(res);
}

//int nct_quat_lattice_equal(const nct_quat_lattice_t *lat1, const nct_quat_lattice_t *lat2);
int nct_quat_test_lattice_equal(){
    int res = 0;
    nct_quat_lattice_t lat, cmp;
    nct_quat_lattice_init(&lat);
    nct_quat_lattice_init(&cmp);

    nct_ibz_mat_4x4_identity(&(lat.basis));
    nct_ibz_mat_4x4_identity(&(cmp.basis));
    res = res || !nct_quat_lattice_equal(&lat,&cmp);
    nct_ibz_set(&(lat.denom),5);
    nct_ibz_set(&(cmp.denom),4);
    res = res || nct_quat_lattice_equal(&lat,&cmp);
    nct_ibz_set(&(lat.denom),1);
    nct_ibz_set(&(cmp.denom),-1);
    res = res || !nct_quat_lattice_equal(&lat,&cmp);
    nct_ibz_set(&(lat.denom),3);
    nct_ibz_set(&(cmp.denom),3);
    res = res || !nct_quat_lattice_equal(&lat,&cmp);
    nct_ibz_set(&(lat.basis[0][0]),1);
    nct_ibz_set(&(lat.basis[0][3]),-1);
    nct_ibz_set(&(lat.basis[1][1]),-2);
    nct_ibz_set(&(lat.basis[2][2]),1);
    nct_ibz_set(&(lat.basis[2][1]),1);
    nct_ibz_set(&(lat.basis[3][3]),-3);
    nct_ibz_set(&(lat.denom),6);
    nct_quat_lattice_hnf(&lat);
    nct_ibz_mat_4x4_copy(&(cmp.basis), &(lat.basis));
    nct_ibz_set(&(cmp.denom),6);
    res = res || !nct_quat_lattice_equal(&lat,&cmp);
    nct_ibz_set(&(cmp.denom),-7);
    res = res || nct_quat_lattice_equal(&lat,&cmp);
    nct_ibz_set(&(cmp.denom),6);
    nct_ibz_set(&(cmp.basis[3][3]),165);
    res = res || nct_quat_lattice_equal(&lat,&cmp);

    if (res != 0){
        printf("Quaternion unit test lattice_equal failed\n");
    }
    nct_quat_lattice_finalize(&lat);
    nct_quat_lattice_finalize(&cmp);
    return(res);
}

//int nct_quat_lattice_inclusion(const nct_quat_lattice_t *sublat, const nct_quat_lattice_t *overlat)
int nct_quat_test_lattice_inclusion(){
    int res = 0;
    nct_quat_lattice_t lat, cmp;
    nct_quat_lattice_init(&lat);
    nct_quat_lattice_init(&cmp);

    nct_ibz_mat_4x4_identity(&(lat.basis));
    nct_ibz_mat_4x4_identity(&(cmp.basis));
    res = res || !nct_quat_lattice_inclusion(&lat,&cmp);
    nct_ibz_set(&(lat.denom),5);
    nct_ibz_set(&(cmp.denom),4);
    res = res || nct_quat_lattice_inclusion(&lat,&cmp);
    nct_ibz_set(&(lat.denom),1);
    nct_ibz_set(&(cmp.denom),3);
    res = res || !nct_quat_lattice_inclusion(&lat,&cmp);
    nct_ibz_set(&(lat.denom),3);
    nct_ibz_set(&(cmp.denom),3);
    res = res || !nct_quat_lattice_inclusion(&lat,&cmp);
    nct_ibz_set(&(lat.basis[0][0]),1);
    nct_ibz_set(&(lat.basis[0][3]),-1);
    nct_ibz_set(&(lat.basis[1][1]),-2);
    nct_ibz_set(&(lat.basis[2][2]),1);
    nct_ibz_set(&(lat.basis[2][1]),1);
    nct_ibz_set(&(lat.basis[3][3]),-3);
    nct_ibz_set(&(lat.denom),6);
    nct_quat_lattice_hnf(&lat);
    nct_ibz_mat_4x4_copy(&(cmp.basis), &(lat.basis));
    nct_ibz_set(&(cmp.denom),6);
    res = res || !nct_quat_lattice_inclusion(&lat,&cmp);
    nct_ibz_set(&(cmp.denom),12);
    res = res || !nct_quat_lattice_inclusion(&lat,&cmp);
    nct_ibz_set(&(cmp.denom),6);
    nct_ibz_set(&(cmp.basis[3][3]),165);
    res = res || nct_quat_lattice_inclusion(&lat,&cmp);

    if (res != 0){
        printf("Quaternion unit test lattice_inclusion failed\n");
    }
    nct_quat_lattice_finalize(&lat);
    nct_quat_lattice_finalize(&cmp);
    return(res);
}

//void nct_quat_lattice_reduce_denom(nct_quat_lattice_t *reduced, const nct_quat_lattice_t *lat);
int nct_quat_test_lattice_reduce_denom(){
    int res = 0;
    int s;
    nct_quat_lattice_t red, lat, cmp;
    nct_quat_lattice_init(&red);
    nct_quat_lattice_init(&cmp);
    nct_quat_lattice_init(&lat);

    s = 15;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_set(&(lat.basis[i][j]),(i+j)*s);
            nct_ibz_set(&(cmp.basis[i][j]),(i+j));
        }
    }
    nct_ibz_set(&(lat.denom),4*s);
    nct_ibz_set(&(cmp.denom), 4);

    nct_quat_lattice_reduce_denom(&red,&lat);
    res = res || (nct_ibz_mat_4x4_equal(&(red.basis),&(cmp.basis)));
    res = res || nct_ibz_cmp(&(red.denom),&(cmp.denom));


    nct_quat_lattice_reduce_denom(&lat,&lat);
    res = res || (nct_ibz_mat_4x4_equal(&(lat.basis),&(cmp.basis)));
    res = res || nct_ibz_cmp(&(lat.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test lattice_reduce_denom failed\n");
    }
    nct_quat_lattice_finalize(&red);
    nct_quat_lattice_finalize(&cmp);
    nct_quat_lattice_finalize(&lat);
    return(res);
}

//void nct_quat_lattice_add(nct_quat_lattice_t *res, const nct_quat_lattice_t *lat1, const nct_quat_lattice_t *lat2);
int nct_quat_test_lattice_add(){
    int res = 0;
    nct_quat_lattice_t lat1, lat2, cmp, sum;
    nct_quat_lattice_init(&lat1);
    nct_quat_lattice_init(&lat2);
    nct_quat_lattice_init(&sum);
    nct_quat_lattice_init(&cmp);
    nct_ibz_mat_4x4_zero(&(lat1.basis));
    nct_ibz_mat_4x4_zero(&(lat2.basis));
    nct_ibz_mat_4x4_zero(&(cmp.basis));
    nct_ibz_set(&(lat1.basis[0][0]),44);
    nct_ibz_set(&(lat1.basis[0][2]),3);
    nct_ibz_set(&(lat1.basis[0][3]),32);
    nct_ibz_set(&(lat2.basis[0][0]),1);
    nct_ibz_set(&(cmp.basis[0][0]),2);
    nct_ibz_set(&(cmp.basis[0][2]),1);
    nct_ibz_set(&(lat1.basis[1][1]),5);
    nct_ibz_set(&(lat2.basis[1][1]),2);
    nct_ibz_set(&(cmp.basis[1][1]),1);
    nct_ibz_set(&(lat1.basis[2][2]),3);
    nct_ibz_set(&(lat2.basis[2][2]),1);
    nct_ibz_set(&(cmp.basis[2][2]),1);
    nct_ibz_set(&(lat1.basis[3][3]),1);
    nct_ibz_set(&(lat2.basis[3][3]),3);
    nct_ibz_set(&(cmp.basis[3][3]),3);
    nct_ibz_set(&(lat1.denom),4);
    nct_ibz_set(&(lat2.denom),6);
    nct_ibz_set(&(cmp.denom),12);

    nct_quat_lattice_add(&sum,&lat1,&lat2);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&(sum.basis)));
    res = res || !nct_quat_lattice_equal(&sum,&cmp);

    // same lattices but not under hnf
    nct_ibz_mat_4x4_zero(&(lat1.basis));
    nct_ibz_mat_4x4_zero(&(lat2.basis));
    nct_ibz_set(&(lat1.basis[0][0]),4);
    nct_ibz_set(&(lat1.basis[0][2]),3);
    nct_ibz_set(&(lat2.basis[0][0]),1);
    nct_ibz_set(&(lat2.basis[0][3]),-1);
    nct_ibz_set(&(lat1.basis[1][1]),5);
    nct_ibz_set(&(lat2.basis[1][1]),-2);
    nct_ibz_set(&(lat1.basis[2][2]),3);
    nct_ibz_set(&(lat2.basis[2][2]),1);
    nct_ibz_set(&(lat2.basis[2][1]),1);
    nct_ibz_set(&(lat1.basis[3][3]),7);
    nct_ibz_set(&(lat2.basis[3][3]),-3);
    nct_ibz_set(&(lat1.denom),4);
    nct_ibz_set(&(lat2.denom),6);

    nct_quat_lattice_add(&sum,&lat1,&lat2);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&(sum.basis)));
    res = res || !nct_quat_lattice_equal(&sum,&cmp);

    //double in place gives hnf
    nct_ibz_mat_4x4_copy(&(cmp.basis),&lat2.basis);
    nct_ibz_copy(&(cmp.denom),&(lat2.denom));
    nct_quat_lattice_hnf(&cmp);
    nct_quat_lattice_add(&lat2,&lat2,&lat2);
    res = res || (!nct_ibz_mat_4x4_is_hnf(&(lat2.basis)));
    res = res || !nct_quat_lattice_equal(&lat2,&cmp);

    if (res != 0){
        printf("Quaternion unit test lattice_add failed\n");
    }
    nct_quat_lattice_finalize(&lat1);
    nct_quat_lattice_finalize(&lat2);
    nct_quat_lattice_finalize(&sum);
    nct_quat_lattice_finalize(&cmp);
    return(res);
}

//int nct_quat_lattice_contains_without_alg(nct_ibz_vec_4_t *coord, const nct_quat_lattice_t *lat, const nct_quat_alg_elem_t *x);
int nct_quat_test_lattice_contains_without_alg(){
    int res = 0;
    nct_quat_alg_elem_t x;
    nct_ibz_vec_4_t coord, cmp;
    nct_quat_lattice_t lat;
    nct_quat_alg_elem_init(&x);
    nct_ibz_vec_4_init(&coord);
    nct_ibz_vec_4_init(&cmp);
    nct_quat_lattice_init(&lat);

    // lattice 1
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_set(&(lat.basis[i][j]),0);
        }
    }
    nct_ibz_set(&(lat.basis[0][0]),4);
    nct_ibz_set(&(lat.basis[0][2]),3);
    nct_ibz_set(&(lat.basis[1][1]),5);
    nct_ibz_set(&(lat.basis[2][2]),3);
    nct_ibz_set(&(lat.basis[3][3]),7);
    nct_ibz_set(&(lat.denom),4);

    // x 1, should fail
    nct_ibz_set(&(x.denom),3);
    nct_ibz_set(&(x.coord[0]),1);
    nct_ibz_set(&(x.coord[1]),-2);
    nct_ibz_set(&(x.coord[2]),26);
    nct_ibz_set(&(x.coord[3]),9);

    res = res || nct_quat_lattice_contains_without_alg(&coord,&lat,&x);
    // again, but with NULL
    res = res || nct_quat_lattice_contains_without_alg(NULL,&lat,&x);


    // lattice 2
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_set(&(lat.basis[i][j]),0);
        }
    }
    nct_ibz_set(&(lat.basis[0][0]),1);
    nct_ibz_set(&(lat.basis[0][3]),-1);
    nct_ibz_set(&(lat.basis[1][1]),-2);
    nct_ibz_set(&(lat.basis[2][2]),1);
    nct_ibz_set(&(lat.basis[2][1]),1);
    nct_ibz_set(&(lat.basis[3][3]),-3);
    nct_ibz_set(&(lat.denom),6);
    nct_quat_lattice_hnf(&lat);
    // x 1, should succeed
    nct_ibz_set(&(x.denom),3);
    nct_ibz_set(&(x.coord[0]),1);
    nct_ibz_set(&(x.coord[1]),-2);
    nct_ibz_set(&(x.coord[2]),26);
    nct_ibz_set(&(x.coord[3]),9);
    nct_ibz_set(&(cmp[0]),2);
    nct_ibz_set(&(cmp[1]),-2);
    nct_ibz_set(&(cmp[2]),52);
    nct_ibz_set(&(cmp[3]),6);

    res = res || (0==nct_quat_lattice_contains_without_alg(&coord,&lat,&x));

    res = res || nct_ibz_cmp(&(coord[0]),&(cmp[0]));
    res = res || nct_ibz_cmp(&(coord[1]),&(cmp[1]));
    res = res || nct_ibz_cmp(&(coord[2]),&(cmp[2]));
    res = res || nct_ibz_cmp(&(coord[3]),&(cmp[3]));
    // again, but with NULL
    res = res || (0==nct_quat_lattice_contains_without_alg(NULL,&lat,&x));


    if (res != 0){
        printf("Quaternion unit test lattice_contains_without_alg failed\n");
    }
    nct_quat_alg_elem_finalize(&x);
    nct_ibz_vec_4_finalize(&coord);
    nct_ibz_vec_4_finalize(&cmp);
    nct_quat_lattice_finalize(&lat);
    return(res);
}


//void nct_quat_lattice_index(nct_ibz_t *index, const nct_quat_lattice_t *sublat, const nct_quat_lattice_t *overlat);
int nct_quat_test_lattice_index(){
    int res = 0;
    nct_quat_lattice_t sublat,overlat;
    nct_ibz_t index;
    nct_ibz_init(&index);
    nct_quat_lattice_init(&sublat);
    nct_quat_lattice_init(&overlat);

    nct_ibz_mat_4x4_zero(&(sublat.basis));
    nct_ibz_mat_4x4_identity(&(overlat.basis));
    nct_ibz_set(&(overlat.denom),2);
    nct_ibz_set(&(sublat.basis[0][0]),2);
    nct_ibz_set(&(sublat.basis[0][1]),0);
    nct_ibz_set(&(sublat.basis[0][2]),1);
    nct_ibz_set(&(sublat.basis[0][3]),0);
    nct_ibz_set(&(sublat.basis[1][0]),0);
    nct_ibz_set(&(sublat.basis[1][1]),4);
    nct_ibz_set(&(sublat.basis[1][2]),2);
    nct_ibz_set(&(sublat.basis[1][3]),3);
    nct_ibz_set(&(sublat.basis[2][0]),0);
    nct_ibz_set(&(sublat.basis[2][1]),0);
    nct_ibz_set(&(sublat.basis[2][2]),1);
    nct_ibz_set(&(sublat.basis[2][3]),0);
    nct_ibz_set(&(sublat.basis[3][0]),0);
    nct_ibz_set(&(sublat.basis[3][1]),0);
    nct_ibz_set(&(sublat.basis[3][2]),0);
    nct_ibz_set(&(sublat.basis[3][3]),1);
    nct_ibz_set(&(sublat.denom),2);
    nct_quat_lattice_index(&index,&sublat,&overlat);

    res = res || (0!=nct_quat_test_helper_nct_ibz_equal_i(&index,8));

    if (res != 0){
        printf("Quaternion unit test lattice_index failed\n");
    }
    nct_quat_lattice_finalize(&sublat);
    nct_quat_lattice_finalize(&overlat);
    nct_ibz_finalize(&index);
    return(res);
}

//void nct_quat_lattice_hnf(nct_quat_lattice_t *lat);
int nct_quat_test_lattice_hnf(){
    int res = 0;
    nct_quat_lattice_t lat, cmp;
    nct_quat_lattice_init(&lat);
    nct_quat_lattice_init(&cmp);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_set(&(lat.basis[i][j]),0);
            nct_ibz_set(&(cmp.basis[i][j]),0);
        }
    }
    nct_ibz_set(&(lat.basis[0][0]),1);
    nct_ibz_set(&(lat.basis[0][3]),-1);
    nct_ibz_set(&(lat.basis[1][1]),-2);
    nct_ibz_set(&(lat.basis[2][2]),1);
    nct_ibz_set(&(lat.basis[2][1]),1);
    nct_ibz_set(&(lat.basis[3][3]),-3);
    nct_ibz_set(&(cmp.basis[0][0]),1);
    nct_ibz_set(&(cmp.basis[1][1]),2);
    nct_ibz_set(&(cmp.basis[2][2]),1);
    nct_ibz_set(&(cmp.basis[3][3]),3);
    nct_ibz_set(&(cmp.denom),6);
    nct_ibz_set(&(lat.denom),6);

    nct_quat_lattice_hnf(&lat);
    res = res || (nct_ibz_mat_4x4_equal(&(lat.basis),&(cmp.basis)));
    res = res || nct_ibz_cmp(&(lat.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test lattice_hnf failed\n");
    }
    nct_quat_lattice_finalize(&lat);
    nct_quat_lattice_finalize(&cmp);
    return(res);
}


// run all algebra tests
int nct_quat_test_algebra(){
    int res = 0;
    res = res | nct_quat_test_init_set_ui();
    res = res | nct_quat_test_alg_sub();
    res = res | nct_quat_test_alg_equal_denom();
    res = res | nct_quat_test_alg_normalize();
    res = res | nct_quat_test_alg_elem_equal();
    res = res | nct_quat_test_alg_elem_is_zero();
    res = res | nct_quat_test_alg_elem_copy_ibz();
    res = res | nct_quat_test_alg_elem_set();
    return(res);
}

// run all lattice tests
int nct_quat_test_lattice(){
    int res = 0;
    res = res | nct_quat_test_lattice_equal();
    res = res | nct_quat_test_lattice_inclusion();
    res = res | nct_quat_test_lattice_reduce_denom();
    res = res | nct_quat_test_lattice_add();
    res = res | nct_quat_test_lattice_contains_without_alg();
    res = res | nct_quat_test_lattice_index();
    res = res | nct_quat_test_lattice_hnf();
    return(res);
}

int tests_helper_tests_all(){
    int res = 0;
    printf("Run test helper tests\n");
    res = res | nct_quat_test_xgcd_non_ct();
    res = res | nct_quat_test_dim4_nct_ibz_vec_4_scalar_div();
    res = res | nct_quat_test_dim4_nct_ibz_mat_4x4_scalar_div();
    res = res | nct_quat_test_integer_nct_ibz_centered_mod();
    res = res | nct_quat_test_integer_nct_ibz_conditional_assign();
    res = res | nct_quat_test_integer_nct_ibz_rounded_div();
    res = res | nct_quat_test_integer_nct_ibz_xgcd_with_u_not_0();
    res = res | nct_quat_test_integer_mod_not_zero();
    res = res | nct_quat_test_dim4_is_hnf();
    res = res | nct_quat_test_dim4_nct_ibz_vec_4_linear_combination();
    res = res | nct_quat_test_dim4_nct_ibz_mat_4xn_hnf_mod_core();
    res = res | nct_quat_test_algebra();
    res = res | nct_quat_test_lattice();
    if(res !=0){
        printf("Some tests failed\n");
    } else{
        printf("All passed\n");
    }
    return(res);
}
