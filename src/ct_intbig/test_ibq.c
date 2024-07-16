#include "test_ct_intbig.h"


int test_ibq_zero_num_denom(){
    int res = 0;
    ibq_t z;
    ibz_t part,cmp;
    ibz_init(&part);
    ibz_init(&cmp);
    ibq_init(&z);

    z[1][0] = (int64_t) 1;
    res = res | !ibq_is_zero(&z);
    z[1][0] = (int64_t) 5;
    res = res | !ibq_is_zero(&z);
    z[1][LIMBNUM-1] = (int64_t) -1;
    cmp[LIMBNUM-1]=(int64_t) -1;
    cmp[0]=5;
    ibq_denom(&part,&z);
    res = res | !(0==ibz_cmp(&part,&cmp));
    res = res | !ibq_is_zero(&z);
    z[1][0] = (int64_t) 0;
    res = res | !ibq_is_zero(&z);
    z[1][0] = (int64_t) 1;
    z[1][LIMBNUM-1] = (int64_t) 0;
    res = res | !ibq_is_zero(&z);


    z[1][0] = (int64_t) 5;
    z[1][1] = (int64_t) 29;
    res = res | !ibq_is_zero(&z);

    z[0][0] = (int64_t) 1024*1024*128;
    res = res | ibq_is_zero(&z);
    
    z[0][0] = (int64_t) 0;
    z[0][2] = (int64_t) 1024;
    cmp[0]=0;
    cmp[2]=1024;
    cmp[LIMBNUM-1]=0;
    ibq_num(&part,&z);
    res = res | !(0==ibz_cmp(&part,&cmp));
    res = res | ibq_is_zero(&z);

    if (res){
        printf("test_ibq_zero_num_denom failed\n");
    }
    ibz_finalize(&part);
    ibz_finalize(&cmp);
    ibq_finalize(&z);
    return(res);
}

int test_ibq_copy_set_swap_red(){
    int res = 0;
    ibz_t n,d, cmp;
    ibq_t a,b,c,z;
    ibz_init(&n);
    ibz_init(&d);
    ibz_init(&cmp);
    ibq_init(&a);
    ibq_init(&b);
    ibq_init(&c);
    ibq_init(&z);

    d[0]=2;
    d[2]=16;
    n[0]=35;
    n[1]=2;
    res = res | !ibq_set(&a,&n,&d);
    ibq_denom(&cmp,&a);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&a);
    res = res | !(0==ibz_cmp(&n,&cmp));
    ibq_copy(&b,&a);
    ibq_denom(&cmp,&b);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&b);
    res = res | !(0==ibz_cmp(&n,&cmp));
    b[0][0]=32;
    ibz_gcd(&cmp,&b[0],&b[1]);
    ibq_red(&b);
    d[0]=1;
    d[2]=8;
    n[0]=16;
    n[1]=1;
    ibq_denom(&cmp,&b);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&b);
    res = res | !(0==ibz_cmp(&n,&cmp));

    d[0]=2;
    d[2]=16;
    n[0]=32;
    n[1]=2;
    res = res | !ibq_set(&a,&n,&d);
    d[0]=1;
    d[2]=8;
    n[0]=16;
    n[1]=1;
    ibq_denom(&cmp,&a);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&a);
    res = res | !(0==ibz_cmp(&n,&cmp));
    b[1][0]=2;
    b[1][2]=16;
    b[0][0]=32;
    b[0][1]=2;
    ibq_copy(&b,&b);
    ibq_denom(&cmp,&b);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&b);
    res = res | !(0==ibz_cmp(&n,&cmp));

    n[1]=3;
    cmp[0]=0;
    cmp[1]=0;
    cmp[2]=0;
    res = res | ibq_set(&c,&n,&cmp);
    c[1][0]=1;

    ibq_swap(&c,&a);
    ibq_denom(&cmp,&c);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&c);
    n[1]=1;
    res = res | !(0==ibz_cmp(&n,&cmp));
    n[1] = 3;
    d[2] = 0;
    d[0] = 1;
    ibq_denom(&cmp,&a);
    res = res | !(0==ibz_cmp(&d,&cmp));
    ibq_num(&cmp,&a);
    res = res | !(0==ibz_cmp(&n,&cmp));

    //Must test red with negative inputs
    ibq_red(&b);
    ibq_num(&n,&b);
    ibz_neg(&n,&n);
    ibq_denom(&d,&b);
    ibq_set(&b,&n,&d);
    ibq_red(&b);
    res = res || (ibz_cmp(&d,&(b[1]))!=0);
    res = res || (ibz_cmp(&n,&(b[0]))!=0);

    if (res){
        printf("test_ibq_copy_set_swap_red failed\n");
    }
    ibz_finalize(&d);
    ibz_finalize(&cmp);
    ibz_finalize(&n);
    ibq_finalize(&c);
    ibq_finalize(&z);
    ibq_finalize(&a);
    ibq_finalize(&b);
    return(res);
}


int test_ibq_neg_abs(){
    int res = 0;
    ibq_t a,b,r;
    ibz_t d,n;
    ibq_init(&a);
    ibq_init(&b);
    ibq_init(&r);
    ibz_init(&d);
    ibz_init(&n);
    d[0] = 2041;
    d[1] = 987;
    n[1] = 2;
    n[2] = 56431;
    ibq_set(&a,&n,&d);
    ibq_copy(&b,&a);
    ibq_abs(&r,&b);
    ibq_red(&b);
    ibq_red(&r);
    ibq_denom(&d,&r);
    ibq_denom(&n,&b);
    res = res || (ibz_cmp(&d,&n)!=0);
    ibq_num(&n,&b);
    ibq_num(&d,&r);
    res = res || (ibz_cmp(&d,&n)!=0);

    ibq_neg(&a,&b);
    ibq_red(&b);
    ibq_red(&a);
    ibq_denom(&d,&a);
    ibq_denom(&n,&b);
    int t1 = (ibz_cmp(&d,&n)!=0);
    ibq_num(&n,&b);
    ibq_num(&d,&a);
    res = res || ((ibz_cmp(&d,&n)!=0) && t1);

    ibq_abs(&r,&a);
    ibq_red(&b);
    ibq_red(&r);
    ibq_denom(&d,&r);
    ibq_denom(&n,&b);
    res = res || (ibz_cmp(&d,&n)!=0);
    ibq_num(&n,&b);
    ibq_num(&d,&r);
    res = res || (ibz_cmp(&d,&n)!=0);

    if (res){
        printf("test_ibq_neg_abs failed\n");
    }
    ibz_finalize(&d);
    ibz_finalize(&n);
    ibq_finalize(&a);
    ibq_finalize(&b);
    ibq_finalize(&r);
    return(res);
}

int test_ibq_addition_substraction_cmp(){
    int res = 0;
    ibq_t a, b, sum, cmp;
    ibz_t n,d;
    ibq_init(&a);
    ibq_init(&b);
    ibq_init(&sum);
    ibq_init(&cmp);
    ibz_init(&d);
    ibz_init(&n);
    ibz_set(&d,1);
    ibq_set(&cmp,&n,&d);
    d[0] = 2;
    d[1] = 20093;
    n[3] = 20;
    n[0] = 6;
    ibq_set(&b,&n,&d);
    ibq_copy(&a,&b);
    res = res || ibq_is_zero(&b);
    res = res || (ibq_cmp(&b,&cmp)<=0);
    ibq_sub(&sum,&b,&b);
    res = res || !ibq_is_zero(&sum);
    res = res || (ibq_cmp(&sum,&cmp)!=0);
    ibq_sub(&sum,&b,&a);
    res = res || !ibq_is_zero(&sum);
    res = res || (ibq_cmp(&sum,&cmp)!=0);
    ibq_sub(&b,&cmp,&a);
    res = res || (ibq_cmp(&b,&cmp)>=0);
    ibq_add(&b,&b,&b);
    res = res || (ibq_cmp(&b,&cmp)>=0);
    ibq_add(&sum,&a,&b);
    res = res || (ibq_cmp(&sum,&cmp)>=0);
    res = res || (ibq_cmp(&b,&sum)>=0);
    ibq_add(&sum,&a,&sum);
    res = res || !ibq_is_zero(&sum);
    res = res || (ibq_cmp(&sum,&cmp)!=0);

    //tests with different denoms
    ibz_set(&(a[0]),-250);
    ibz_set(&(a[1]),210);
    ibz_set(&(b[0]),7);
    ibz_set(&(b[1]),50);
    ibz_set(&(cmp[0]),-1103);
    ibz_set(&(cmp[1]),1050);
    ibq_add(&sum,&a,&b);
    res = res || (ibq_cmp(&sum,&cmp)!=0);
    ibq_sub(&sum,&sum,&a);
    res = res || (ibq_cmp(&sum,&b)!=0);
    ibq_add(&sum,&b,&a);
    res = res || (ibq_cmp(&sum,&cmp)!=0);
    ibq_sub(&sum,&sum,&b);
    res = res || (ibq_cmp(&sum,&a)!=0);

    if (res){
        printf("test_ibq_addition_substraction_cmp failed\n");
    }
    ibz_finalize(&d);
    ibz_finalize(&n);
    ibq_finalize(&a);
    ibq_finalize(&b);
    ibq_finalize(&sum);
    ibq_finalize(&cmp);
    return(res);
}

int test_ibq_multiplication(){
    int res = 0;
    ibq_t a, b, prod, cmp;
    ibq_init(&a);
    ibq_init(&b);
    ibq_init(&cmp);
    ibq_init(&prod);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),-10);
    ibz_set(&(b[1]),1);
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),1);
    ibq_mul(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),-10);
    ibz_set(&(b[1]),3);
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),3);
    ibq_mul(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),-5);
    ibz_set(&(b[0]),10);
    ibz_set(&(b[1]),3);
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),3);
    ibq_mul(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),3);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),10);
    ibz_set(&(b[1]),9);
    ibz_set(&(cmp[0]),2);
    ibz_set(&(cmp[1]),3);
    ibq_mul(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);


    ibz_set(&(a[0]),-3);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),10);
    ibz_set(&(b[1]),-9);
    ibz_set(&(cmp[0]),2);
    ibz_set(&(cmp[1]),3);
    ibq_mul(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);
    if (res){
        printf("test_ibq_multiplication failed\n");
    }
    ibq_finalize(&a);
    ibq_finalize(&b);
    ibq_finalize(&cmp);
    ibq_finalize(&prod);
    return(res);
}

int test_ibq_division(){
    int res = 0;
    ibq_t a, b, prod, cmp;
    ibq_init(&a);
    ibq_init(&b);
    ibq_init(&cmp);
    ibq_init(&prod);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),1);
    ibz_set(&(b[1]),-10);
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),1);
    ibq_div(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),-3);
    ibz_set(&(b[1]),10);
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),3);
    ibq_div(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),-5);
    ibz_set(&(b[0]),3);
    ibz_set(&(b[1]),10);
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),3);
    ibq_div(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),3);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),9);
    ibz_set(&(b[1]),10);
    ibz_set(&(cmp[0]),2);
    ibz_set(&(cmp[1]),3);
    ibq_div(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    ibz_set(&(a[0]),-3);
    ibz_set(&(a[1]),5);
    ibz_set(&(b[0]),-9);
    ibz_set(&(b[1]),10);
    ibz_set(&(cmp[0]),2);
    ibz_set(&(cmp[1]),3);
    ibq_div(&prod,&a,&b);
    res = res || (ibq_cmp(&prod,&cmp)!=0);

    a[0][0]=2;
    a[0][1]=2312;
    a[0][02]=9987;
    a[1][0]=112;
    a[1][1]=21092;
    ibz_set(&(cmp[0]),1);
    ibz_set(&(cmp[1]),1);
    ibq_div(&a,&a,&a);
    res = res || (ibq_cmp(&a,&cmp)!=0);

    a[0][0]=2;
    a[0][1]=2312;
    a[0][02]=9987;
    a[1][0]=1123;
    a[1][1]=21092;

    b[0][0]=2;
    b[0][1]=220;
    b[0][02]=2672;
    b[1][0]=11872;
    b[1][1]=152;
    ibq_div(&prod,&a,&b);
    ibq_mul(&prod,&prod,&b);
    res = res || (ibq_cmp(&prod,&a)!=0);
    ibq_neg(&b,&b);
    ibq_div(&prod,&a,&b);
    ibq_mul(&prod,&b,&prod);
    res = res || (ibq_cmp(&prod,&a)!=0);
    ibq_neg(&a,&a);
    ibq_div(&prod,&a,&b);
    ibq_mul(&prod,&b,&prod);
    res = res || (ibq_cmp(&prod,&a)!=0);
    ibq_neg(&b,&b);
    ibq_div(&prod,&a,&b);
    ibq_mul(&prod,&b,&prod);
    res = res || (ibq_cmp(&prod,&a)!=0);

    if (res){
        printf("test_ibq_division failed\n");
    }
    ibq_finalize(&a);
    ibq_finalize(&b);
    ibq_finalize(&cmp);
    ibq_finalize(&prod);
    return(res);
}

int test_ibq_to_ibz() {
    int res = 0;
    int mynums[10] = {1, 8, 9, 5, 137137,-2,-2,2,3,0};
    int mydens[10] = {1, 4, 10, 4, 9876,-1,1,-1,-2,1};
    int close_ints[10] = {1, 2, 0, 1, 13,2,-2,-2,-1,0};

    ibz_t quo, rem, a, b, myint;

    ibz_init(&quo);
    ibz_init(&rem);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&myint);

    for (int i=0; i < 10; ++i) {
        ibz_set(&a, mynums[i]);
        ibz_set(&b, mydens[i]);
        ibz_set(&myint, close_ints[i]);
        ibz_div_floor(&quo, &rem, &a, &b);
        res = res || !(0 == ibz_cmp(&quo, &myint));
    }

    if (res) {
        printf("test_ibq_to_ibz failed\n");
    }

    ibz_finalize(&quo);
    ibz_finalize(&rem);
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&myint);

    return res;
}

int ibq_tests(){
    int res = 0;
    printf("Run ibq tests:\n");
    res = res | test_ibq_zero_num_denom();
    res = res | test_ibq_copy_set_swap_red();
    res = res | test_ibq_neg_abs();
    res = res | test_ibq_addition_substraction_cmp();
    res = res | test_ibq_multiplication();
    res = res | test_ibq_division();
    res = res | test_ibq_to_ibz();
    if(!res){
        printf("All passed\n");
    } else{
        printf("Some tests failed\n");
    }
    return res;
}