#include "test_ct_intbig.h"

int test_ibz_zero_cmp_init(){
    int res = 0;
    ibz_t z,a,b;
    ibz_init(&z);
    ibz_init(&a);
    ibz_init(&b);

    res = res | !ibz_is_zero(&z);

    z[0] = (int64_t) 1024*1024*128;
    res = res | ibz_is_zero(&z);
    
    z[0] = (int64_t) 0;
    z[2] = (int64_t) 1024;
    res = res | ibz_is_zero(&z); 

    a[2] = 1024;
    res = res | ibz_cmp(&a,&z);

    z[2] = (int64_t) 1024*1024*128;
    res = res | !(ibz_cmp(&z,&b)>0);

    z[2]=0;
    res = res | !(ibz_cmp(&a,&z)>0);

    b[2] = 1024;
    b[1] = 200;
    res = res | !(ibz_cmp(&a,&b)<0);


    a[0] = 5000;
    res = res | !(ibz_cmp(&a,&b)<0);

    for(int i =0; i<LIMBNUM;i++){
        a[i] = ((int64_t)-1) ^ a[i];
    }
    res = res | !(ibz_cmp(&a,&z)<0);

    for(int i =0; i<LIMBNUM;i++){
        b[i] = ((int64_t)-1) ^ b[i];
    }
    res = res | !(ibz_cmp(&a,&b)>0);


    if (res){
        printf("test_ibz_zero_cmp_init failed\n");
    }
    ibz_finalize(&z);
    ibz_finalize(&a);
    ibz_finalize(&b);
    return(res);
}

int test_ibz_copy_set(){
    int res = 0;
    ibz_t a,b,c, z;
    ibz_init(&z);
    ibz_init(&c);
    ibz_init(&a);
    ibz_init(&b);

    ibz_set(&a,50);
    ibz_set(&b,50);
    res = res || ibz_cmp(&a,&b);
    res = res || !(ibz_cmp(&a,&z)>0);

    ibz_set(&a, -5);
    res = res || !(ibz_cmp(&a,&z)<0);

    b[2]=2;
    ibz_copy(&b,&a);
    res = res || ibz_cmp(&a,&b);
    ibz_set(&c,1);
    ibz_copy(&a,&c);
    res = res || ibz_cmp(&a,&c);
    res = res || !ibz_cmp(&a,&b);

    // test sing sub for negative sets
    ibz_set(&b,4);
    ibz_sub(&c,&a,&b);
    ibz_set(&a,-3);
    res = res || ibz_cmp(&a,&c);

    ibz_set(&b,-1);
    for(int i = 1; i<LIMBNUM; i++){
        b[i] = 0;
    }
    ibz_add(&b,&b,&b);
    ibz_set_from_str(&a,"36893488147419103230",10);
    res = res || ibz_cmp(&a,&b);

    ibz_set(&b,10);
    ibz_set_from_str(&a,"10",10);
    res = res || ibz_cmp(&a,&b);

    if (res){
        printf("test_ibz_copy_set failed\n");
    }
    ibz_finalize(&c);
    ibz_finalize(&z);
    ibz_finalize(&a);
    ibz_finalize(&b);
    return(res);
}

int test_ibz_addition_substraction(){
    int res = 0;
    ibz_t a,b,c, z;
    ibz_init(&z);
    ibz_init(&c);
    ibz_init(&a);
    ibz_init(&b);

    ibz_set(&a, 1024*1024*256);
    ibz_set(&b, 1024*512*256);
    ibz_add(&b,&b,&b);
    res = res || ibz_cmp(&b,&a);


    ibz_set(&a, 1024*1024*256);
    ibz_set(&b, 1024*512*256);
    ibz_sub(&a,&a,&b);
    res = res || ibz_cmp(&b,&a);


    ibz_set(&a, 1024*1024*256);
    ibz_set(&b, 1024*512*256);
    ibz_sub(&c,&b,&a);
    ibz_sub(&b,&z,&b);
    res = res || ibz_cmp(&b,&c);

    ibz_add(&b,&b,&b);
    ibz_sub(&a,&z,&a);
    res = res || ibz_cmp(&b,&a);

    ibz_set(&c, 200);
    ibz_add(&b,&b,&c);
    ibz_set(&c, 220);
    ibz_add(&a,&a,&c);
    ibz_set(&c, -20);
    ibz_add(&a,&a,&c);
    res = res || ibz_cmp(&b,&a);


    ibz_set(&a, 1024*1024*256+2000);
    ibz_set(&b, 1024*512*256+1000);
    ibz_add(&b,&b,&b);
    ibz_add(&b,&b,&b);
    ibz_add(&b,&b,&b);
    ibz_add(&a,&a,&a);
    ibz_add(&a,&a,&a);
    res = res || ibz_cmp(&b,&a);

    ibz_set(&a,0);
    res = res || !ibz_is_zero(&a);
    res = res || ibz_cmp(&a,&z);

    if (res){
        printf("test_ibz_addition_substraction failed\n");
    }
    ibz_finalize(&c);
    ibz_finalize(&z);
    ibz_finalize(&a);
    ibz_finalize(&b);
    return(res);
}

int test_ibz_neg_abs(){
    int res = 0;
    ibz_t a, b, z;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&z);

    a[2]=100;
    a[0]=10;
    ibz_neg(&b,&a);
    ibz_add(&z,&b,&a);
    res = res || !ibz_is_zero(&z);

    ibz_abs(&b,&b);
    res = res || ibz_cmp(&a,&b);

    a[LIMBNUM-1] = -1;
    ibz_copy(&b,&a);
    ibz_neg(&a,&a);
    ibz_add(&z,&b,&a);
    res = res || !ibz_is_zero(&z);
    ibz_abs(&b,&b);
    res = res || ibz_cmp(&a,&b);

    if (res){
        printf("test_ibz_neg_abs failed\n");
    }
    ibz_finalize(&b);
    ibz_finalize(&a);
    ibz_finalize(&z);
    return(res);
}

int test_ibz_multiplication(){
    int res = 0;
    ibz_t a,b,c, z;
    ibz_init(&z);
    ibz_init(&c);
    ibz_init(&a);
    ibz_init(&b);

    a[0]=5;
    b[0]=7;
    b[4]= 12;
    b[2]=10000001;
    ibz_mul(&c,&a,&b);
    ibz_mul(&z,&b,&a);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_add(&c,&b,&b);
    ibz_add(&c,&c,&c);
    ibz_add(&c,&c,&b);
    res = res || !(0 == ibz_cmp(&c,&z));
    a[4]=100000112;
    ibz_set(&z,0);
    z[4]=a[4];
    ibz_mul(&z,&b,&z);
    ibz_add(&z,&c,&z);
    ibz_mul(&c,&a,&b);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_neg(&a,&a);
    ibz_neg(&b,&b);
    ibz_mul(&c,&a,&b);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_mul(&c,&b,&a);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_neg(&b,&b);
    ibz_neg(&z,&z);
    ibz_mul(&c,&a,&b);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_mul(&c,&b,&a);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_neg(&a,&a);
    ibz_neg(&b,&b);
    ibz_mul(&c,&a,&b);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_mul(&c,&b,&a);
    res = res || !(0 == ibz_cmp(&c,&z));
    ibz_copy(&c,&a);
    ibz_mul(&a,&b,&a);
    res = res || !(0 == ibz_cmp(&a,&z));
    ibz_copy(&a,&c);
    ibz_mul(&a,&a,&b);
    res = res || !(0 == ibz_cmp(&a,&z));
    ibz_set(&a,1);
    ibz_mul(&c,&a,&b);
    res = res || !(0 == ibz_cmp(&c,&b));
    ibz_mul(&c,&b,&a);
    res = res || !(0 == ibz_cmp(&c,&b));


    if (res){
        printf("test_ibz_multiplication failed\n");
    }
    ibz_finalize(&c);
    ibz_finalize(&z);
    ibz_finalize(&a);
    ibz_finalize(&b);
    return(res);
}

int test_ibz_division(){
    int res = 0;
    ibz_t a,b,c,d,e, z,q,r;
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&e);
    ibz_init(&z);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&q);
    ibz_init(&r);

    a[2] = 10;
    a[1] = 5;
    a[0] = 1098;
    b[0] = 7;
    b[2] = 600;
    b[3] = 10001;
    c[0] = 1;
    ibz_mul(&d,&a,&b);
    ibz_add(&d,&d,&c);
    ibz_div_floor(&q,&r,&d,&b);
    res = res || !(0 == ibz_cmp(&q,&a));
    res = res || !(0 == ibz_cmp(&c,&r));
    ibz_neg(&a,&a);
    ibz_mul(&d,&a,&b);
    ibz_add(&d,&d,&c);
    ibz_div_floor(&q,&r,&d,&b);
    res = res || !(0 == ibz_cmp(&q,&a));
    res = res || !(0 == ibz_cmp(&c,&r));

    ibz_set(&a,0);
    ibz_set(&b,0);
    ibz_set(&c,0);
    a[2] = 11271;
    a[1] = 5767658;
    a[0] = 1098;
    b[0] = 97;
    b[2] = 602290;
    b[3] = 9020031;
    c[0] = 40119222;
    c[2] = 1009087;
    ibz_mul(&d,&a,&b);
    ibz_add(&d,&d,&c);
    
    ibz_neg(&d,&d);
    ibz_neg(&b,&b);
    ibz_copy(&e,&d);
    ibz_div_floor(&q,&r,&d,&b);
    ibz_mul(&c,&q,&b);
    ibz_add(&c,&c,&r);
    res = res || !(0 == ibz_cmp(&e,&d));
    res = res || !(0 == ibz_cmp(&d,&c));
    res = res || !(ibz_cmp(&z,&r)<=0);
    ibz_abs(&e,&b);
    res = res || !(ibz_cmp(&r,&e)<0);
    ibz_neg(&d,&d);
    ibz_copy(&e,&d);
    ibz_div_floor(&q,&r,&d,&b);
    ibz_mul(&c,&q,&b);
    ibz_add(&c,&c,&r);
    res = res || !(0 == ibz_cmp(&d,&c));
    res = res || !(0 == ibz_cmp(&e,&c));
    res = res || !(ibz_cmp(&z,&r)<=0);
    ibz_abs(&e,&b);
    res = res || !(ibz_cmp(&r,&e)<0);
    ibz_neg(&b,&b);
    ibz_copy(&e,&d);
    ibz_div_floor(&q,&r,&d,&b);
    ibz_mul(&c,&q,&b);
    ibz_add(&c,&c,&r);
    res = res || !(0 == ibz_cmp(&d,&c));
    res = res || !(0 == ibz_cmp(&e,&c));
    res = res || !(ibz_cmp(&z,&r)<=0);
    ibz_abs(&e,&b);
    res = res || !(ibz_cmp(&r,&e)<0);
    ibz_neg(&d,&d);
    ibz_copy(&e,&d);
    ibz_div_floor(&q,&r,&d,&b);
    ibz_mul(&c,&q,&b);
    ibz_add(&c,&c,&r);
    res = res || !(0 == ibz_cmp(&d,&c));
    res = res || !(0 == ibz_cmp(&e,&c));
    res = res || !(ibz_cmp(&z,&r)<=0);
    ibz_abs(&e,&b);
    res = res || !(ibz_cmp(&r,&e)<0);

    ibz_set(&d,24);
    ibz_set(&b,6);
    ibz_copy(&e,&d);
    ibz_div_floor(&q,&r,&d,&b);
    ibz_mul(&c,&q,&b);
    ibz_add(&c,&c,&r);
    res = res || !(0 == ibz_cmp(&e,&d));
    res = res || !(0 == ibz_cmp(&d,&c));
    res = res || !(ibz_cmp(&z,&r)<=0);
    ibz_abs(&e,&b);
    res = res || !(ibz_cmp(&r,&e)<0);

    ibz_set(&d,2);
    ibz_set(&b,24);
    ibz_div_floor(&q,&r,&b,&d);
    ibz_mul(&c,&q,&d);
    res = res || !(0 == ibz_cmp(&b,&c));
    res = res || !(ibz_is_zero(&r));

    ibz_set(&d,2);
    ibz_set(&b,-24);
    ibz_div_floor(&q,&r,&b,&d);
    ibz_mul(&c,&q,&d);
    res = res || !(0 == ibz_cmp(&b,&c));
    res = res || !(ibz_is_zero(&r));

    if (res){
        printf("test_ibz_division failed\n");
    }
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&e);
    ibz_finalize(&z);
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&q);
    ibz_finalize(&r);
    return(res);
}

int test_ibz_gcd(){
    int res = 0;
    ibz_t a,b,c,d,e, z,q,r;
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&e);
    ibz_init(&z);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&q);
    ibz_init(&r);

    ibz_set(&a,24);
    ibz_set(&b,54);
    ibz_gcd(&d,&a,&b);
    ibz_set(&c,24);
    res = res || !(0 == ibz_cmp(&a,&c));
    ibz_set(&c,54);
    res = res || !(0 == ibz_cmp(&b,&c));
    ibz_set(&c,6);
    res = res || !(0 == ibz_cmp(&d,&c));
    ibz_div_floor(&q,&r,&a,&d);
    res = res || !(ibz_is_zero(&r));
    ibz_div_floor(&c,&r,&b,&d);
    res = res || !(ibz_is_zero(&r));
    ibz_gcd(&r,&q,&d);
    ibz_gcd(&z,&c,&d);
    ibz_gcd(&r,&r,&z);
    ibz_set(&z,1);
    res = res || !(0 == ibz_cmp(&r,&z));
    ibz_set(&z,0);

    a[2]=5556;
    a[1]=55561;
    a[0]=100;
    b[3]=1000864;
    b[1]=1276;
    b[0]=220;
    ibz_gcd(&d,&a,&b);
    ibz_div_floor(&q,&r,&a,&d);
    res = res || !(ibz_is_zero(&r));
    ibz_div_floor(&c,&r,&b,&d);
    res = res || !(ibz_is_zero(&r));
    ibz_gcd(&r,&q,&d);
    ibz_gcd(&z,&c,&d);
    ibz_gcd(&r,&r,&z);
    ibz_set(&z,1);
    res = res || !(0 == ibz_cmp(&r,&z));
    ibz_set(&z,0);

    ibz_set(&a,0);
    ibz_set(&b,0);
    a[0]=2;
    a[2]=16;
    b[0]=34;
    b[1]=2;
    ibz_gcd(&d,&a,&b);
    ibz_set(&r,6);
    ibz_abs(&d,&d);
    res = res | !(0==ibz_cmp(&r,&d));
    
    if (res){
        printf("test_ibz_gcd failed\n");
    }
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&e);
    ibz_finalize(&z);
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&q);
    ibz_finalize(&r);
    return(res);
}

int test_ibz_bitwise_and() {
    int res = 0;

    int mya[6] = {40, 99, 80, 87, 2, 3};
    int myb[6] = {32, 27, 60, 73, 1, 1};
    int myres[6] = {32, 3, 16, 65, 0, 1};

    ibz_t a, b, andres, band;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&andres);
    ibz_init(&band);

    for (int i = 0; i < 6; ++i) {
        ibz_set(&a, mya[i]);
        ibz_set(&b, myb[i]);
        ibz_set(&andres, myres[i]);

        ibz_bitwise_and(&band, &a, &b);

        res = res || !(0 == ibz_cmp(&band, &andres));
    }

    ibz_set_from_str(&a, "340282366920938463481821351505477763082", 10);
    ibz_set_from_str(&b, "18446744073709551626", 10);
    ibz_set_from_str(&andres, "18446744073709551626", 10);
    ibz_bitwise_and(&band, &a, &b);

    res = res || !(0 == ibz_cmp(&band, &andres));

    if (res) {
        printf("test_ibz_bitwise_and failed\n");
    }

    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&andres);
    ibz_finalize(&band);

    return res;
}

int test_ibz_internal_truncate() {
    int res = 0;
    ibz_t a, tres;

    int mya[5] = {143, 13, -65, 0,124};
    int myb[5] = {35, 34, 34, 28,8};
    int myres[5] = {143, 13, -65, 0,124};

    ibz_init(&a);
    ibz_init(&tres);

    for (int i = 0; i < 5; ++i) {
        ibz_set(&a, mya[i]);
        ibz_set(&tres, myres[i]);

        ibz_internal_truncate(&a, myb[i]);

        res = res || !(0 == ibz_cmp(&tres, &a));
    }

    if (res) {
        printf("Test ibz_internal_truncate failed\n");
    }

    ibz_finalize(&a);
    ibz_finalize(&tres);

    return res;
}

int test_ibz_internal_divsteps() {
    int res = 0;

    //int n[4] = {36, 36, 27, 27};
    //int t[4] = {47, 47, 35, 35};
    int a[4] = {1781, 1794, 137, 138};
    int b[4] = {13, 13, 139, 139};
    int myres[4] = {13, 13, 1, 1};

    ibz_t dres, result, f, g;

    ibz_init(&dres);
    ibz_init(&result);
    ibz_init(&f);
    ibz_init(&g);

    for (int i = 0; i < 4; ++i) {
        ibz_set(&f, a[i]);
        ibz_set(&g, b[i]);
        ibz_set(&result, myres[i]);

        ibz_internal_divsteps(&f, &g);
        ibz_abs(&f,&f);
        ibz_abs(&g,&g);

        res = res || !(0 == ibz_cmp(&result, &f));
    }

    if (res) {
        printf("Test_ibz_internal_divsteps failed\n");
    }

    ibz_finalize(&dres);
    ibz_finalize(&result);
    ibz_finalize(&f);
    ibz_finalize(&g);

    return res;
}

int test_gcd() {
    int res = 0;
    int a[3] = {137, 2, 1781};
    int b[3] = {139, 138, 137};
    int mygcd[3] = {1, 2, 13};

    ibz_t gcd, f, g, myres;
    ibz_init(&gcd);
    ibz_init(&f);
    ibz_init(&g);
    ibz_init(&myres);

    for (int i = 0; i < 3; ++i) {
        ibz_set(&f, a[i]);
        ibz_set(&g, b[i]);
        ibz_set(&myres, mygcd[i]);

        ibz_gcd(&gcd, &f, &g);

        res = res || !(0 == ibz_cmp(&gcd, &myres));
    }

    if (res) {
        printf("Test_gcd failed\n");
    }

    ibz_finalize(&gcd);
    ibz_finalize(&f);
    ibz_finalize(&g);

    return res;
}


int ibz_tests(){
    int res = 0;
    printf("Run ibz tests:\n");
    res = res | test_ibz_zero_cmp_init();
    res = res | test_ibz_copy_set();
    res = res | test_ibz_addition_substraction();
    res = res | test_ibz_neg_abs();
    res = res | test_ibz_multiplication();
    res = res | test_ibz_division();
    res = res | test_ibz_bitwise_and();
    res = res | test_ibz_internal_truncate();
    res = res | test_ibz_internal_divsteps();
    res = res | test_ibz_gcd();
    if(!res){
        printf("All passed\n");
    } else{
        printf("Some tests failed\n");
    }
    return res;
}