#include "test_bkz.h"
#include <stdlib.h>
#include <stdio.h>

//return 0 if equal, 1 otherwise
int quat_bkz_matrix_equal(const quat_bkz_matrix_t *a,const quat_bkz_matrix_t *b){
    int res = 0;
    res = res | ibz_mat_4x4_equal(&(a->basis),&(b->basis));
    res = res | ibz_mat_4x4_equal(&(a->G),&(b->G));
    res = res | ibq_mat_4x4_equal(&(a->mu),&(b->mu));
    res = res | ibq_mat_4x4_equal(&(a->r),&(b->r));
    return(res !=0);
}

//return 0 if equal to bkz matrix setup from is basis, 1 otherwise
int quat_test_helper_bkz_test_integrity(const quat_bkz_matrix_t *g, const quat_alg_t *alg){
    ibz_mat_4x4_t mat;
    quat_bkz_matrix_t b;
    ibz_mat_4x4_init(&mat);
    quat_bkz_matrix_init(&b);
    ibz_mat_4x4_transpose(&mat,&(g->basis));
    quat_bkz_matrix_set(&b,&mat,alg);
    int res = quat_bkz_matrix_equal(g,&b);
    ibz_mat_4x4_finalize(&mat);
    quat_bkz_matrix_finalize(&b);
    return(res);
}

//void quat_bkz_matrix_init(bkz_matrix_t *g);
//void quat_bkz_matrix_finalize(bkz_matrix_t *g);
int quat_test_bkz_matrix_finit(){
    int res = 0;
    ibz_t test;
    quat_bkz_matrix_t mat;
    ibz_init(&test);
    quat_bkz_matrix_init(&mat);
    ibz_set(&test,1);
    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            ibz_set(&(mat.basis[i][j]),i+4*j);
            ibz_set(&(mat.G[i][j]),i+4*j);
            ibq_set(&(mat.mu[i][j]),&(mat.G[i][j]),&test);
            ibq_set(&(mat.r[i][j]),&(mat.G[i][j]),&test);
        }
    }
    for(int i = 0; i<4; i++){
        for(int j = 0; j<4; j++){
            ibz_set(&test,i+4*j);
            res = res || (0!=ibz_cmp(&test,&(mat.basis[i][j])));
            res = res || (0!=ibz_cmp(&test,&(mat.G[i][j])));
            ibq_num(&test,&(mat.mu[i][j]));
            res = res || (0!=ibz_cmp(&test,&(mat.basis[i][j])));
            ibq_num(&test,&(mat.r[i][j]));
            res = res || (0!=ibz_cmp(&test,&(mat.basis[i][j])));
            ibq_denom(&test,&(mat.mu[i][j]));
            res = res || (0!=ibz_cmp(&test,&(mat.basis[1][0])));
            ibq_denom(&test,&(mat.r[i][j]));
            res = res || (0!=ibz_cmp(&test,&(mat.basis[1][0])));
        }
    }
    if(res){
        printf("    Test bkz_matrix_finit failed\n");
    }
    quat_bkz_matrix_finalize(&mat);
    ibz_finalize(&test);
    return(res);
}


//void ibq_round(ibz_t *rounded, ibq_t *q);
int quat_test_bkz_ibq_round(){
    int res = 0;
    ibz_t cmp, out;
    ibq_t x;
    ibz_init(&cmp);
    ibz_init(&out);
    ibq_init(&x);
    quat_test_helper_ibq_set_i(&x,1,3);
    ibz_set(&cmp,0);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,-1,3);
    ibz_set(&cmp,0);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,1,2);
    ibz_set(&cmp,1);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,9,5);
    ibz_set(&cmp,2);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,-9,5);
    ibz_set(&cmp,-2);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,-11,5);
    ibz_set(&cmp,-2);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,23,7);
    ibz_set(&cmp,3);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,21,7);
    ibz_set(&cmp,3);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,3,-1);
    ibz_set(&cmp,-3);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    quat_test_helper_ibq_set_i(&x,0,7);
    ibz_set(&cmp,0);
    ibq_round(&out,&x);
    res = res || (ibz_cmp(&out,&cmp)!=0);
    if(res){
        printf("    Test bkz_ibq_round failed\n");
    }
    ibz_finalize(&cmp);
    ibz_finalize(&out);
    ibq_finalize(&x);
    return(res);
}

//void ibz_vec_2x4_mul(ibz_vec_4_t *a, ibz_vec_4_t *b,const ibz_mat_2x2_t *U);
int quat_test_bkz_vec_2x4_mul(){
    int res = 0;
    ibz_vec_4_t a,b,c,d;
    ibz_mat_2x2_t U;
    ibz_mat_2x2_init(&U);
    ibz_vec_4_init(&a);
    ibz_vec_4_init(&b);
    ibz_vec_4_init(&c);
    ibz_vec_4_init(&d);
    ibz_vec_4_set(&a,1,2,-1,0);
    ibz_vec_4_set(&b,0,-3,5,2);
    quat_test_helper_ibz_mat_2x2_set(&U,1,0,0,1);
    ibz_vec_4_set(&c,1,2,-1,0);
    ibz_vec_4_set(&d,0,-3,5,2);
    ibz_vec_2x4_mul(&a,&b,&U);
    res = res || quat_test_helper_ibz_vec_4_equal(&a,&c);
    res = res || quat_test_helper_ibz_vec_4_equal(&b,&d);
    quat_test_helper_ibz_mat_2x2_set(&U,0,1,1,0);
    ibz_vec_2x4_mul(&a,&b,&U);
    res = res || quat_test_helper_ibz_vec_4_equal(&a,&d);
    res = res || quat_test_helper_ibz_vec_4_equal(&b,&c);
    ibz_vec_4_set(&a,1,2,-1,0);
    ibz_vec_4_set(&b,0,-3,5,2);
    quat_test_helper_ibz_mat_2x2_set(&U,1,-1,2,0);
    ibz_vec_4_set(&c,1,5,-6,-2);
    ibz_vec_4_set(&d,2,4,-2,0);
    ibz_vec_2x4_mul(&a,&b,&U);
    res = res || quat_test_helper_ibz_vec_4_equal(&a,&c);
    res = res || quat_test_helper_ibz_vec_4_equal(&b,&d);
    if(res){
        printf("    Test bkz_vec_2x4_mul failed\n");
    }
    ibz_mat_2x2_finalize(&U);
    ibz_vec_4_finalize(&a);
    ibz_vec_4_finalize(&b);
    ibz_vec_4_finalize(&c);
    ibz_vec_4_finalize(&d);
    return(res);
}

//void ibz_vec_4_bilinear(ibz_t *prod, const ibz_vec_4_t *a, const ibz_vec_4_t *b, const quat_alg_t *alg);
int quat_test_bkz_bilinear(){
    int res = 0;
    ibz_vec_4_t a, b;
    ibz_t prod, cmp;
    quat_alg_t alg;
    ibz_vec_4_init(&a);
    ibz_vec_4_init(&b);
    ibz_init(&prod);
    ibz_init(&cmp);
    ibz_set(&prod,11);
    quat_alg_init_set(&alg,&prod);
    ibz_vec_4_set(&a,1,2,-3,-1);
    ibz_vec_4_set(&b,-2,3,5,1);
    ibz_vec_4_bilinear(&prod,&a,&b,&alg);
    ibz_set(&cmp, -2*1+6-15*11-11);
    res = res || (0!=ibz_cmp(&cmp,&prod));
    ibz_vec_4_set(&a,0,1,1,0);
    ibz_vec_4_set(&b,1,2,3,1);
    ibz_vec_4_bilinear(&prod,&a,&b,&alg);
    ibz_set(&cmp, 2+33);
    res = res || (0!=ibz_cmp(&cmp,&prod));

    if(res){
        printf("    Test bkz_bilinear failed\n");
    }
    ibz_vec_4_finalize(&a);
    ibz_vec_4_finalize(&b);
    ibz_finalize(&prod);
    ibz_finalize(&cmp);
    quat_alg_finalize(&alg);
    return(res);
}


//void quat_bkz_matrix_set(quat_bkz_matrix_t *g, const ibz_mat_4x4_t *mat, quat_alg_t *alg);
int quat_test_bkz_set_matrix(){
    int res = 0;
    quat_bkz_matrix_t g,cmp;
    quat_alg_t alg;
    ibz_t p;
    ibz_mat_4x4_t mat;
    ibz_mat_4x4_init(&mat);
    quat_bkz_matrix_init(&g);
    quat_bkz_matrix_init(&cmp);
    ibz_init(&p);
    ibz_set(&p,1);
    quat_alg_init_set(&alg,&p);
    ibz_set(&(mat[0][0]),1);
    ibz_set(&(mat[1][1]),2);
    ibz_set(&(mat[2][2]),15);
    ibz_set(&(mat[3][3]),-2);
    ibz_set(&(mat[3][2]),3);
    ibz_set(&(mat[0][2]),-6);
    ibz_mat_4x4_transpose(&(cmp.basis),&mat);
    ibz_set(&(cmp.G[0][0]),1);
    ibz_set(&(cmp.G[1][1]),4);
    ibz_set(&(cmp.G[2][0]),-6);
    ibz_set(&(cmp.G[2][2]),270);
    ibz_set(&(cmp.G[3][2]),-6);
    ibz_set(&(cmp.G[3][3]),4);
    quat_test_helper_ibq_set_i(&(cmp.mu[2][0]),-6,1);
    quat_test_helper_ibq_set_i(&(cmp.mu[3][2]),-1,39);
    quat_test_helper_ibq_set_i(&(cmp.r[0][0]),1,1);
    quat_test_helper_ibq_set_i(&(cmp.r[1][1]),4,1);
    quat_test_helper_ibq_set_i(&(cmp.r[2][0]),-6,1);
    quat_test_helper_ibq_set_i(&(cmp.r[2][2]),234,1);
    quat_test_helper_ibq_set_i(&(cmp.r[3][2]),-6,1);
    quat_test_helper_ibq_set_i(&(cmp.r[3][3]),50,13);
    quat_bkz_matrix_set(&g,&mat,&alg);
    res = res || quat_bkz_matrix_equal(&cmp,&g);
    //compare to sage
    //G: [[1], [0, 4], [-6, 0, 270], [0, 0, -6, 4]]
    //Mu: [[], [0], [-6, 0], [0, 0, -1/39]]
    //r: [[1], [0, 4], [-6, 0, 234], [0, 0, -6, 50/13]]

    ibz_set(&(alg.p),11);
    ibz_set(&(alg.gram[2][2]),11);
    ibz_set(&(alg.gram[3][3]),11);
    ibz_set(&(cmp.G[2][2]),2610);
    ibz_set(&(cmp.G[3][2]),-66);
    ibz_set(&(cmp.G[3][3]),44);
    quat_test_helper_ibq_set_i(&(cmp.r[2][2]),2574,1);
    quat_test_helper_ibq_set_i(&(cmp.r[3][2]),-66,1);
    quat_test_helper_ibq_set_i(&(cmp.r[3][3]),550,13);
    quat_bkz_matrix_set(&g,&mat,&alg);
    res = res | quat_bkz_matrix_equal(&cmp,&g);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res = res || (0!=ibz_cmp(&(g.basis[i][j]),&(mat[j][i])));
            if(i<j){
                res = res || (0!=quat_test_helper_ibz_equal_i(&(g.G[i][j]),0));
                res = res || (0!=quat_test_helper_ibq_equal_i(&(g.r[i][j]),0,1));
                res = res || (0!=quat_test_helper_ibq_equal_i(&(g.mu[i][j]),0,1));
            }
        }
    }
    //[[1], [0, 4], [-6, 0, 2610], [0, 0, -66, 44]]
    //[[], [0], [-6, 0], [0, 0, -1/39]]
    //[[1], [0, 4], [-6, 0, 2574], [0, 0, -66, 550/13]]
    if(res){
        printf("    Test bkz_set_matrix failed\n");
    }
    ibz_finalize(&p);
    quat_alg_finalize(&alg);
    ibz_mat_4x4_finalize(&mat);
    quat_bkz_matrix_finalize(&g);
    quat_bkz_matrix_finalize(&cmp);
    return(res);
};

//void quat_bkz_update_matrix(quat_bkz_matrix_t *g, int start, int end, quat_alg_t *alg);
int quat_test_bkz_update_matrix(){
    int res = 0;
    ibz_t p;
    quat_alg_t alg;
    quat_bkz_matrix_t g,b;
    ibz_mat_4x4_t mat;
    int start = 1;
    int end = 3;
    quat_bkz_matrix_init(&g);
    quat_bkz_matrix_init(&b);
    ibz_mat_4x4_init(&mat);
    ibz_init(&p);
    ibz_set(&p,19);
    quat_alg_init_set(&alg,&p);
    ibz_set(&(mat[0][1]),1);
    ibz_set(&(mat[1][2]),7);
    ibz_set(&(mat[0][2]),-42);
    ibz_set(&(mat[2][3]),14);
    ibz_set(&(mat[1][3]),14);
    ibz_set(&(mat[0][3]),-4);
    ibz_set(&(mat[3][0]),24);
    quat_bkz_matrix_set(&g,&mat,&alg);
    ibz_set(&(g.basis[1][0]),2);
    ibz_set(&(g.basis[2][0]),42);
    ibz_set(&(g.basis[2][1]),-21);
    ibz_vec_4_bilinear(&(g.G[2][1]),&(g.basis[1]),&(g.basis[2]),&alg);
    ibz_vec_4_bilinear(&(g.G[3][1]),&(g.basis[1]),&(g.basis[3]),&alg);
    ibz_vec_4_bilinear(&(g.G[3][2]),&(g.basis[2]),&(g.basis[3]),&alg);
    ibz_vec_4_bilinear(&(g.G[1][1]),&(g.basis[1]),&(g.basis[1]),&alg);
    ibz_vec_4_bilinear(&(g.G[2][2]),&(g.basis[2]),&(g.basis[2]),&alg);
    quat_bkz_update_matrix(&g,start,end);
    ibz_set(&(mat[0][1]),2);
    ibz_set(&(mat[0][2]),42);
    ibz_set(&(mat[1][2]),-21);
    quat_bkz_matrix_set(&b,&mat,&alg);
    for(int i =0; i < end; i++){
        for(int j =0; j < i; j++){
            res = res | (0!=ibz_cmp(&(g.basis[i][j]),&(b.basis[i][j])));
            res = res | (0!=ibz_cmp(&(g.G[i][j]),&(b.G[i][j])));
            res = res | (0!=ibq_cmp(&(g.mu[i][j]),&(b.mu[i][j])));
            res = res | (0!=ibq_cmp(&(g.r[i][j]),&(b.r[i][j])));
        }
    }
    if(res){
        printf("    Test bkz_update_matrix failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    quat_bkz_matrix_finalize(&g);
    quat_bkz_matrix_finalize(&b);
    ibz_finalize(&p);
    quat_alg_finalize(&alg);
    return(res);
};

//void quat_bkz_size_reduce(quat_bkz_matrix_t *g, int start, int end, quat_alg_t *alg);
int quat_test_bkz_size_reduce(){
    int res = 0;
    ibz_t p;
    ibq_t tmp, half;
    quat_alg_t alg;
    quat_bkz_matrix_t g,b;
    ibz_mat_4x4_t mat;
    int start = 1;
    int end = 3;
    quat_bkz_matrix_init(&g);
    quat_bkz_matrix_init(&b);
    ibz_mat_4x4_init(&mat);
    ibz_init(&p);
    ibq_init(&tmp);
    ibq_init(&half);
    ibz_set(&p,19);
    quat_alg_init_set(&alg,&p);
    ibz_set(&(mat[0][1]),1);
    ibz_set(&(mat[1][2]),7);
    ibz_set(&(mat[0][2]),-42);
    ibz_set(&(mat[2][3]),14);
    ibz_set(&(mat[1][3]),14);
    ibz_set(&(mat[0][3]),-4);
    ibz_set(&(mat[3][0]),24);
    ibz_set(&(mat[3][1]),42);
    quat_bkz_matrix_set(&g,&mat,&alg);
    quat_bkz_size_reduce(&g,0,2);
    res = res | quat_test_helper_bkz_test_integrity(&g,&alg);
    quat_test_helper_ibq_set_i(&half,1,2);
    for(int i = 0; i<2;i++){
        for(int j = 0; j<i;j++){
            ibq_abs(&tmp,&(g.mu[i][j]));
            res = res | (0<ibq_cmp(&tmp,&half));
        }
    }
    quat_bkz_matrix_set(&g,&mat,&alg);
    quat_bkz_size_reduce(&g,0,4);
    res = res | quat_test_helper_bkz_test_integrity(&g,&alg);
    quat_test_helper_ibq_set_i(&half,1,2);
    for(int i = 0; i<4;i++){
        for(int j = 0; j<i;j++){
            ibq_abs(&tmp,&(g.mu[i][j]));
            res = res | (0<ibq_cmp(&tmp,&half));
        }
    }
    if(res){
        printf("    Test bkz_size_reduce failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    quat_bkz_matrix_finalize(&g);
    quat_bkz_matrix_finalize(&b);
    ibz_finalize(&p);
    quat_alg_finalize(&alg);
    ibq_finalize(&tmp);
    ibq_finalize(&half);
    return(res);
};

//void quat_bkz_update_after_lagrange(quat_bkz_matrix_t *g, const ibz_mat_2x2_t *U, int index);
int quat_test_bkz_update_after_lagrange(){
    int res = 0;
    ibz_mat_2x2_t U;
    ibz_mat_4x4_t mat;
    quat_bkz_matrix_t g;
    quat_alg_t alg;
    ibz_t p;
    ibz_init(&p);
    ibz_mat_4x4_init(&mat);
    ibz_mat_2x2_init(&U);
    quat_bkz_matrix_init(&g);
    ibz_set(&p,19);
    quat_alg_init_set(&alg,&p);

    ibz_vec_4_set(&(mat[0]),4,0,2,0);
    ibz_vec_4_set(&(mat[1]),2,0,4,30);
    ibz_vec_4_set(&(mat[2]),0,0,8,24);
    ibz_vec_4_set(&(mat[3]),10,3,0,72);
    quat_bkz_matrix_set(&g,&mat,&alg);
    ibz_set(&(U[0][0]),1);
    ibz_set(&(U[0][1]),1);
    ibz_set(&(U[1][1]),1);
    ibz_vec_2x4_mul(&(g.basis[1]),&(g.basis[2]),&U);
    quat_bkz_update_after_lagrange(&g,&U,1);
    quat_bkz_update_matrix(&g,1,4);
    res = res | quat_test_helper_bkz_test_integrity(&g,&alg);

    if(res){
        printf("    Test bkz_update_after_lagrange failed\n");
    }
    ibz_finalize(&p);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_2x2_init(&U);
    quat_bkz_matrix_finalize(&g);
    quat_alg_finalize(&alg);
    return(res);
};

//void quat_bkz_lagrange_reduction_gram(ibz_mat_2x2_t *U, ibz_mat_2x2_t *G, int lagrange_tours);
int quat_test_bkz_lagrange_reduction_gram(){
    int res = 0;
    ibz_mat_2x2_t U;
    ibq_mat_2x2_t G;
    ibz_vec_4_t a,b;
    quat_alg_t alg;
    ibz_t p;
    ibz_t one;
    ibq_t cmp;
    ibq_init(&cmp);
    ibz_init(&p);
    ibz_init(&one);
    ibq_mat_2x2_init(&G);
    ibz_mat_2x2_init(&U);
    ibz_vec_4_init(&a);
    ibz_vec_4_init(&b);
    ibz_set(&p,19);
    quat_alg_init_set(&alg,&p);
    ibz_set(&one,1);
    ibz_vec_4_set(&a,19,23,0,-3);
    ibz_vec_4_set(&b,-70,2,1,0);
    ibz_vec_4_bilinear(&p,&a,&a,&alg);
    ibq_set(&(G[0][0]),&p,&one);
    ibz_vec_4_bilinear(&p,&a,&b,&alg);
    ibq_set(&(G[1][0]),&p,&one);
    ibz_vec_4_bilinear(&p,&b,&b,&alg);
    ibq_set(&(G[1][1]),&p,&one);
    quat_bkz_lagrange_reduction_gram(&U,&G,20);
    ibz_vec_2x4_mul(&a,&b,&U);
    //test if consistent
    ibz_vec_4_bilinear(&p,&a,&a, &alg);
    ibq_set(&cmp,&p,&one);
    res = res | (0!=ibq_cmp(&cmp,&(G[0][0])));
    ibz_vec_4_bilinear(&p,&a,&b, &alg);
    ibq_set(&cmp,&p,&one);
    res = res | (0!=ibq_cmp(&cmp,&(G[1][0])));
    ibz_vec_4_bilinear(&p,&b,&b, &alg);
    ibq_set(&cmp,&p,&one);
    res = res | (0!=ibq_cmp(&cmp,&(G[1][1])));
    //test if reduced by comparing to output of sage version
    // [1061, 0], [-223, 3416]
    res = res | (0<quat_test_helper_ibq_equal_i(&(G[0][0]),1061,1));
    res = res | (0<quat_test_helper_ibq_equal_i(&(G[1][1]),3416,1));
    res = res | (0<ibq_cmp(&(G[0][0]),&(G[1][1])));

    if(res){
        printf("    Test bkz_lagrange_reduction_gram failed\n");
    }
    ibz_finalize(&p);
    ibz_finalize(&one);
    ibq_finalize(&cmp);
    ibq_mat_2x2_finalize(&G);
    ibz_mat_2x2_finalize(&U);
    quat_alg_finalize(&alg);
    ibz_vec_4_finalize(&a);
    ibz_vec_4_finalize(&b);
    return(res);
};

//void quat_lattice_bkz(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, int bkz_tours, int lagrange_tours, quat_alg_t *alg);
int quat_test_bkz_lattice_bkz(){
    int res = 0;
    ibz_t p;
    quat_alg_t alg;
    quat_lattice_t lat;
    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);
    ibz_mat_4x4_init(&(lat.basis));
    ibz_init(&(lat.denom));
    ibz_init(&p);
    ibz_set(&p,19);
    quat_alg_init_set(&alg,&p);
//ideal generated by x=1+3i+j+k and 8
//[1 0 6 7]
//[0 1 1 6]
//[0 0 8 0]
//[0 0 0 8]
    ibz_vec_4_set(&(lat.basis[0]),1,0,0,0);
    ibz_vec_4_set(&(lat.basis[1]),0,1,0,0);
    ibz_vec_4_set(&(lat.basis[2]),6,1,8,0);
    ibz_vec_4_set(&(lat.basis[3]),7,6,0,8);
    ibz_set(&(lat.denom),1);
    quat_lattice_bkz(&red,&lat,100,30,&alg);
    ibz_set(&p,8);
    res = res | quat_test_helper_is_reduced(&red,&p,&alg);
//[ 1/2  3/2  5/2 27/2]
//[   0    2    2   12]
//[   0    0    4    4]
//[   0    0    0   16]
//norm 16
    ibz_vec_4_set(&(lat.basis[0]),1,0,0,0);
    ibz_vec_4_set(&(lat.basis[1]),3,4,0,0);
    ibz_vec_4_set(&(lat.basis[2]),5,4,8,0);
    ibz_vec_4_set(&(lat.basis[3]),27,24,8,32);
    ibz_set(&(lat.denom),2);
    quat_lattice_bkz(&red,&lat,100,30,&alg);
    ibz_set(&p,16*2*2);//multiply norm by 2^2 for denom
    res = res | quat_test_helper_is_reduced(&red,&p,&alg);
    
    if(res){
        printf("    Test bkz_lattice_bkz failed\n");
    }
    ibz_mat_4x4_finalize(&red);
    ibz_mat_4x4_finalize(&(lat.basis));
    ibz_finalize(&(lat.denom));
    ibz_finalize(&p);
    quat_alg_finalize(&alg);
    return(res);
};

int bkz_tests(){
    int res = 0;
    printf("Run bkz tests:\n");
    res = res | quat_test_bkz_bilinear();
    res = res | quat_test_bkz_matrix_finit();
    res = res | quat_test_bkz_ibq_round();
    res = res | quat_test_bkz_vec_2x4_mul();
    res = res | quat_test_bkz_set_matrix();
    res = res | quat_test_bkz_update_matrix();
    res = res | quat_test_bkz_size_reduce();
    res = res | quat_test_bkz_update_after_lagrange();
    res = res | quat_test_bkz_lagrange_reduction_gram();
    res = res | quat_test_bkz_lattice_bkz();
    if(!res){
        printf("All passed\n");
    } else{
        printf("Some tests failed\n");
    }
    return res;
}
