
#include "nct_helpers.h"
#include <stdio.h>

// functions whose equivalent is in ct_intbig

void nct_ibz_mat_4x4_transpose(nct_ibz_mat_4x4_t *transposed, const nct_ibz_mat_4x4_t *mat){
    nct_ibz_mat_4x4_t work;
    nct_ibz_mat_4x4_init(&work);
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            nct_ibz_copy(&(work[i][j]),&((*mat)[j][i]));
        }
    }
    nct_ibz_mat_4x4_copy(transposed,&work);
    nct_ibz_mat_4x4_finalize(&work);
}

void nct_ibz_mat_4x4_copy(nct_ibz_mat_4x4_t *copy, const nct_ibz_mat_4x4_t *copied){
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            nct_ibz_copy(&((*copy)[i][j]),&((*copied)[i][j]));
        }
    }
}

void nct_ibz_mat_4x4_init(nct_ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_init(&(*mat)[i][j]);
        }
    }
}
void nct_ibz_mat_2x2_init(nct_ibz_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            nct_ibz_init(&(*x)[i][j]);
        }
    }
}
void nct_ibq_mat_2x2_init(nct_ibq_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            nct_ibq_init(&(*x)[i][j]);
        }
    }
}
void nct_ibq_mat_4x4_init(nct_ibq_mat_4x4_t *x){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibq_init(&(*x)[i][j]);
        }
    }
}
void nct_ibz_vec_4_init(nct_ibz_vec_4_t *x){
    for(int j = 0; j < 4; j++){
        nct_ibz_init(&(*x)[j]);
    }
}

void nct_ibz_mat_4x4_finalize(nct_ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_finalize(&(*mat)[i][j]);
        }
    }
}
void nct_ibz_mat_2x2_finalize(nct_ibz_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            nct_ibz_finalize(&(*x)[i][j]);
        }
    }
}
void nct_ibq_mat_2x2_finalize(nct_ibq_mat_2x2_t *x){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            nct_ibq_finalize(&(*x)[i][j]);
        }
    }
}
void nct_ibq_mat_4x4_finalize(nct_ibq_mat_4x4_t *x){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibq_finalize(&(*x)[i][j]);
        }
    }
}
void nct_ibz_vec_4_finalize(nct_ibz_vec_4_t *x){
    for(int j = 0; j < 4; j++){
        nct_ibz_finalize(&(*x)[j]);
    }
}

// functions whose equivalent is in ct_helpers

void nct_quat_alg_init_set(nct_quat_alg_t *alg, const nct_ibz_t *p){
    nct_ibz_init(&(*alg).p);
    nct_ibz_mat_4x4_init(&(*alg).gram);
    nct_ibz_copy(&(*alg).p, p);
    nct_ibz_set(&(*alg).gram[0][0], 1);
    nct_ibz_set(&(*alg).gram[1][1], 1);
    nct_ibz_copy(&(*alg).gram[2][2], p);
    nct_ibz_copy(&(*alg).gram[3][3], p);
}

void nct_quat_alg_init_set_ui(nct_quat_alg_t *alg, unsigned int p) {
    nct_ibz_t bp;
    nct_ibz_init(&bp);
    nct_ibz_set(&bp, p);
    nct_quat_alg_init_set(alg, &bp);
    nct_ibz_finalize(&bp);
}

void nct_quat_alg_finalize(nct_quat_alg_t *alg){
    nct_ibz_finalize(&(*alg).p);
    nct_ibz_mat_4x4_finalize(&(*alg).gram);
}

void nct_quat_alg_elem_init(nct_quat_alg_elem_t *elem){
    nct_ibz_vec_4_init(&(*elem).coord);
    nct_ibz_init(&(*elem).denom);
    nct_ibz_set(&(*elem).denom, 1);
}
void nct_quat_alg_elem_finalize(nct_quat_alg_elem_t *elem){
    nct_ibz_vec_4_finalize(&(*elem).coord);
    nct_ibz_finalize(&(*elem).denom);
}

void nct_quat_lattice_init(nct_quat_lattice_t *lat){
    nct_ibz_mat_4x4_init(&(*lat).basis);
    nct_ibz_init(&(*lat).denom);
    nct_ibz_set(&(*lat).denom, 1);
}
void nct_quat_lattice_finalize(nct_quat_lattice_t *lat){
    nct_ibz_finalize(&(*lat).denom);
    nct_ibz_mat_4x4_finalize(&(*lat).basis);
}

// set integer vector coefficients to given integers
void nct_ibz_vec_4_set(nct_ibz_vec_4_t *vec, int a, int b, int c, int d){
    nct_ibz_set(&((*vec)[0]),a);
    nct_ibz_set(&((*vec)[1]),b);
    nct_ibz_set(&((*vec)[2]),c);
    nct_ibz_set(&((*vec)[3]),d);
}

void nct_ibz_vec_4_bilinear(nct_ibz_t *prod, const nct_ibz_vec_4_t *a, const nct_ibz_vec_4_t *b, const nct_quat_alg_t *alg){
    nct_ibz_t tmp;
    nct_ibz_init(&tmp);
    nct_ibz_mul(prod,&((*a)[0]),&((*b)[0]));
    nct_ibz_mul(&tmp,&((*a)[1]),&((*b)[1]));
    nct_ibz_add(prod,prod,&tmp);
    nct_ibz_mul(&tmp,&((*a)[2]),&((*b)[2]));
    nct_ibz_mul(&tmp,&tmp,&(alg->p));
    nct_ibz_add(prod,prod,&tmp);
    nct_ibz_mul(&tmp,&((*a)[3]),&((*b)[3]));
    nct_ibz_mul(&tmp,&tmp,&(alg->p));
    nct_ibz_add(prod,prod,&tmp);
    nct_ibz_finalize(&tmp);
}

// set rational to fraction of two small integers
void nct_quat_test_helper_nct_ibq_set_i(nct_ibq_t *q, int n, int d){
    nct_ibz_t num,denom;
    nct_ibz_init(&num);
    nct_ibz_init(&denom);
    nct_ibz_set(&num,n);
    nct_ibz_set(&denom,d);
    nct_ibq_set(q,&num,&denom);
    nct_ibz_finalize(&num);
    nct_ibz_finalize(&denom);
}

// test equality of integer vectors, return 0 if equal, 1 if not
int nct_quat_test_helper_nct_ibz_vec_4_equal(const nct_ibz_vec_4_t *vec, const nct_ibz_vec_4_t *cmp){
    int res = 0;
    res = res | (nct_ibz_cmp(&((*vec)[0]),&((*cmp)[0]))!=0);
    res = res | (nct_ibz_cmp(&((*vec)[1]),&((*cmp)[1]))!=0);
    res = res | (nct_ibz_cmp(&((*vec)[2]),&((*cmp)[2]))!=0);
    res = res | (nct_ibz_cmp(&((*vec)[3]),&((*cmp)[3]))!=0);
    return(res);
}

// compares integer to small int via nct_ibz_cmp
int nct_quat_test_helper_nct_ibz_equal_i(const nct_ibz_t *x, int cmp){
    nct_ibz_t c;
    nct_ibz_init(&c);
    nct_ibz_set(&c,cmp);
    int res = nct_ibz_cmp(x,&c);
    nct_ibz_finalize(&c);
    return(res);
}

//return 0 if equal, 1 otherwise
int nct_ibz_mat_4x4_equal(const nct_ibz_mat_4x4_t *a, const nct_ibz_mat_4x4_t *b){
    int res = 0;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res = res | nct_ibz_cmp(&((*a)[i][j]),&((*b)[i][j]));
        }
    }
    return(res!=0);
}

//return 0 if reduced enough, 1 otherwise
int nct_quat_test_helper_is_reduced(const nct_ibz_mat_4x4_t *mat, const nct_ibz_t *norm_ideal, const nct_quat_alg_t *alg, int print_flag){
    nct_ibz_mat_4x4_t trans;
    nct_ibz_t prod, norm, ref;
    nct_ibz_init(&prod);
    nct_ibz_init(&norm);
    nct_ibz_init(&ref);
    nct_ibz_mat_4x4_init(&trans);
    nct_ibz_mat_4x4_transpose(&trans,mat);
    nct_ibz_set(&prod,1);
    for (int i = 0; i < 4; i++){
        nct_ibz_vec_4_bilinear(&(norm),&(trans[i]),&(trans[i]),alg);
        nct_ibz_mul(&prod,&prod,&norm);
    }
    nct_ibz_mul(&ref,&(alg->p),&(alg->p));
    nct_ibz_set(&norm,4);
    nct_ibz_mul(&ref,&ref,&norm);
    nct_ibz_mul(&norm,norm_ideal,norm_ideal);
    nct_ibz_mul(&norm,&norm,&norm);
    nct_ibz_mul(&ref,&norm,&ref);
    int res = (0<nct_ibz_cmp(&prod,&ref));
    if (print_flag){
        printf("Minkowski\n");
        nct_ibz_printf("rhs = %Zd\n", &prod);
        nct_ibz_printf("lsh = %Zd\n", &ref);
    }
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&norm);
    nct_ibz_finalize(&ref);
    nct_ibz_mat_4x4_finalize(&trans);
    return(res);
}

void nct_ibz_content(nct_ibz_t *content, const nct_ibz_vec_4_t *v) {
  nct_ibz_gcd(content, &((*v)[0]), &((*v)[1]));
  nct_ibz_gcd(content, &((*v)[2]), content);
  nct_ibz_gcd(content, &((*v)[3]), content);
}

void nct_ibz_vec_4_set_ibz(nct_ibz_vec_4_t *vec, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3){
    nct_ibz_set(&((*vec)[0]),coord0);
    nct_ibz_set(&((*vec)[1]),coord1);
    nct_ibz_set(&((*vec)[2]),coord2);
    nct_ibz_set(&((*vec)[3]),coord3);
}

void nct_ibz_vec_4_copy(nct_ibz_vec_4_t *new, const nct_ibz_vec_4_t  *vec){
    for (int i = 0; i <4; i++){
        nct_ibz_copy(&((*new)[i]),&((*vec)[i]));
    }
}

void nct_ibz_vec_4_copy_ibz(nct_ibz_vec_4_t *res, const nct_ibz_t *coord0,const nct_ibz_t *coord1,const nct_ibz_t *coord2,const nct_ibz_t *coord3){
  nct_ibz_copy(&((*res)[0]),coord0);
  nct_ibz_copy(&((*res)[1]),coord1);
  nct_ibz_copy(&((*res)[2]),coord2);
  nct_ibz_copy(&((*res)[3]),coord3);
}

void nct_ibz_vec_4_negate(nct_ibz_vec_4_t *neg, const nct_ibz_vec_4_t  *vec){
    for (int i = 0; i <4; i++){
        nct_ibz_neg(&((*neg)[i]),&((*vec)[i]));
    }
}

void nct_ibz_vec_4_add(nct_ibz_vec_4_t *res, const nct_ibz_vec_4_t *a, const nct_ibz_vec_4_t *b){
  nct_ibz_add(&((*res)[0]),&((*a)[0]),&((*b)[0]));
  nct_ibz_add(&((*res)[1]),&((*a)[1]),&((*b)[1]));
  nct_ibz_add(&((*res)[2]),&((*a)[2]),&((*b)[2]));
  nct_ibz_add(&((*res)[3]),&((*a)[3]),&((*b)[3]));
}

void nct_ibz_vec_4_sub(nct_ibz_vec_4_t *res, const nct_ibz_vec_4_t *a, const nct_ibz_vec_4_t *b){
  nct_ibz_sub(&((*res)[0]),&((*a)[0]),&((*b)[0]));
  nct_ibz_sub(&((*res)[1]),&((*a)[1]),&((*b)[1]));
  nct_ibz_sub(&((*res)[2]),&((*a)[2]),&((*b)[2]));
  nct_ibz_sub(&((*res)[3]),&((*a)[3]),&((*b)[3]));
}

int nct_ibz_vec_4_is_zero(const nct_ibz_vec_4_t *x){
  int res = 1;
  for (int i = 0; i < 4; i++){
    res &= nct_ibz_is_zero(&((*x)[i]));
  }
  return(res);
}

void nct_ibz_mat_4x4_zero(nct_ibz_mat_4x4_t *zero){
    for(int i = 0; i <4; i++){
        for(int j = 0;  j<4; j++){
            nct_ibz_set(&((*zero)[i][j]),0);
        }
    }
}

void nct_ibz_mat_4x4_identity(nct_ibz_mat_4x4_t *id){
    for(int i = 0; i <4; i++){
        for(int j = 0;  j<4; j++){
            nct_ibz_set(&((*id)[i][j]),0);
        }
        nct_ibz_set(&((*id)[i][i]),1);
    }
}

int nct_ibz_mat_4x4_is_identity(const nct_ibz_mat_4x4_t *mat){
    int res = 1;
    for(int i = 0; i <4; i++){
        for(int j = 0;  j<4; j++){
            res = res && (0==nct_quat_test_helper_nct_ibz_equal_i(&((*mat)[i][j]),(i==j)));
        }
    }
    return(res);
}

void nct_ibz_vec_4_scalar_mul(nct_ibz_vec_4_t *prod, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec){
    for(int i = 0; i < 4; i++){
        nct_ibz_mul(&((*prod)[i]),&((*vec)[i]),scalar);
    }
}

void nct_ibz_mat_4x4_scalar_mul(nct_ibz_mat_4x4_t *prod, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_mul(&((*prod)[i][j]),&((*mat)[i][j]),scalar);
        }
    }
}

void nct_ibz_mat_4x4_eval(nct_ibz_vec_4_t  *res, const nct_ibz_mat_4x4_t *mat, const nct_ibz_vec_4_t *vec){
    nct_ibz_vec_4_t sum;
    nct_ibz_t prod;
    nct_ibz_init(&prod);
    nct_ibz_vec_4_init(&sum);
    //suppose initialization to 0
    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            nct_ibz_mul(&prod, &(*mat)[i][j], &(*vec)[j]);
            nct_ibz_add(&(sum[i]),&(sum[i]), &prod);
        }
    }
    nct_ibz_vec_4_copy(res,&sum);
    nct_ibz_finalize(&prod);
    nct_ibz_vec_4_finalize(&sum);
}

void nct_ibz_mat_4x4_gcd(nct_ibz_t *gcd, const nct_ibz_mat_4x4_t *mat){
    nct_ibz_t d;
    nct_ibz_init(&d);
    nct_ibz_copy(&d, &((*mat)[0][0]));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_gcd(&d,&d,&((*mat)[i][j]));
        }
    }
    nct_ibz_copy(gcd,&d);
    nct_ibz_finalize(&d);
}

// Triangular matrix operations
int nct_ibz_mat_4x4_triangular_equal(const nct_ibz_mat_4x4_t *mat1, const nct_ibz_mat_4x4_t *mat2){
    int res = 1;
    for(int i = 0; i < 4; i++){
        for(int j = i; j < 4; j++){
            res = res && !nct_ibz_cmp(&((*mat1)[i][j]),&((*mat2)[i][j]));
        }
    }
    return(res);
}

void nct_ibz_mat_4x4_triangular_scalar_mul(nct_ibz_mat_4x4_t *prod, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat){
    nct_ibz_t s;
    nct_ibz_init(&s);
    nct_ibz_copy(&s,scalar);
    for(int i = 0; i < 4; i++){
        for(int j = i; j < 4; j++){
            nct_ibz_mul(&((*prod)[i][j]),&s,&((*mat)[i][j]));
        }
    }
    nct_ibz_finalize(&s);
}

void nct_ibz_mat_4x4_triangular_eval(nct_ibz_vec_4_t  *res, const nct_ibz_mat_4x4_t *mat, const nct_ibz_vec_4_t *vec){
    nct_ibz_vec_4_t sum;
    nct_ibz_t prod;
    nct_ibz_init(&prod);
    nct_ibz_vec_4_init(&sum);
    for (int i = 0; i <4; i++){
        for (int j = i; j <4; j++){
            nct_ibz_mul(&prod, &(*mat)[i][j], &(*vec)[j]);
            nct_ibz_add(&(sum[i]),&(sum[i]), &prod);
        }
    }
    nct_ibz_vec_4_copy(res,&sum);
    nct_ibz_finalize(&prod);
    nct_ibz_vec_4_finalize(&sum);
}

//internal helper functions
void nct_ibz_mat_4x4_mul(nct_ibz_mat_4x4_t *res, const nct_ibz_mat_4x4_t *a, const nct_ibz_mat_4x4_t *b){
    nct_ibz_mat_4x4_t mat;
    nct_ibz_t prod;
    nct_ibz_init(&prod);
    nct_ibz_mat_4x4_init(&mat);
    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            nct_ibz_set(&(mat[i][j]),0);
            for (int k = 0; k <4; k++){
                nct_ibz_mul(&prod,&((*a)[i][k]), &((*b)[k][j]));
                nct_ibz_add(&(mat[i][j]), &(mat[i][j]), &prod);
            }
        }
    }
    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            nct_ibz_copy(&((*res)[i][j]),&(mat[i][j]));
        }
    }
    nct_ibz_mat_4x4_finalize(&mat);
    nct_ibz_finalize(&prod);
}

//4x4 inversion helper functions
void nct_ibz_mat_2x2_det_from_ibz(nct_ibz_t *det, const nct_ibz_t *a11, const nct_ibz_t *a12, const nct_ibz_t *a21, const nct_ibz_t *a22){
    nct_ibz_t prod;
    nct_ibz_init(&prod);
    nct_ibz_mul(&prod,a12,a21);
    nct_ibz_mul(det,a11,a22);
    nct_ibz_sub(det,det,&prod);
    nct_ibz_finalize(&prod);
}

void nct_ibz_inv_dim4_make_coeff_pmp(nct_ibz_t *coeff, const nct_ibz_t *a1, const nct_ibz_t *a2, const nct_ibz_t *b1, const nct_ibz_t *b2, const nct_ibz_t *c1, const nct_ibz_t *c2){
    nct_ibz_t prod, sum;
    nct_ibz_init(&prod);
    nct_ibz_init(&sum);
    nct_ibz_mul(&sum,a1,a2);
    nct_ibz_mul(&prod,b1,b2);
    nct_ibz_sub(&sum,&sum,&prod);
    nct_ibz_mul(&prod,c1,c2);
    nct_ibz_add(coeff,&sum,&prod);
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&sum);
}

void nct_ibz_inv_dim4_make_coeff_mpm(nct_ibz_t *coeff, const nct_ibz_t *a1, const nct_ibz_t *a2, const nct_ibz_t *b1, const nct_ibz_t *b2, const nct_ibz_t *c1, const nct_ibz_t *c2){
    nct_ibz_t prod, sum;
    nct_ibz_init(&prod);
    nct_ibz_init(&sum);
    nct_ibz_mul(&sum,b1,b2);
    nct_ibz_mul(&prod,a1,a2);
    nct_ibz_sub(&sum,&sum,&prod);
    nct_ibz_mul(&prod,c1,c2);
    nct_ibz_sub(coeff,&sum,&prod);
    nct_ibz_finalize(&prod);
    nct_ibz_finalize(&sum);
}

//Method from https://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf 3rd of May 2023, 16h15 CEST
int nct_ibz_mat_4x4_inv_with_det_as_denom(nct_ibz_mat_4x4_t *inv, nct_ibz_t *det, const nct_ibz_mat_4x4_t *mat){
    int invertible;
    nct_ibz_t prod,work_det;
    nct_ibz_mat_4x4_t work;
    nct_ibz_t s[6];
    nct_ibz_t c[6];
    for (int i = 0; i < 6; i++){
        nct_ibz_init(&(s[i]));
        nct_ibz_init(&(c[i]));
    }
    nct_ibz_mat_4x4_init(&work);
    nct_ibz_init(&prod);
    nct_ibz_init(&work_det);

    //compute some 2x2 minors, store them in s and c
    for (int i = 0; i < 3; i++){
        nct_ibz_mat_2x2_det_from_ibz(&(s[i]),&((*mat)[0][0]),&((*mat)[0][i+1]),&((*mat)[1][0]),&((*mat)[1][i+1]));
        nct_ibz_mat_2x2_det_from_ibz(&(c[i]),&((*mat)[2][0]),&((*mat)[2][i+1]),&((*mat)[3][0]),&((*mat)[3][i+1]));
    }
    for (int i = 0; i < 2; i++){
        nct_ibz_mat_2x2_det_from_ibz(&(s[3+i]),&((*mat)[0][1]),&((*mat)[0][2+i]),&((*mat)[1][1]),&((*mat)[1][2+i]));
        nct_ibz_mat_2x2_det_from_ibz(&(c[3+i]),&((*mat)[2][1]),&((*mat)[2][2+i]),&((*mat)[3][1]),&((*mat)[3][2+i]));
    }
    nct_ibz_mat_2x2_det_from_ibz(&(s[5]),&((*mat)[0][2]),&((*mat)[0][3]),&((*mat)[1][2]),&((*mat)[1][3]));
    nct_ibz_mat_2x2_det_from_ibz(&(c[5]),&((*mat)[2][2]),&((*mat)[2][3]),&((*mat)[3][2]),&((*mat)[3][3]));

    //compute det
    nct_ibz_set(&work_det,0);
    for (int i = 0; i < 6; i++){
        nct_ibz_mul(&prod,&(s[i]),&(c[5-i]));
        if ((i != 1) && (i != 4)){
            nct_ibz_add(&work_det,&work_det,&prod);
        } else {
            nct_ibz_sub(&work_det,&work_det,&prod);
        }
    }
    //compute transposed adjugate
    for (int j = 0; j < 4; j++){
        for (int k = 0; k < 2; k++){
            if ((k + j + 1) % 2 == 1){
                nct_ibz_inv_dim4_make_coeff_pmp(&(work[j][k]), &((*mat)[1-k][(j==0)]), &(c[6-j-(j==0)]), &((*mat)[1-k][2-(j>1)]), &(c[4-j-(j==1)]), &((*mat)[1-k][3-(j==3)]), &(c[3-j-(j==1)-(j==2)]));
            } else {
               nct_ibz_inv_dim4_make_coeff_mpm(&(work[j][k]), &((*mat)[1-k][(j==0)]), &(c[6-j-(j==0)]), &((*mat)[1-k][2-(j>1)]), &(c[4-j-(j==1)]), &((*mat)[1-k][3-(j==3)]), &(c[3-j-(j==1)-(j==2)]));
            }
        }
        for (int k = 2; k < 4; k++){
            if ((k + j + 1) % 2 == 1){
                nct_ibz_inv_dim4_make_coeff_pmp(&(work[j][k]), &((*mat)[3-(k==3)][(j==0)]), &(s[6-j-(j==0)]),&((*mat)[3-(k==3)][2-(j>1)]), &(s[4-j-(j==1)]),&((*mat)[3-(k==3)][3-(j==3)]), &(s[3-j-(j==1)-(j==2)]));
            } else {
               nct_ibz_inv_dim4_make_coeff_mpm(&(work[j][k]), &((*mat)[3-(k==3)][(j==0)]), &(s[6-j-(j==0)]),&((*mat)[3-(k==3)][2-(j>1)]), &(s[4-j-(j==1)]),&((*mat)[3-(k==3)][3-(j==3)]), &(s[3-j-(j==1)-(j==2)]));
            }
        }
    }
    if(inv != NULL){
        // put transposed adjugate in result, or 0 if no inverse
        nct_ibz_set(&prod,!nct_ibz_is_zero(&work_det));
        nct_ibz_mat_4x4_scalar_mul(inv,&prod,&work);
    }
    //output det
    invertible = !nct_ibz_is_zero(&work_det);
    if(det != NULL)
        nct_ibz_copy(det,&work_det);
    for (int i = 0; i < 6; i++){
        nct_ibz_finalize(&s[i]);
        nct_ibz_finalize(&c[i]);
    }
    nct_ibz_mat_4x4_finalize(&work);
    nct_ibz_finalize(&work_det);
    nct_ibz_finalize(&prod);
    return(invertible);
}
