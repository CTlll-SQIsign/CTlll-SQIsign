#include "bkz.h"
#include <assert.h>

void quat_bkz_matrix_init(quat_bkz_matrix_t *g){
    ibz_mat_4x4_init(&(g->basis));
    ibz_mat_4x4_init(&(g->G));
    ibq_mat_4x4_init(&(g->mu));
    ibq_mat_4x4_init(&(g->r));
}

void quat_bkz_matrix_finalize(quat_bkz_matrix_t *g){
    ibz_mat_4x4_finalize(&(g->basis));
    ibz_mat_4x4_finalize(&(g->G));
    ibq_mat_4x4_finalize(&(g->mu));
    ibq_mat_4x4_finalize(&(g->r));
}

void ibz_vec_2x4_mul(ibz_vec_4_t *a, ibz_vec_4_t *b,const ibz_mat_2x2_t *U){
    ibz_t tmp, prod;
    ibz_init(&tmp);
    ibz_init(&prod);
    for(int i = 0; i < 4; i++){
        ibz_mul(&tmp,&((*a)[i]),&((*U)[0][0]));
        ibz_mul(&prod,&((*b)[i]),&((*U)[0][1]));
        ibz_add(&tmp,&tmp,&prod);
        ibz_mul(&prod,&((*a)[i]),&((*U)[1][0]));
        ibz_mul(&((*b)[i]),&((*b)[i]),&((*U)[1][1]));
        ibz_add(&((*b)[i]),&prod,&((*b)[i]));
        ibz_copy(&((*a)[i]),&tmp);
    }
    ibz_finalize(&tmp);
    ibz_finalize(&prod);
}
#include <stdio.h>
void quat_bkz_matrix_set(quat_bkz_matrix_t *g, const ibz_mat_4x4_t *mat, const quat_alg_t *alg){
    ibz_t den, num;
    ibq_t tmp;
    quat_lattice_t lat;
    ibq_init(&tmp);
    ibz_init(&den);
    ibz_init(&num);
    quat_lattice_init(&lat);
    for(int i = 0; i<4;i++){
        for(int j = 0; j<4;j++)
            ibz_copy(&(lat.basis[i][j]),&((*mat)[i][j]));
    }
    // transpose to be able to use algebra elemnts (which are then rows) more easily as vectors
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_copy(&(g->basis[i][j]),&((*mat)[j][i]));
        }
    }
    // G is the gram matrix: Symmetric, so only compute half
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < i+1; j++){
            ibz_vec_4_bilinear(&(g->G[i][j]),&(g->basis[i]),&(g->basis[j]),alg);
        }
    }
#ifndef NO_DIVISION
    ibz_t norm;
    ibz_init(&norm);
    ibz_copy(&norm,&(g->G[0][0]));
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < i+1; j++){
            ibz_gcd(&norm,&(g->G[i][j]),&norm);
        }
    }
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < i+1; j++){
            ibz_div_floor(&(g->G[i][j]),&den,&(g->G[i][j]),&norm);
            assert(ibz_is_zero(&den));
        }
    }
    ibz_finalize(&norm);
#endif
    //adapted from GSO_compute in sage version
    ibz_set(&den,1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < i; j++){
            ibq_set(&(g->r[i][j]),&(g->G[i][j]),&den);
            for(int k = 0; k < j+1; k++){
                ibq_mul(&tmp,&(g->mu[j][k]),&(g->r[i][k]));
                ibq_sub(&(g->r[i][j]),&(g->r[i][j]),&tmp);
            }
            ibq_div(&(g->mu[i][j]),&(g->r[i][j]),&(g->r[j][j]));
        }
        ibq_set(&(g->r[i][i]),&(g->G[i][i]),&den);
        for(int j = 0; j < i+1; j++){
                ibq_mul(&tmp,&(g->mu[i][j]),&(g->r[i][j]));
                ibq_sub(&(g->r[i][i]),&(g->r[i][i]),&tmp);
        }
    }
    ibz_finalize(&den);
    ibz_finalize(&num);
    ibq_finalize(&tmp);
    quat_lattice_finalize(&lat);
}

void quat_bkz_update_matrix(quat_bkz_matrix_t *g, int start, int end){
    ibz_t den, num;
    ibq_t tmp;
    ibq_init(&tmp);
    ibz_init(&den);
    ibz_init(&num);
    ibz_set(&den,1);
    for(int i = start; (i < end) & (i<4); i++){
        for(int j = 0; j < i; j++){
            ibq_set(&(g->r[i][j]),&(g->G[i][j]),&den);
            for(int k = 0; k < j+1; k++){
                ibq_mul(&tmp,&(g->mu[j][k]),&(g->r[i][k]));
                ibq_sub(&(g->r[i][j]),&(g->r[i][j]),&tmp);
            }
            ibq_div(&(g->mu[i][j]),&(g->r[i][j]),&(g->r[j][j]));
        }
        ibq_set(&(g->r[i][i]),&(g->G[i][i]),&den);
        for(int j = 0; j < i+1; j++){
                ibq_mul(&tmp,&(g->mu[i][j]),&(g->r[i][j]));
                ibq_sub(&(g->r[i][i]),&(g->r[i][i]),&tmp);
        }
    }
    ibz_finalize(&den);
    ibz_finalize(&num);
    ibq_finalize(&tmp);
}

void quat_bkz_size_reduce(quat_bkz_matrix_t *g, int start, int end){
    ibz_t delta_nf, gij,gjj, tmp, temp;
    ibz_init(&delta_nf);
    ibz_init(&gij);
    ibz_init(&gjj);
    ibz_init(&tmp);
    ibz_init(&temp);
    //for i in range(start,end):
    for(int i = start; i<end;i++){
        //for j in range(i-1,-1,-1):
        for(int j = i-1; j > -1; j--){
            //self.compute_GSO( start=0, end=i+1 )
            quat_bkz_update_matrix(g,0,i+1);
            //delta_nf = round( self.Mu[i][j] )
            ibq_round(&delta_nf,&(g->mu[i][j]));
            //gij, gjj = self.G[i][j], self.G[j][j]
            ibz_copy(&gij,&(g->G[i][j]));
            ibz_copy(&gjj,&(g->G[j][j]));
            //for kappa in range(i):
            for(int kappa = 0; kappa<i;kappa++){
                //self.G[i][kappa] -=  delta_nf*self.G[j][kappa] if kappa<=j else delta_nf*self.G[kappa][j]  #<b_i-delta_nf*b_j,b_kappa> = G_i_kappa - delta_nf*G_j_kappa
                int minkappaj = j*(kappa>j)+kappa*(1-(kappa>j));
                int maxkappaj = j*(kappa<=j)+kappa*(1-(kappa<=j));
                ibz_mul(&tmp,&delta_nf,&(g->G[maxkappaj][minkappaj]));
                ibz_sub(&(g->G[i][kappa]),&(g->G[i][kappa]),&tmp);
            }
            //for kappa in range(i+1,n):
            for(int kappa = i+1; kappa<4;kappa++){
                //self.G[kappa][i] -=  delta_nf*self.G[kappa][j] #<b_kappa,b_i-delta_nf*b_j> = G_kappa_i - delta_nf*G_kappa_j
                ibz_mul(&tmp,&delta_nf,&(g->G[kappa][j]));
                ibz_sub(&(g->G[kappa][i]),&(g->G[kappa][i]),&tmp);
            }
            //self.G[i][i] -= ( 2*delta_nf*gij-delta_nf*delta_nf*gjj ) #<b_i-delta_nf*b_j,b_i-delta_nf*b_j> = G_i_i -2*delta_nf*G_i_j + delta_nf^2*G_j_j
            ibz_mul(&tmp,&delta_nf,&gjj);
            ibz_set(&temp,2);
            ibz_mul(&temp,&temp,&gij);
            ibz_sub(&tmp,&temp,&tmp);
            ibz_mul(&tmp,&delta_nf,&tmp);
            ibz_sub(&(g->G[i][i]),&(g->G[i][i]),&tmp);
            //B[i] -= delta_nf *B[j] #update the basis
            ibz_mul(&tmp,&delta_nf,&(g->basis[j][0]));
            ibz_sub(&(g->basis[i][0]),&(g->basis[i][0]),&tmp);
            ibz_mul(&tmp,&delta_nf,&(g->basis[j][1]));
            ibz_sub(&(g->basis[i][1]),&(g->basis[i][1]),&tmp);
            ibz_mul(&tmp,&delta_nf,&(g->basis[j][2]));
            ibz_sub(&(g->basis[i][2]),&(g->basis[i][2]),&tmp);
            ibz_mul(&tmp,&delta_nf,&(g->basis[j][3]));
            ibz_sub(&(g->basis[i][3]),&(g->basis[i][3]),&tmp);
            //self.compute_GSO( start=i, end=i+1 )
            quat_bkz_update_matrix(g,i,i+1);
        }
    }
    ibz_finalize(&delta_nf);
    ibz_finalize(&gij);
    ibz_finalize(&gjj);
    ibz_finalize(&tmp);
    ibz_finalize(&temp);
}

void quat_bkz_update_after_lagrange(quat_bkz_matrix_t *g, const ibz_mat_2x2_t *U, int index){
    ibz_t g0,g1,g2, tmp, prod, temp;
    ibz_init(&g0);
    ibz_init(&g1);
    ibz_init(&g2);
    ibz_init(&tmp);
    ibz_init(&prod);
    ibz_init(&temp);
    //u0,u1 = U
    //g0, g1, g2 = self.G[i][i], self.G[i+1][i+1], self.G[i+1][i]
    ibz_copy(&g0,&(g->G[index][index]));
    ibz_copy(&g1,&(g->G[index+1][index+1]));
    ibz_copy(&g2,&(g->G[index+1][index]));
    //for ell in range(i):
    for(int ell = 0; ell < index; ell++){
        //tmp = U_[0][0]*self.G[i][ell] + U_[0][1]*self.G[i+1][ell]
        ibz_mul(&tmp,&((*U)[0][0]),&(g->G[index][ell]));
        ibz_mul(&prod,&((*U)[0][1]),&(g->G[index+1][ell]));
        ibz_add(&tmp,&tmp,&prod);
        //self.G[i+1][ell] = U_[1][0]*self.G[i][ell] + U_[1][1]*self.G[i+1][ell]
        ibz_mul(&prod,&((*U)[1][1]),&(g->G[index+1][ell]));
        ibz_mul(&(g->G[index+1][ell]),&((*U)[1][0]),&(g->G[index][ell]));
        ibz_add(&(g->G[index+1][ell]),&prod,&(g->G[index+1][ell]));
        //self.G[i][ell] = tmp
        ibz_copy(&(g->G[index][ell]),&tmp);
    }
    //for kappa in range(i+2,n):
    for(int kappa = index+2; kappa < 4; kappa++){
        //tmp = U_[0][0]*self.G[kappa][i] + U_[0][1]*self.G[kappa][i+1]
        ibz_mul(&tmp,&((*U)[0][0]),&(g->G[kappa][index]));
        ibz_mul(&prod,&((*U)[0][1]),&(g->G[kappa][index+1]));
        ibz_add(&tmp,&tmp,&prod);
        //self.G[kappa][i+1] = U_[1][0]*self.G[kappa][i] + U_[1][1]*self.G[kappa][i+1]
        ibz_mul(&prod,&((*U)[1][1]),&(g->G[kappa][index+1]));
        ibz_mul(&(g->G[kappa][index+1]),&((*U)[1][0]),&(g->G[kappa][index]));
        ibz_add(&(g->G[kappa][index+1]),&prod,&(g->G[kappa][index+1]));
        //self.G[kappa][i] = tmp
        ibz_copy(&(g->G[kappa][index]),&tmp);
    }
    //self.G[i][i] =      u0[0]*u0[0]*g0 + (u0[0]*u0[1]+u0[1]*u0[0])*g2+u0[1]*u0[1]*g1
    ibz_mul(&prod,&((*U)[0][0]),&((*U)[0][1]));
    ibz_set(&tmp,2);
    ibz_mul(&tmp,&g2,&tmp);
    ibz_mul(&tmp,&prod,&tmp);
    ibz_mul(&prod,&((*U)[0][1]),&((*U)[0][1]));
    ibz_mul(&prod,&prod,&g1);
    ibz_add(&tmp,&prod,&tmp);
    ibz_mul(&prod,&((*U)[0][0]),&((*U)[0][0]));
    ibz_mul(&prod,&prod,&g0);
    ibz_add(&(g->G[index][index]),&prod,&tmp); //self.G[i][i] =      u0[0]*u0[0]*g0 + (u0[0]*u0[1]+u0[1]*u0[0])*g2+u0[1]*u0[1]*g1
    //self.G[i+1][i+1] =  u1[0]*u1[0]*g0 + (u1[0]*u1[1]+u1[1]*u1[0])*g2+u1[1]*u1[1]*g1
    ibz_mul(&prod,&((*U)[1][0]),&((*U)[1][1]));
    ibz_set(&tmp,2);
    ibz_mul(&tmp,&g2,&tmp);
    ibz_mul(&tmp,&prod,&tmp);
    ibz_mul(&prod,&((*U)[1][1]),&((*U)[1][1]));
    ibz_mul(&prod,&prod,&g1);
    ibz_add(&tmp,&prod,&tmp);
    ibz_mul(&prod,&((*U)[1][0]),&((*U)[1][0]));
    ibz_mul(&prod,&prod,&g0);
    ibz_add(&(g->G[index+1][index+1]),&prod,&tmp); //self.G[i+1][i+1] =  u1[0]*u1[0]*g0 + (u1[0]*u1[1]+u1[1]*u1[0])*g2+u1[1]*u1[1]*g1
    //self.G[i+1][i] =    u1[0]*u0[0]*g0 + (u1[0]*u0[1]+u1[1]*u0[0])*g2+u1[1]*u0[1]*g1
    ibz_mul(&prod,&((*U)[1][0]),&((*U)[0][1]));
    ibz_mul(&tmp,&((*U)[0][0]),&((*U)[1][1]));
    ibz_add(&tmp,&prod,&tmp);
    ibz_mul(&tmp,&g2,&tmp);
    ibz_mul(&prod,&((*U)[1][1]),&((*U)[0][1]));
    ibz_mul(&prod,&prod,&g1);
    ibz_add(&tmp,&prod,&tmp);
    ibz_mul(&prod,&((*U)[1][0]),&((*U)[0][0]));
    ibz_mul(&prod,&prod,&g0);
    ibz_add(&(g->G[index+1][index]),&prod,&tmp); //self.G[i+1][i] =    u1[0]*u0[0]*g0 + (u1[0]*u0[1]+u1[1]*u0[0])*g2+u1[1]*u0[1]*g1
    ibz_finalize(&g0);
    ibz_finalize(&g1);
    ibz_finalize(&g2);
    ibz_finalize(&tmp);
    ibz_finalize(&temp);
    ibz_finalize(&prod);
}

void quat_bkz_lagrange_reduction_gram(ibz_mat_2x2_t *U, ibq_mat_2x2_t *G, int lagrange_tours){
    ibq_t mu;
    ibz_t tmp, muround, prod,one;
    ibq_t tmpq, prodq;
    ibq_init(&tmpq);
    ibq_init(&prodq);
    ibq_init(&mu);
    ibz_init(&tmp);
    ibz_init(&one);
    ibz_init(&prod);
    ibz_init(&muround);
    ibz_set(&one,1);
    ibq_div(&mu,&((*G)[1][0]),&((*G)[0][0]));
    ibz_set(&((*U)[0][0]),1);
    ibz_set(&((*U)[1][1]),1);
    ibz_set(&((*U)[0][1]),0);
    ibz_set(&((*U)[1][0]),0);
    for(int cntr = 0; cntr < lagrange_tours;cntr++){
        ibq_div(&mu,&((*G)[1][0]),&((*G)[0][0]));
        ibq_round(&muround,&mu);
        ibq_set(&mu,&muround,&one);
        ibz_set(&tmp,2);
        ibz_mul(&tmp,&tmp,&muround);
        ibq_set(&tmpq,&tmp,&one);
        ibq_mul(&tmpq,&tmpq,&((*G)[1][0]));
        ibz_mul(&prod,&muround,&muround);
        ibq_set(&prodq,&prod,&one);
        ibq_mul(&prodq,&prodq,&((*G)[0][0]));
        ibq_sub(&tmpq,&tmpq,&prodq);
        ibq_sub(&((*G)[1][1]),&((*G)[1][1]),&tmpq); //G[1][1] -= ( 2*muround*G[1][0]-muround**2 * G[0][0] )
        ibq_mul(&prodq,&mu,&((*G)[0][0]));
        ibq_sub(&((*G)[1][0]),&((*G)[1][0]),&prodq);
        ibz_set(&tmp,-1);
        ibq_set(&mu,&muround,&tmp);
        ibz_mul(&tmp,&muround,&((*U)[0][0]));
        ibz_sub(&((*U)[1][0]),&((*U)[1][0]),&tmp);
        ibz_mul(&tmp,&muround,&((*U)[0][1]));
        ibz_sub(&((*U)[1][1]),&((*U)[1][1]),&tmp); //U[1] -= muround*(U[0])
        ibq_copy(&tmpq,&((*G)[0][0]));
        ibq_copy(&((*G)[0][0]),&((*G)[1][1]));
        ibq_copy(&((*G)[1][1]),&tmpq);//G[0][0], G[1][1] = G[1][1], G[0][0] #swap
        ibz_copy(&tmp,&((*U)[0][0]));
        ibz_copy(&((*U)[0][0]),&((*U)[1][0]));
        ibz_copy(&((*U)[1][0]),&tmp);
        ibz_copy(&tmp,&((*U)[0][1]));
        ibz_copy(&((*U)[0][1]),&((*U)[1][1]));
        ibz_copy(&((*U)[1][1]),&tmp);//U[0], U[1] = U[1], U[0]
    }
    ibq_div(&mu,&((*G)[1][0]),&((*G)[0][0]));
    ibz_conditional_swap(&((*U)[0][0]), &((*U)[1][0]), ibq_cmp(&((*G)[1][1]) , &((*G)[0][0])) < 0);
    ibz_conditional_swap(&((*U)[0][1]), &((*U)[1][1]), ibq_cmp(&((*G)[1][1]) , &((*G)[0][0])) < 0);
    ibq_finalize(&mu);
    ibz_finalize(&tmp);
    ibq_finalize(&tmpq);
    ibz_finalize(&prod);
    ibq_finalize(&prodq);
    ibz_finalize(&muround);
    ibz_finalize(&one);
}

void quat_lattice_bkz(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, int bkz_tours, int lagrange_tours, const quat_alg_t *alg){
    quat_bkz_matrix_t g;
    ibz_mat_2x2_t U;
    ibq_mat_2x2_t M;
    ibq_t tmp;
    ibq_init(&tmp);
    ibz_mat_2x2_init(&U);
    ibq_mat_2x2_init(&M);
    quat_bkz_matrix_init(&g);
//G = GSO_obj( B )
//G.compute_GSO(0,B.nrows())
//n = len(list(B))
    quat_bkz_matrix_set(&g,&(lattice->basis),alg);
//for t in range(lll_max_tours):
    for(int t =0; t < bkz_tours;t++){
        //for i in range(n-1):
        for(int i = 0; i < 3; i++){
            //B = G.size_reduce(B,i+1,i+2)
            quat_bkz_size_reduce(&g,i+1,i+2);
            //M =  [ [G.Rr[i][i]], [G.Mu[i+1][i]*G.Rr[i][i],G.Mu[i+1][i]**2*G.Rr[i][i] + G.Rr[i+1][i+1]] ]
            ibq_copy(&(M[0][0]),&(g.r[i][i]));
            ibq_mul(&(M[1][0]),&(g.mu[i+1][i]),&(g.r[i][i]));
            ibq_mul(&tmp,&(M[1][0]),&(g.mu[i+1][i]));
            ibq_add(&(M[1][1]),&tmp,&(g.r[i+1][i+1]));

            //U = lagrange_reduction_gram(M)
            quat_bkz_lagrange_reduction_gram(&U,&M,lagrange_tours);
            //B[i:i+2]=U*B[i:i+2]
            //set g.basis i, i+1 to UB i,i+1: B[i:i+2]=U*B[i:i+2];
            ibz_vec_2x4_mul(&(g.basis[i]),&(g.basis[i+1]),&U);
            //G.update_after_svp_oracle_rank_2(U,i)
            quat_bkz_update_after_lagrange(&g,&U,i);
            //G.compute_GSO(start=i, end=i+2);
            quat_bkz_update_matrix(&g,i,i+2);
            //B = G.size_reduce(B,i,i+2)
            quat_bkz_size_reduce(&g,i,i+2);
        }
    }
    ibz_mat_4x4_transpose(red,&(g.basis));
    quat_bkz_matrix_finalize(&g);
    ibz_mat_2x2_finalize(&U);
    ibq_mat_2x2_finalize(&M);
    ibq_finalize(&tmp);
}