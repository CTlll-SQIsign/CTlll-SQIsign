#include "nct_intbig_lll/lll.h"
#include "ct_helpers/ct_helpers.h"
#include "lll/lll.h"
#include "bench.h"
#include "bkz_constants.h"
#include "inputs.h"
#include "read_data.h"
#include <stdio.h>
#include <inttypes.h>

void nct_quat_lattice_init(nct_quat_lattice_t *lat){
    nct_ibz_mat_4x4_init(&(*lat).basis);
    nct_ibz_init(&(*lat).denom);
    nct_ibz_set(&(*lat).denom, 1);
}
void nct_quat_lattice_finalize(nct_quat_lattice_t *lat){
    nct_ibz_finalize(&(*lat).denom);
    nct_ibz_mat_4x4_finalize(&(*lat).basis);
}

void ct_to_nct_integer_translation(nct_ibz_t *nct, const ibz_t *ct){
    char str_ct[PRINTSIZE];
    ibz_t a;
    ibz_init(&a);
    ibz_abs(&a,ct);
    ibz_to_str(str_ct,&a);
    nct_ibz_set_from_str(nct,str_ct,10);
    if(ibz_cmp(&a,ct)!=0) nct_ibz_neg(nct,nct);
    ibz_finalize(&a);
}

//write init and finit for nct lattices and matrices
void nct_quat_lattice_translate(nct_quat_lattice_t* nct_lat, const quat_lattice_t *lat){
    for (int i = 0; i<4;i++){
        for (int j = 0; j<4;j++){
            ct_to_nct_integer_translation(&(nct_lat->basis[i][j]), &(lat->basis[i][j]));
        }
    }
    ct_to_nct_integer_translation(&(nct_lat->denom), &(lat->denom));
}

//write translate for lattices

uint64_t quat_bench_nct_lll_from_file(){
    printf("Run benchmarks for LLL with variable-size integers on %d lattices with %d iterations:\n",LATTICE_DATA_NUMBER, BENCHMARK_ITERATIONS);
    printf("Lattice file: %s\n\n",DATA_FILE);
    uint64_t cycles_list[LATTICE_DATA_NUMBER];
    uint64_t start;
    nct_ibz_mat_4x4_t red;
    nct_ibz_t palg;
    quat_lattice_t lat;
    nct_quat_lattice_t nct_lat;
    nct_ibz_mat_4x4_init(&red);
    quat_lattice_init(&lat);
    nct_quat_lattice_init(&nct_lat);
    nct_ibz_init(&palg);
    ct_to_nct_integer_translation(&palg,&QUATALG_P);
    
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++)
        cycles_list[i] = 0;

    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        read_lattice_from_file(&lat,i);
        nct_quat_lattice_translate(&nct_lat,&lat);
        //need to translate lat into right lattice type
        start= cpucycles();
        for (int j = 0; j < BENCHMARK_ITERATIONS; j++){
            nct_quat_lattice_lll(&red, &nct_lat,&palg,0);
        }
        cycles_list[i] = (cycles_list[i]+cpucycles()-start)/BENCHMARK_ITERATIONS;
    }
    nct_ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    nct_quat_lattice_finalize(&nct_lat);
    nct_ibz_finalize(&palg);
    printf("Average cycles per lattice:\n");
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        printf("%d. %"PRIu64"\n", i,cycles_list[i]);
    }
    uint64_t avg = 0;
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        avg = avg + cycles_list[i];
    }
    avg = avg/((0==LATTICE_DATA_NUMBER)+LATTICE_DATA_NUMBER);
    return(avg);
}

int main(int argc, char* argv[]){
    quat_bench_nct_lll_from_file();
    return(0);
}
