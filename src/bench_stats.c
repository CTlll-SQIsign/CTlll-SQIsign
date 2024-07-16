#include "lll/lll.h"
#include "bkz/bkz.h"
#include "bench.h"
#include "inputs.h"
#include "bkz_constants.h"
#include "test_helpers/test_code.h"
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>

#ifndef NTESTS
#define NTESTS 0
#endif

#if defined(DATA_FILE) && defined(LATTICE_DATA_NUMBER)
#include "read_data.h"

void quat_bench_bkz_from_file(){
    printf("Run benchmarks for %d lattices with %d iterations:\n",LATTICE_DATA_NUMBER, BENCHMARK_ITERATIONS);
    printf("BKZ parameters: Lagrange tours %d, outer tours %d, integer bitsize %d\n",LAGRANGE_TOURS, BKZ_TOURS,64*LIMBNUM);
    printf("Lattice file: %s\n\n",DATA_FILE);
    ibz_mat_4x4_t red;
    quat_lattice_t lat;
    quat_alg_t alg;
    quat_alg_init_set(&alg,&QUATALG_P);
    quat_lattice_init(&lat);
    ibz_mat_4x4_init(&red);

    uint64_t start;
    uint64_t cycles = 0;


    FILE *f;
    char* filename;
    filename = malloc(strlen(DATA_FILE)+strlen("bench_")+1);
    strcpy(filename, "bench_");
    strcat(filename, DATA_FILE);
    f = fopen(filename,"w");


    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        read_lattice_from_file(&lat,i);
        start = cpucycles();
        quat_lattice_bkz(&red, &lat,BKZ_TOURS,LAGRANGE_TOURS,&alg);
        cycles = cpucycles()-start;
        fprintf(f,"%"PRIu64"\n", cycles);
        if(i%100==0){
            printf("%d lattices are bkz reduced \n ", i);
        }
    }

    fclose(f);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);

}

int main(int argc, char* argv[]){

    if (NTESTS>LATTICE_DATA_NUMBER){
        printf("NTESTS = %d > LATTICE_DATA_NUMBER = %d \n", NTESTS,LATTICE_DATA_NUMBER);
    }
    quat_lattice_t lat;
    quat_alg_t alg;

    quat_lattice_init(&lat);
    quat_alg_init_set(&alg,&QUATALG_P);


    //run a few tests on the input lattice to check that the size of ints was set appropriately
    //if any of the test fails, benchmarking appears useless
    int res = 0;
    int bkz_res = 0;
    int ideal = 1;


    for (int i = 0; i<NTESTS; i++){
        read_lattice_from_file(&lat,i);
        res = res | quat_test_bkz_on_lattice(&lat,&alg,BKZ_TOURS,LAGRANGE_TOURS,0);
    }
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);

    if (res){
        printf("Test bkz_on_list failed\n");
    }
    else{
        printf("Test bkz_on_list passed\n");
        quat_bench_bkz_from_file();
    }
    return(0);

}

#else
int main(int argc, char* argv[]){
    printf("input.h not correctly set for bench_stats.c\n");
    return(0);
}
#endif