#include "lll/lll.h"
#include "bkz/bkz.h"
#include "bench.h"
#include "inputs.h"
#include "bkz_constants.h"
#include <stdio.h>
#include <inttypes.h>

#ifndef LLL_BENCH
#define LLL_BENCH 0
#endif

#ifndef BKZ_BENCH
#define BKZ_BENCH 1
#endif

#if defined(DATA_FILE) && defined(LATTICE_DATA_NUMBER)
#include "read_data.h"

uint64_t quat_bench_bkz_from_file(){
    printf("Run benchmarks for BKZ on %d lattices with %d iterations:\n",LATTICE_DATA_NUMBER, BENCHMARK_ITERATIONS);
    printf("BKZ parameters: Lagrange tours %d, outer tours %d, integer bitsize %d\n",LAGRANGE_TOURS, BKZ_TOURS,64*LIMBNUM);
    printf("Lattice file: %s\n\n",DATA_FILE);
    ibz_mat_4x4_t red;
    uint64_t cycles_list[LATTICE_DATA_NUMBER];
    uint64_t start;
    quat_lattice_t lat;
    quat_alg_t alg;
    quat_alg_init_set(&alg,&QUATALG_P);
    quat_lattice_init(&lat);
    ibz_mat_4x4_init(&red);
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++)
        cycles_list[i] = 0;

    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        read_lattice_from_file(&lat,i);
        start= cpucycles();
        for (int j = 0; j < BENCHMARK_ITERATIONS; j++){
            quat_lattice_bkz(&red, &lat,BKZ_TOURS,LAGRANGE_TOURS,&alg);
        }
        cycles_list[i] = (cycles_list[i]+cpucycles()-start) /BENCHMARK_ITERATIONS;
    }
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
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


uint64_t quat_bench_lll_from_file(){
    printf("Run benchmarks for LLL on %d lattices with %d iterations:\n",LATTICE_DATA_NUMBER, BENCHMARK_ITERATIONS);
    printf("LLL parameters: integer bitsize %d\n",64*LIMBNUM);
    printf("Lattice file: %s\n\n",DATA_FILE);
    ibz_mat_4x4_t red;
    uint64_t cycles_list[LATTICE_DATA_NUMBER];
    uint64_t start;
    quat_lattice_t lat;
    quat_alg_t alg;
    quat_alg_init_set(&alg,&QUATALG_P);
    quat_lattice_init(&lat);
    ibz_mat_4x4_init(&red);
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++)
        cycles_list[i] = 0;

    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        read_lattice_from_file(&lat,i);
        start= cpucycles();
        for (int j = 0; j < BENCHMARK_ITERATIONS; j++){
            quat_lattice_lll(&red, &lat,&QUATALG_P,0);
        }
        cycles_list[i] = (cycles_list[i]+cpucycles()-start)/BENCHMARK_ITERATIONS;
    }
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
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
    if(BKZ_BENCH != 0)
        quat_bench_bkz_from_file();
    if(LLL_BENCH != 0)
        quat_bench_lll_from_file();
    if((BKZ_BENCH==0)&&(LLL_BENCH==0))
        printf("No benchmarks to run: Pleae enable benchmarks for LLL or BKZ\n");
    return(0);
}

#elif defined(LATTICE_NUMBER)

uint64_t quat_bench_bkz_on_list(int iterations, int lattice_number,const quat_lattice_t lat[lattice_number],const quat_alg_t *alg){
    printf("Run benchmarks for BKZ on %d lattices with %d iterations:\n",lattice_number,iterations);
    printf("BKZ parameters: Lagrange tours %d, outer tours %d, integer bitsize %d\n",LAGRANGE_TOURS, BKZ_TOURS,64*LIMBNUM);
    ibz_mat_4x4_t red;
    uint64_t cycles_list[LATTICE_NUMBER];
    uint64_t start;
    for(int i = 0; i<lattice_number;i++)
        cycles_list[i] = 0;
    ibz_mat_4x4_init(&red);
    for (int j = 0; j < iterations; j++){
        for(int i = 0; i<lattice_number;i++){
            start= cpucycles();
            quat_lattice_bkz(&red, &(lat[i]),BKZ_TOURS,LAGRANGE_TOURS,alg);
            cycles_list[i] = cycles_list[i]+cpucycles()-start;
        }
    }
    ibz_mat_4x4_finalize(&red);
    printf("Average cycles per lattice:\n");
    for(int i = 0; i<lattice_number;i++){
        printf("%d. %"PRIu64"\n", i,(cycles_list[i]/iterations));
    }
    uint64_t avg = 0;
    for(int i = 0; i<lattice_number;i++){
        avg = avg + cycles_list[i]/iterations;
    }
    avg = avg/((0==lattice_number)+lattice_number);
    return(avg);
}


uint64_t quat_bench_lll_on_list(int iterations, int lattice_number,const quat_lattice_t lat[lattice_number],const quat_alg_t *alg){
    printf("Run benchmarks for LLL on %d lattices with %d iterations:\n",lattice_number,iterations);
    printf("LLL parameters:  integer bitsize %d\n",64*LIMBNUM);
    ibz_mat_4x4_t red;
    uint64_t start;
    uint64_t cycles_list[LATTICE_NUMBER];
    for(int i = 0; i<lattice_number;i++)
        cycles_list[i] = 0;
    ibz_mat_4x4_init(&red);
    for (int j = 0; j < iterations; j++){
        for(int i = 0; i<lattice_number;i++){
            start = cpucycles();
            quat_lattice_lll(&red, &(lat[i]),&(alg->p),0);
            cycles_list[i] = cycles_list[i]+cpucycles()-start;
        }
    }
    ibz_mat_4x4_finalize(&red);
    printf("average cycles per lattice:\n");
    for(int i = 0; i<lattice_number;i++){
        printf("%d. %"PRIu64"\n", i,(cycles_list[i]/iterations));
    }
    uint64_t avg = 0;
    for(int i = 0; i<lattice_number;i++){
        avg = avg + cycles_list[i]/iterations;
    }
    avg = avg/((0==lattice_number)+lattice_number);;
    return(avg);
}

//run test on input 1, ct runnner on input 0 or no input
int main(int argc, char* argv[]){
    //maybe take parametrers from command line and make LLL bench optional
    quat_alg_t alg;
    quat_alg_init_set(&alg,&QUATALG_P);
    if(BKZ_BENCH!=0)
        quat_bench_bkz_on_list(BENCHMARK_ITERATIONS, LATTICE_NUMBER,LATTICE_LIST,&alg);
    if(LLL_BENCH!=0)
        quat_bench_lll_on_list(BENCHMARK_ITERATIONS, LATTICE_NUMBER,LATTICE_LIST,&alg);
    if((BKZ_BENCH==0)&&(LLL_BENCH==0))
        printf("No benchmarks to run: Pleae enable benchmarks for LLL or BKZ\n");
    quat_alg_finalize(&alg);
}

#else
int main(int argc, char* argv[]){
    printf("input.h not correctly set for bench.c\n");
    return(0);
}
#endif