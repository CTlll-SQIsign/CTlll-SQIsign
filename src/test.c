#include "bkz/bkz.h"
#include "lll/lll.h"
#include "inputs.h"
#include "bkz_constants.h"
#include "test_helpers/test_code.h"
#include <stdio.h>
#include <assert.h>

#ifndef LLL_TEST
#define LLL_TEST 0
#endif
#ifndef BKZ_TEST
#define BKZ_TEST 1
#endif

#ifndef PRINT_FLAG
#define PRINT_FLAG 0
#endif

#if defined(LATTICE_NUMBER)

int quat_test_on_list(int lattice_number,const quat_lattice_t lat[lattice_number], const quat_alg_t *alg){
    int lll_res = 0;
    int bkz_res = 0;
    int res = 0;
    for(int i = 0; i<lattice_number;i++){
        if(BKZ_TEST!=0){
            bkz_res = quat_test_bkz_on_lattice(&((lat)[i]),alg,BKZ_TOURS,LAGRANGE_TOURS,PRINT_FLAG);
            if(!bkz_res) printf("BKZ test %d passed\n",i+1); else printf("BKZ test %d failed, continue test\n",i+1);
        }
        if(LLL_TEST!=0){
            lll_res = quat_test_lll_on_lattice(&((lat)[i]),alg);
            if(!lll_res) printf("LLL test %d passed\n",i+1); else printf("LLL test %d failed, continue test\n",i+1);
        }
        res = res | lll_res | bkz_res;
    }
    if(res){
        printf("Test on_list failed\n");
    }
    return(res);
}


int main(int argc, char* argv[]){
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set(&alg,&QUATALG_P);
    printf("Run reduction tests on %d lattices\n",LATTICE_NUMBER);

    res = res | quat_test_on_list(LATTICE_NUMBER,LATTICE_LIST,&alg);
    if((BKZ_TEST==0)&&(LLL_TEST==0)){
        printf("No test to be run: Please enable test for LLL or BKZ\n");
    } else {
        if(res)
            printf("Some test failed\n");
        else
            printf("All tests passed\n");
    }
    quat_alg_finalize(&alg);
}

#elif (defined(DATA_FILE)&&defined(LATTICE_DATA_NUMBER))
#include "read_data.h"

int quat_test_from_file(ibz_t *norm, quat_lattice_t *lat, const quat_alg_t *alg){
    int res = 0;
    int lll_res = 0;
    int bkz_res = 0;
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        read_lattice_from_file(lat,i);
        if(BKZ_TEST!=0){
            bkz_res = quat_test_bkz_on_lattice(lat,alg,BKZ_TOURS,LAGRANGE_TOURS, PRINT_FLAG);
            if(!PRINT_FLAG){
                if(!bkz_res) printf("BKZ test %d passed\n",i+1); else printf("BKZ test %d failed, continue test\n",i+1);
            }
        }
        if(LLL_TEST!=0){
            lll_res = quat_test_lll_on_lattice(lat,alg);
            if(!lll_res) printf("LLL test %d passed\n",i+1); else printf("LLL test %d failed, continue test\n",i+1);
        }
        res = res | lll_res | bkz_res;
    }
    if(res){
        printf("Test from_file failed\n");
    }
    return(res);
}

int main(int argc, char* argv[]){
    int res = 0;
    quat_alg_t alg;
    ibz_t norm;
    quat_lattice_t lat;
    ibz_init(&norm);
    quat_lattice_init(&lat);
    quat_alg_init_set(&alg,&QUATALG_P);
    printf("Run reduction tests on %d lattices\n",LATTICE_DATA_NUMBER);
    res = quat_test_from_file(&norm,&lat,&alg);
    if((BKZ_TEST==0)&&(LLL_TEST==0)){
        printf("No test to be run: Please enable test for LLL or BKZ\n");
    } else {
        if(res)
            printf("Some test failed\n");
        else
            printf("All tests passed\n");
    }
    quat_lattice_finalize(&lat);
    ibz_finalize(&norm);
    quat_alg_finalize(&alg);
}

#else
int main(int argc, char* argv[]){
    printf("input.h not correctly set for bench.h\n");
    return(0);
}
#endif