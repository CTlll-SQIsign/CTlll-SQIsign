#include "ct_intbig/test_ct_intbig.h"
#include "test_helpers/test_test_helpers.h"
#include "bkz/test_bkz.h"
#include "lll/test_lll.h"
#include <stdio.h>

int main(){
    int res = 0;
    printf("Run tests\n");
    res = res | ibz_tests();
    res = res | ibq_tests();
    res = res | tests_helper_tests_all();
    res = res | bkz_tests();
    res = res | lll_tests();
    if(res !=0){
        printf("Some tests failed\n");
    } else{
        printf("All passed\n");
    }
    return(res);
}