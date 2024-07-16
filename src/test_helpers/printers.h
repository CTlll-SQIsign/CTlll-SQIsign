#ifndef NCT_Printers_H
#define NCT_Printers_H

#include "nct_helpers.h"


static inline void nct_ibz_mat_2x2_print(const nct_ibz_mat_2x2_t *mat){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            nct_ibz_printf("%Zd, ",&((*mat)[i][j]));
        }
        nct_ibz_printf("\n");
    }
}

static inline void nct_ibz_mat_4x4_print(const nct_ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_printf("%Zd, ",&((*mat)[i][j]));
        }
        nct_ibz_printf("\n");
    }
}

static inline void nct_ibq_mat_4x4_print(const nct_ibq_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_printf("%Qd, ",&((*mat)[i][j]));
        }
        nct_ibz_printf("\n");
    }
}

static inline void nct_quat_lattice_print(const nct_quat_lattice_t *lat)
{
    nct_ibz_printf("lattice\n");
    nct_ibz_printf("denominator: %Zd\n", (lat->denom));
    nct_ibz_printf("basis: ");
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            nct_ibz_printf("%Zd ", &((lat->basis)[i][j]));
        }
        nct_ibz_printf("\n       ");
    }
    nct_ibz_printf("\n");
}
#endif