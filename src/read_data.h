#include <stdio.h>
#include "ct_intbig/ct_intbig.h"

static inline void read_lattice_from_file(quat_lattice_t *lat, int lattice_index){
    FILE *f;
    f = fopen(DATA_FILE,"r");
    fseek(f,sizeof(quat_lattice_t)*lattice_index,SEEK_SET);
    fread(lat,sizeof(quat_lattice_t),1,f);
    fclose(f);
}