//
// Created by nzvi on 19.05.2020.
//

#ifndef INC_2D_CALC_PHYSICS_H
#define INC_2D_CALC_PHYSICS_H

#include <math.h>
#include "chemistry.h"
#include "global.h"

void get_left_boundary(prim_t *par_m, prim_t *par_p);

void get_right_boundary(prim_t *par_m, prim_t *par_p);

void get_bottom_boundary(prim_t *par_m, prim_t *par_p);

void get_top_boundary(prim_t *par_m, prim_t *par_p);

void calc_prim_params(prim_t *par);

#endif //INC_2D_CALC_PHYSICS_H
