//
// Created by nzvi on 19.05.2020.
//

#ifndef INC_2D_CALC_PARAMETERS_H
#define INC_2D_CALC_PARAMETERS_H

#include <math.h>
#include <stdio.h>
#include "basis.h"
#include "chemistry.h"
#include "global.h"

double get_field_ro(int i, int j, double x, double y);

double get_field_ru(int i, int j, double x, double y);

double get_field_rv(int i, int j, double x, double y);

double get_field_re(int i, int j, double x, double y);

double get_field_rc(int i, int j, double x, double y, int i_com);

double get_field_T(int i, int j, double x, double y);

void cons_to_prim(int i, int j, double x, double y, prim_t *prim);

double get_cell_cp(int i, int j);

double get_cell_M(int i, int j);

#endif //INC_2D_CALC_PARAMETERS_H
