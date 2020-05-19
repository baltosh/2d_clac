//
// Created by nzvi on 19.05.2020.
//

#ifndef INC_2D_CALC_PARAMETERS_H
#define INC_2D_CALC_PARAMETERS_H

#include "global.h"

double get_field_ro(int i, int j, double x, double y);

double get_field_ru(int i, int j, double x, double y);

double get_field_rv(int i, int j, double x, double y);

double get_field_re(int i, int j, double x, double y);

double get_field_rc(int i, int j, double x, double y, int e);

void cons_to_prim(int i, int j, double x, double y, prim_t *prim);

double get_component_cp(int id);

double get_component_M(int id);

double get_cell_cp(int i, int j);

double get_cell_M(int i, int j);

#endif //INC_2D_CALC_PARAMETERS_H
