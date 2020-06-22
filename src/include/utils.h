//
// Created by nzvi on 19.05.2020.
//

#ifndef INC_2D_CALC_UTILS_H
#define INC_2D_CALC_UTILS_H

#include "global.h"
#include <string.h>
#include "parameters.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

void zero_r();

void inverse_matr(double a_src[BASE_FN_COUNT][BASE_FN_COUNT],
                  double am[BASE_FN_COUNT][BASE_FN_COUNT], int N);

void mult_matr_vec(double matr[BASE_FN_COUNT][BASE_FN_COUNT], double vec[BASE_FN_COUNT],
                   double res[BASE_FN_COUNT]);

double _minmod(double a, double b, double c);

void save_vtk(int step);

#endif //INC_2D_CALC_UTILS_H
