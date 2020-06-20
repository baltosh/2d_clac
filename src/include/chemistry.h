//
// Created by nzvi on 19.05.2020.
//

#ifndef INC_2D_CALC_CHEMISTRY_H
#define INC_2D_CALC_CHEMISTRY_H

#include "global.h"

double reactionSpeed0(double[], double);

double reactionSpeed1(double[], double);

double reactionSpeed2(double[], double);

double reactionSpeed3(double[], double);

double reactionSpeed4(double[], double);

extern double (*reactionSpeeds[COMPONENTS_COUNT])(double[], double);

double k1(double);

double k2(double);

double arrheniusEq(double, double, double);

#endif //INC_2D_CALC_CHEMISTRY_H
