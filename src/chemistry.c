//
// Created by nzvi on 19.05.2020.
//

#include "include/chemistry.h"

double k1 = 0.503;
double k2 = 0.073;

double (*reactionSpeeds[COMPONENTS_COUNT])(double[]) = {reactionSpeed0,
                                                        reactionSpeed1,
                                                        reactionSpeed2,
                                                        reactionSpeed3,
                                                        reactionSpeed4};

//todo change
double reactionSpeed0(double rc[]) {
//    double k1 = 928.979;
//    double k2 = 244.927;
    return -k1 * rc[0] - 2 * k2 * rc[0] * rc[0];
};

double reactionSpeed1(double rc[]) {
//    double k1 = 928.979;
//    double k2 = 244.927;
    return k1 * rc[0] + k2 * rc[0] * rc[0];
};

double reactionSpeed2(double rc[]) {
//    double k1 = 928.979;
//    double k2 = 244.927;
    return k1 * rc[0];
};

double reactionSpeed3(double rc[]) {
//    double k1 = 928.979;
//    double k2 = 244.927;
    return 2 * k2 * rc[0] * rc[0];
};