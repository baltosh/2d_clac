//
// Created by nzvi on 19.05.2020.
//


#include "include/chemistry.h"

//double k1 = 0.503;
//double k2 = 0.073;

double (*reactionSpeeds[COMPONENTS_COUNT])(double[], double) = {reactionSpeed0,
                                                        reactionSpeed1,
                                                        reactionSpeed2,
                                                        reactionSpeed3};

double arrheniusEq(double A, double E, double T)
{
    return A * exp(-E / (GAS_CONSTANT * T));
}

double k1(double T){
    double A = 1.08e+16;
    double E = 2.5e+05;
    return arrheniusEq(A, E, T);
}

double k2(double T){
    double A = 3.16e+16;
    double E = 2.7e+05;
    return arrheniusEq(A, E, T);
}

//todo change
double reactionSpeed0(double rc[], double T) {
    return -k1(T) * rc[0] - 2 * k2(T) * rc[0] * rc[0];
}

double reactionSpeed1(double rc[], double T) {
    return k1(T) * rc[0] + k2(T) * rc[0] * rc[0];
}

double reactionSpeed2(double rc[], double T) {
    return k1(T) * rc[0];
}

double reactionSpeed3(double rc[], double T) {
    return 2 * k2(T) * rc[0] * rc[0];
}

double get_component_cp(int id) {
    //этан, этилен, водород, метан
    double cps[COMPONENTS_COUNT] = {469.64/**1014.16*/, 572.57, 14274.97, 925.55};
    return cps[id];
}

double get_component_M(int id) {
    //этан, этилен, водород, метан
    double Ms[COMPONENTS_COUNT] = {0.03007012/**0.02869409*/, 0.02805418, 0.00201594, 0.01604303};
    return Ms[id];
}