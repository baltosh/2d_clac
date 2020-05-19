//
// Created by nzvi on 19.05.2020.
//

#include "include/physics.h"
#include <math.h>

void calc_prim_params(double cp, double m_mol, prim_t *par) {
    double cv = cp - GAS_CONSTANT / m_mol;
    double gam = cp / cv;

    par->e = par->p / (par->r * (gam - 1.0));
    par->e_tot = par->e + (par->u * par->u + par->v * par->v) * 0.5;
    par->cz = sqrt(gam * par->p / par->r);
    par->t = par->e / cv;
}

void get_left_boundary(prim_t *par_m, prim_t *par_p) {
    double cp = 1014.16;
    double m_mol = 0.02869409;

    par_m->r = par_p->r;
    par_m->u = -par_p->u;
    par_m->v = par_p->v;
    par_m->p = par_p->p;
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_m->c[i_comp] = par_p->c[i_comp];
    }

//    double cp = 469.64;
//    double m_mol = 0.03007012;
//    double cv = cp - R_GAS / m_mol;
//    double gam = cp / cv;
//
//    par_m->r = 1.293;
//    par_m->u = 10;
//    par_m->v = 0.0;
//    par_m->p = par_m->r * R_GAS * T / m_mol;
    calc_prim_params(cp, m_mol, par_m);
}

void get_right_boundary(prim_t *par_m, prim_t *par_p) {
    double cp = 1014.16;
    double m_mol = 0.02869409;

    par_p->r = par_m->r;
    par_p->u = -par_m->u;
    par_p->v = par_m->v;
    par_p->p = par_m->p;
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_p->c[i_comp] = par_m->c[i_comp];
    }

//    double cp = 469.64;
//    double m_mol = 0.03007012;
//    double cv = cp - R_GAS / m_mol;
//    double gam = cp / cv;
//
//    par_p->r = par_m->r;
//    par_p->u = par_m->u;
//    par_p->v = par_m->v;
//    par_p->p = par_m->p;
    calc_prim_params(cp, m_mol, par_p);
}

void get_bottom_boundary(prim_t *par_m, prim_t *par_p) {

    double cp = 1014.16;
    double m_mol = 0.02869409;
    double cv = cp - GAS_CONSTANT / m_mol;
    double gam = cp / cv;

    par_m->r = 12.090;
    par_m->u = 0.0;
    par_m->v = 97.76;
    par_m->p = 2.152e+5;
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_m->c[i_comp] = 0.0;
    }
    par_m->c[0] = 1.0;

//    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
//        par_m->c[i_comp] = par_p->c[i_comp];
//    }
//    double cp = 469.64;
//    double m_mol = 0.03007012;
//    double cv = cp - R_GAS / m_mol;
//    double gam = cp / cv;
//
//    par_m->r = par_p->r;
//    par_m->p = par_p->p;
//    par_m->u = par_p->u;
//    par_m->v = -par_p->v;
    calc_prim_params(cp, m_mol, par_m);
}

void get_top_boundary(prim_t *par_m, prim_t *par_p) {
    double cp = 1014.16;
    double m_mol = 0.02869409;

    par_p->r = par_m->r;
    par_p->u = par_m->u;
    par_p->v = -par_m->v;
    par_p->p = par_m->p;
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_p->c[i_comp] = par_m->c[i_comp];
    }
//    double cp = 469.64;
//    double m_mol = 0.03007012;
//    double cv = cp - R_GAS / m_mol;
//    double gam = cp / cv;
//
//    par_p->r = par_m->r;
//    par_p->p = par_m->p;
//    par_p->u = par_m->u;
//    par_p->v = -par_m->v;
    calc_prim_params(cp, m_mol, par_p);
}
