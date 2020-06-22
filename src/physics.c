//
// Created by nzvi on 19.05.2020.
//

#include "include/physics.h"


void calc_prim_params(prim_t *par) {
    double cv = par->cp - GAS_CONSTANT / par->m_mol;
    double gam = par->cp / cv;

    par->e = par->p / (par->r * (gam - 1.0));
    par->e_tot = par->e + (par->u * par->u + par->v * par->v) * 0.5;
    par->cz = sqrt(gam * par->p / par->r);
    par->t = par->e / cv;
}

void get_left_boundary(prim_t *par_m, prim_t *par_p) {
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_m->c[i_comp] = 0.0;
    }
    par_m->c[0] = 1.0;

    par_m->cp = get_component_cp(0);
    par_m->m_mol = get_component_M(0);

    par_m->r = 12.293;
    par_m->u = 10;
    par_m->v = 0.0;
    par_m->p = 30.e+5;//par_m->r * GAS_CONSTANT * par_p->t / par_m->m_mol;
    calc_prim_params(par_m);
}

void get_right_boundary(prim_t *par_m, prim_t *par_p) {
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_p->c[i_comp] = par_m->c[i_comp];
    }
    par_p->cp = par_m->cp;
    par_p->m_mol = par_m->m_mol;

    par_p->r = par_m->r;
    par_p->u = par_m->u;
    par_p->v = par_m->v;
    par_p->p = par_m->p;
    calc_prim_params(par_p);
}

void get_bottom_boundary(prim_t *par_m, prim_t *par_p) {

//    double cp = 1014.16;
//    double m_mol = 0.02869409;
//
//    par_m->r = 12.090;
//    par_m->u = 0.0;
//    par_m->v = 97.76;
//    par_m->p = 2.152e+5;
//    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
//        par_m->c[i_comp] = 0.0;
//    }
//    par_m->c[0] = 1.0;

    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_m->c[i_comp] = par_p->c[i_comp];
    }
    par_m->cp = par_p->cp;
    par_m->m_mol = par_p->m_mol;

    par_m->r = par_p->r;
    par_m->p = par_p->p;
    par_m->u = par_p->u;
    par_m->v = -par_p->v;
    calc_prim_params(par_m);
}

void get_top_boundary(prim_t *par_m, prim_t *par_p) {
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        par_p->c[i_comp] = par_m->c[i_comp];
    }
    par_p->cp = par_m->cp;
    par_p->m_mol = par_m->m_mol;

    par_p->r = par_m->r;
    par_p->p = par_m->p;
    par_p->u = par_m->u;
    par_p->v = -par_m->v;
    calc_prim_params(par_p);
}
