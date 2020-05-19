//
// Created by nzvi on 19.05.2020.
//

#include <math.h>
#include "include/basis.h"
#include "include/parameters.h"

void cons_to_prim(int i, int j, double x, double y, prim_t *prim) {
    double cp = get_cell_cp(i, j);
    double m_mol = get_cell_M(i, j);

    double cv = cp - GAS_CONSTANT / m_mol;
    double gam = cp / cv;

    double ro = 0;
    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        ro += get_field_rc(i, j, x, y, i_comp);
    }
    double ru = get_field_ru(i, j, x, y);
    double rv = get_field_rv(i, j, x, y);
    double re = get_field_re(i, j, x, y);

    prim->r = ro;
    prim->u = ru / ro;
    prim->v = rv / ro;
    prim->e_tot = re / ro;
    prim->e = prim->e_tot - (prim->u * prim->u + prim->v * prim->v) * 0.5;
    prim->p = prim->r * prim->e * (gam - 1.0);
    prim->cz = sqrt(gam * prim->p / prim->r);
    prim->t = prim->e / cv;

    prim->cp = cp;
    prim->m_mol = m_mol;

    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        prim->c[i_comp] = get_field_rc(i, j, x, y, i_comp) / ro;
    }
}

double get_component_cp(int id) {
    //этан, этилен, водород, метан
    double cps[COMPONENTS_COUNT] = {1014.16/**469.64*/, 572.57, 14274.97, 925.55};
    return cps[id];
}

double get_component_M(int id) {
    //этан, этилен, водород, метан
    double Ms[COMPONENTS_COUNT] = {0.02869409/**0.03007012*/, 0.02805418, 0.00201594, 0.01604303};
    return Ms[id];
}

double get_cell_cp(int i, int j) {
    double cp = 0;
    double gjw;

    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (int igp = 0; igp < GP_CELL_COUNT; igp++) {
            gjw = gj_cell[i][j] * gw_cell[i][j][igp] / (HX * HY);
            cp += get_field_rc(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y, i_com) /
                  get_field_ro(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y) * get_component_cp(i_com) *
                  gjw;
        }
    }

    return cp;
}

double get_cell_M(int i, int j) {

    double M = 0;
    double gjw;

    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (int igp = 0; igp < GP_CELL_COUNT; igp++) {
            gjw = gj_cell[i][j] * gw_cell[i][j][igp] / (HX * HY);
            M += get_field_rc(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y, i_com) /
                 get_field_ro(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y) / get_component_M(i_com) *
                 gjw;
        }
    }

    return 1 / M;
}

double get_field_ro(int i, int j, double x, double y) {
    double result = 0.;
    int i_func;
    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
            result += data[i][j].rc[i_com][i_func] * bf(i_func, i, j, x, y);
        }
    }
    return result;
}

double get_field_ru(int i, int j, double x, double y) {
    double result = 0.;
    int i_func;

    for (i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
        result += data[i][j].ru[i_func] * bf(i_func, i, j, x, y);
    }
    return result;
}

double get_field_rv(int i, int j, double x, double y) {
    double result = 0.;
    int i_func;

    for (i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
        result += data[i][j].rv[i_func] * bf(i_func, i, j, x, y);
    }
    return result;
}

double get_field_re(int i, int j, double x, double y) {
    double result = 0.;
    int i_func;

    for (i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
        result += data[i][j].re[i_func] * bf(i_func, i, j, x, y);
    }
    return result;
}

double get_field_rc(int i, int j, double x, double y, int k) {
    double result = 0.;
    int i_func;

    for (i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
        result += data[i][j].rc[i_func][k] * bf(i_func, i, j, x, y);
    }
    return result;
}