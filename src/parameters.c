//
// Created by nzvi on 19.05.2020.
//

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

double get_cell_cp(int i, int j) {
    double cp = 0;

    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (int igp = 0; igp < GP_CELL_COUNT; igp++) {
            cp += get_field_rc(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y, i_com) /
                  get_field_ro(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y) * get_component_cp(i_com) *
                  gw_cell[i][j][igp];
        }
    }

    cp *= gj_cell[i][j] / (HX * HY);
    return cp;
}

double get_cell_M(int i, int j) {
    double M = 0;

    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (int igp = 0; igp < GP_CELL_COUNT; igp++) {
            M += get_field_rc(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y, i_com) /
                 get_field_ro(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y) / get_component_M(i_com) *
                 gw_cell[i][j][igp];
        }
    }

    M *= gj_cell[i][j] / (HX * HY);
    return 1 / M;
}

double get_field_T(int i, int j, double x, double y) {
    double cp = get_cell_cp(i, j);
    double m_mol = get_cell_M(i, j);
    double ro = get_field_ro(i, j, x, y);
    double re = get_field_re(i, j, x, y);
    double u = get_field_ru(i, j, x, y) / ro;
    double v = get_field_rv(i, j, x, y) / ro;

    double cv = cp - GAS_CONSTANT / m_mol;
    double e_tot = re / ro;
    double e = e_tot - (u * u + v * v) * 0.5;
    double t = e / cv;
    return t;
}

double get_field_ro(int i, int j, double x, double y) {
    double result = 0.;

    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (int i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
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

double get_field_rc(int i, int j, double x, double y, int i_com) {
    double result = 0.;
    int i_func;

    for (i_func = 0; i_func < BASE_FN_COUNT; i_func++) {
        result += data[i][j].rc[i_com][i_func] * bf(i_func, i, j, x, y);
    }
    return result;
}