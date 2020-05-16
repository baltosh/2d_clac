#include <stdio.h>
#include <math.h>
#include <string.h>
#include "global.h"

int main(int argc, char **argv) {
    double t = 0.0;
    int step = 0;

    init();
    calc_cnc();
    calc_cnc();
    save_vtk(step);
    while (t < TIME_MAX) {
        t += TAU;
        step++;
        zero_r();
        calc_cnc();
        calc_int();
        calc_flx();
        calc_new();
        calc_lim();
        if (step % SAVE_STEP == 0) save_vtk(step);
    }
    return 0;
}

void inverse_matr(double a_src[BASE_FN_COUNT][BASE_FN_COUNT],
                  double am[BASE_FN_COUNT][BASE_FN_COUNT], int N) {
    double a[BASE_FN_COUNT][BASE_FN_COUNT];
    for (int i = 0; i < BASE_FN_COUNT; i++) {
        for (int j = 0; j < BASE_FN_COUNT; j++) {
            a[i][j] = a_src[i][j];
        }
    }
    double detA = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[0][1] * a[1][2] * a[2][0]
                  - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];

    double m[3][3];
    m[0][0] = a[1][1] * a[2][2] - a[2][1] * a[1][2];
    m[0][1] = a[2][0] * a[1][2] - a[1][0] * a[2][2];
    m[0][2] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
    m[1][0] = a[2][1] * a[0][2] - a[0][1] * a[2][2];
    m[1][1] = a[0][0] * a[2][2] - a[2][0] * a[0][2];
    m[1][2] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
    m[2][0] = a[0][1] * a[1][2] - a[1][1] * a[0][2];
    m[2][1] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
    m[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];

    am[0][0] = m[0][0] / detA;
    am[0][1] = m[1][0] / detA;
    am[0][2] = m[2][0] / detA;
    am[1][0] = m[0][1] / detA;
    am[1][1] = m[1][1] / detA;
    am[1][2] = m[2][1] / detA;
    am[2][0] = m[0][2] / detA;
    am[2][1] = m[1][2] / detA;
    am[2][2] = m[2][2] / detA;
}

void mult_matr_vec(double matr[BASE_FN_COUNT][BASE_FN_COUNT], double vec[BASE_FN_COUNT],
                   double res[BASE_FN_COUNT]) {
    for (int i = 0; i < BASE_FN_COUNT; i++) {
        res[i] = 0.0;
        for (int j = 0; j < BASE_FN_COUNT; j++) {
            res[i] += matr[i][j] * vec[j];
        }
    }
}

void zero_r() {
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            memset(&r_int[i][j], 0, sizeof(data_t));
        }
    }
}

void init() {

    HX = (X_MAX - X_MIN) / CELLS_X_COUNT;
    HY = (Y_MAX - Y_MIN) / CELLS_Y_COUNT;

    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            point_t p[GP_CELL_COUNT];
            p[0].x = -1.0 / sqrt(3.0);
            p[0].y = -1.0 / sqrt(3.0);
            p[1].x = 1.0 / sqrt(3.0);
            p[1].y = -1.0 / sqrt(3.0);
            p[2].x = 1.0 / sqrt(3.0);
            p[2].y = 1.0 / sqrt(3.0);
            p[3].x = -1.0 / sqrt(3.0);
            p[3].y = 1.0 / sqrt(3.0);
            double xmin = X_MIN + i * HX;
            double ymin = Y_MIN + j * HY;
            double xmax = xmin + HX;
            double ymax = ymin + HY;
            for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                gp_cell[i][j][i_gp].x = 0.5 * (xmin + xmax) + p[i_gp].x * (xmax - xmin) * 0.5;
                gp_cell[i][j][i_gp].y = 0.5 * (ymin + ymax) + p[i_gp].y * (ymax - ymin) * 0.5;
                gw_cell[i][j][i_gp] = 1.0;
                gj_cell[i][j] = 0.25 * (xmax - xmin) * (ymax - ymin);
            }

            centerCell[i][j].x = 0.5 * (xmin + xmax);
            centerCell[i][j].y = 0.5 * (ymin + ymax);
        }
    }

    for (int i = 0; i <= CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            double xmin = X_MIN + i * HX;
            double ymin = Y_MIN + j * HY;
            double ymax = ymin + HY;
            double p[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                gp_edge_x[i][j][i_gp].x = xmin;
                gp_edge_x[i][j][i_gp].y = 0.5 * (ymin + ymax) + p[i_gp] * (ymax - ymin) * 0.5;
                gw_edge_x[i][j][i_gp] = 1.0;
                gj_edge_x[i][j] = 0.5 * (ymax - ymin);
            }
        }
    }

    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j <= CELLS_Y_COUNT; j++) {
            double xmin = X_MIN + i * HX;
            double ymin = Y_MIN + j * HY;
            double xmax = xmin + HX;
            double p[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                gp_edge_y[i][j][i_gp].x = 0.5 * (xmin + xmax) + p[i_gp] * (xmax - xmin) * 0.5;;
                gp_edge_y[i][j][i_gp].y = ymin;
                gw_edge_y[i][j][i_gp] = 1.0;
                gj_edge_y[i][j] = 0.5 * (xmax - xmin);
            }
        }
    }

    double matr[BASE_FN_COUNT][BASE_FN_COUNT];
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            for (int m = 0; m < BASE_FN_COUNT; m++) {
                for (int l = 0; l < BASE_FN_COUNT; l++) {
                    matr[m][l] = 0.0;
                    for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                        matr[m][l] += gw_cell[i][j][i_gp] *
                                      bf(m, i, j, gp_cell[i][j][i_gp].x, gp_cell[i][j][i_gp].y)
                                      * bf(l, i, j, gp_cell[i][j][i_gp].x, gp_cell[i][j][i_gp].y);
                    }
                    matr[m][l] *= gj_cell[i][j];
                }
            }
            inverse_matr(matr, matr_a[i][j], BASE_FN_COUNT);
        }
    }
    //todo change
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            double x = X_MIN + (i + 0.5) * HX;
            double y = Y_MIN + (j + 0.5) * HY;
            memset(&data[i][j], 0, sizeof(data_t));
            double r, p, u, v, cp_mix, m_mol_mix;
            double rc[COMPONENTS_COUNT], c[COMPONENTS_COUNT], cp[COMPONENTS_COUNT], m_mol[COMPONENTS_COUNT];
            if (y < Y_MAX / 2.0) {
                p = 2.152e+5;
                u = 0.0;
                v = 97.76;

                rc[0] = 1.2601;
                rc[1] = 0;
                rc[2] = 0;
                rc[3] = 0;

                cp[0] = 1014.16;
                cp[1] = 1014.16;
                cp[2] = 1014.16;
                cp[3] = 1014.16;

                m_mol[0] = 0.03007;
                m_mol[1] = 0.02805;
                m_mol[2] = 0.001;
                m_mol[3] = 0.657;
            } else {
                p = 2.152e+5;
                u = 0.0;
                v = -97.76;

                rc[0] = 1.2601;
                rc[1] = 0;
                rc[2] = 0;
                rc[3] = 0;

                cp[0] = 1014.16;
                cp[1] = 1014.16;
                cp[2] = 1014.16;
                cp[3] = 1014.16;

                m_mol[0] = 0.03007;
                m_mol[1] = 0.02805;
                m_mol[2] = 0.001;
                m_mol[3] = 0.657;
            }

            r = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                r += rc[i_comp];
            }

            cp_mix = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                cp_mix += rc[i_comp] / r * cp[i_comp];
            }

            r = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                r += c[i_comp] * cp[i_comp];
            }

            double tmp = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                tmp += c[i_comp] / m_mol[i_comp];
            }
            m_mol_mix = 1 / tmp;

            double cv = cp_mix - GAS_CONSTANT / m_mol_mix;
            double gam = cp_mix / cv;
            data[i][j].rc[0][0] = rc[0];
            data[i][j].rc[1][0] = rc[1];
            data[i][j].rc[2][0] = rc[2];
            data[i][j].rc[3][0] = rc[3];
            data[i][j].ru[0] = r * u;
            data[i][j].rv[0] = r * v;
            data[i][j].re[0] = p / (gam - 1.0) + r * (u * u + v * v) * 0.5;
        }
    }
}

void calc_cnc() {

    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            point_t center = centerCell[i][j];
            q_point[0][0].x = center.x - HX / 2;
            q_point[0][0].y = center.y;
            q_point[0][1].x = center.x;
            q_point[0][1].y = center.y + HY / 2;

            q_point[1][0].x = center.x;
            q_point[1][0].y = center.y;
            q_point[1][1].x = center.x + HX / 2;
            q_point[1][1].y = center.y + HY / 2;

            q_point[2][0].x = center.x - HX / 2;
            q_point[2][0].y = center.y - HY / 2;
            q_point[2][1].x = center.x;
            q_point[2][1].y = center.y;

            q_point[3][0].x = center.x;
            q_point[3][0].y = center.y - HY / 2;
            q_point[3][1].x = center.x + HX / 2;
            q_point[3][1].y = center.y;
//считаем интегральное среднее
            for (int i_quad = 0; i_quad < 4; i_quad++) {
                point_t p[GP_CELL_COUNT];
                p[0].x = -1.0 / sqrt(3.0);
                p[0].y = -1.0 / sqrt(3.0);
                p[1].x = 1.0 / sqrt(3.0);
                p[1].y = -1.0 / sqrt(3.0);
                p[2].x = 1.0 / sqrt(3.0);
                p[2].y = 1.0 / sqrt(3.0);
                p[3].x = -1.0 / sqrt(3.0);
                p[3].y = 1.0 / sqrt(3.0);
                double xmin = q_point[i_quad][0].x;
                double ymin = q_point[i_quad][0].y;
                double xmax = q_point[i_quad][1].x;
                double ymax = q_point[i_quad][1].y;
                for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                    q_gp_cell[i_quad][i_gp].x = 0.5 * (xmin + xmax) + p[i_gp].x * (xmax - xmin) * 0.5;
                    q_gp_cell[i_quad][i_gp].y = 0.5 * (ymin + ymax) + p[i_gp].y * (ymax - ymin) * 0.5;
                    q_gw_cell[i_quad][i_gp] = 1.0;
                    q_gj_cell[i_quad] = 0.25 * (xmax - xmin) * (ymax - ymin);
                }
            }
            for (int i_quad = 0; i_quad < 4; i_quad++) {
                double old_rc[COMPONENTS_COUNT], new_rc[COMPONENTS_COUNT];
                for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                    old_rc[i_com] = 0.0;
                    for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                        old_rc[i_com] +=
                                get_field_rc(i, j, q_gp_cell[i_quad][i_gp].x, q_gp_cell[i_quad][i_gp].y, i_com) *
                                q_gw_cell[i_quad][i_gp];
                    }
                    old_rc[i_com] *= q_gj_cell[i_quad] / (HX * HY / 4);
                }
//                memset(old_rc, 0, sizeof(double) * COMPONENTS_COUNT);
//                memset(new_rc, 0, sizeof(double) * COMPONENTS_COUNT);
                double t = 0.0;
                double dt = TAU / 10;

                while (t < TAU) {
                    t += dt;
                    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                        new_rc[i_com] = old_rc[i_com] + dt * reactionSpeeds[i_com](oldCnc);
                    }

                    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                        old_rc[i_com] = new_rc[i_com];
                    }
                }

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {


                }
            }

            for (int q = 0; q < 4; q++) {

                double oldCnc[COMPONENTS_COUNT], newCnc[COMPONENTS_COUNT];

                memset(oldCnc, 0, sizeof(double) * COMPONENTS_COUNT);
                memset(newCnc, 0, sizeof(double) * COMPONENTS_COUNT);

                prim_t par;
                cons_to_prim(i, j, centerQuad[q].x, centerQuad[q].y, &par);
                double *concentration = concentrations[i][j];

                double t = 0.0;
                double dt = TAU / 10;

                for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                    oldCnc[i_com] = par.r * concentration[i_com];
                }

                while (t < TAU) {
                    t += dt;
                    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                        newCnc[i_com] = oldCnc[i_com] + dt * reactionSpeeds[i_com](oldCnc);
                    }

                    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                        oldCnc[i_com] = newCnc[i_com];
                    }
                }

                for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                    newConcentration[i_com] += oldCnc[i_com] - TAU * (
                            (newCnc[i_com] - concentration[i_com] * par.r) * par.u / HX +
                            (newCnc[i_com] - concentration[i_com] * par.r) * par.v / HY);
                }
            }
            for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
                newConcentration[i_com] /= 4;
            }
        }
    }
};

void calc_int() {
    data_t fint;
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                point_t pt = gp_cell[i][j][i_gp];
                prim_t par;
                cons_to_prim(i, j, pt.x, pt.y, &par);
                double f[3 + COMPONENTS_COUNT], g[3 + COMPONENTS_COUNT];

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    f[i_comp] = par.r * par.c[i_comp] * par.u;
                }

                f[COMPONENTS_COUNT] = par.r * par.u * par.u + par.p;
                f[COMPONENTS_COUNT + 1] = par.r * par.u * par.v;
                f[COMPONENTS_COUNT + 2] = par.u * (par.r * par.e_tot + par.p);

                g[0] = par.r * par.v;
                g[COMPONENTS_COUNT] = par.r * par.v * par.u;
                g[COMPONENTS_COUNT + 1] = par.r * par.v * par.v + par.p;
                g[COMPONENTS_COUNT + 2] = par.v * (par.r * par.e_tot + par.p);

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        fint.fields[i_fld][k] += gw_cell[i][j][i_gp] * f[i_fld] * bf_dx(k, i, j, pt.x, pt.y);
                        fint.fields[i_fld][k] += gw_cell[i][j][i_gp] * g[i_fld] * bf_dy(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    fint.fields[i_fld][k] *= gj_cell[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i][j].fields[i_fld][k] += fint.fields[i_fld][k];
                }
            }


        }
    }
}

void calc_flx() {
    // X direction
    for (int i = 1; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_x[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i - 1, j, pt.x, pt.y, &par_m);
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    flx[i_comp] = 0.5 * ((par_p.r * par_p.c[i_comp] * par_p.u + par_m.r * par_m.c[i_comp] * par_m.u)
                                         - alpha * (par_p.r * par_p.c[i_comp] - par_m.r * par_m.c[i_comp]));
                }

                flx[COMPONENTS_COUNT] =
                        0.5 * ((par_p.r * par_p.u * par_p.u + par_p.p + par_m.r * par_m.u * par_m.u + par_m.p)
                               - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[COMPONENTS_COUNT + 1] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                                   - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[COMPONENTS_COUNT + 2] = 0.5 *
                                            (((par_p.r * par_p.e_tot + par_p.p) * par_p.u +
                                              (par_m.r * par_m.e_tot + par_m.p) * par_m.u)
                                             - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        int_m.fields[i_fld][k] += gw_edge_x[i][j][i_gp] * flx[i_fld] * bf(k, i - 1, j, pt.x, pt.y);
                        int_p.fields[i_fld][k] += gw_edge_x[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    int_m.fields[i_fld][k] *= gj_edge_x[i][j];
                    int_p.fields[i_fld][k] *= gj_edge_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i - 1][j].fields[i_fld][k] -= int_m.fields[i_fld][k];
                    r_int[i][j].fields[i_fld][k] += int_p.fields[i_fld][k];
                }
            }
        }
    }

    // Y direction
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 1; j < CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_y[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j - 1, pt.x, pt.y, &par_m);
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    flx[i_comp] = 0.5 * ((par_p.r * par_p.c[i_comp] * par_p.v + par_m.r * par_m.c[i_comp] * par_m.v)
                                         - alpha * (par_p.r * par_p.c[i_comp] - par_m.r * par_m.c[i_comp]));
                }

                flx[COMPONENTS_COUNT] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                               - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[COMPONENTS_COUNT + 1] =
                        0.5 * ((par_p.r * par_p.v * par_p.v + par_p.p + par_m.r * par_m.v * par_m.v + par_m.p)
                               - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[COMPONENTS_COUNT + 2] = 0.5 *
                                            (((par_p.r * par_p.e_tot + par_p.p) * par_p.v +
                                              (par_m.r * par_m.e_tot + par_m.p) * par_m.v)
                                             - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        int_m.fields[i_fld][k] += gw_edge_y[i][j][i_gp] * flx[i_fld] * bf(k, i, j - 1, pt.x, pt.y);
                        int_p.fields[i_fld][k] += gw_edge_y[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    int_m.fields[i_fld][k] *= gj_edge_y[i][j];
                    int_p.fields[i_fld][k] *= gj_edge_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i][j - 1].fields[i_fld][k] -= int_m.fields[i_fld][k];
                    r_int[i][j].fields[i_fld][k] += int_p.fields[i_fld][k];
                }
            }
        }
    }

    // Left
    for (int i = 0; i <= 0; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_x[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GAS_CONSTANT / m_mol;
                    double gam = cp / cv;

                    par_m.r = par_p.r;
                    par_m.u = -par_p.u;
                    par_m.v = par_p.v;
                    par_m.p = par_p.p;
                    par_m.e = par_m.p / (par_m.r * (gam - 1.0));
                    par_m.e_tot = par_m.e + (par_m.u * par_m.u + par_m.v * par_m.v) * 0.5;
                    par_m.cz = sqrt(gam * par_m.p / par_m.r);
                    par_m.t = par_m.e / cv;

                }
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    flx[i_comp] = 0.5 * ((par_p.r * par_p.c[i_comp] * par_p.u + par_m.r * par_m.c[i_comp] * par_m.u)
                                         - alpha * (par_p.r * par_p.c[i_comp] - par_m.r * par_m.c[i_comp]));
                }

                flx[COMPONENTS_COUNT] =
                        0.5 * ((par_p.r * par_p.u * par_p.u + par_p.p + par_m.r * par_m.u * par_m.u + par_m.p)
                               - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[COMPONENTS_COUNT + 1] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                                   - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[COMPONENTS_COUNT + 2] = 0.5 *
                                            (((par_p.r * par_p.e_tot + par_p.p) * par_p.u +
                                              (par_m.r * par_m.e_tot + par_m.p) * par_m.u)
                                             - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        int_p.fields[i_fld][k] += gw_edge_x[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    int_p.fields[i_fld][k] *= gj_edge_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i][j].fields[i_fld][k] += int_p.fields[i_fld][k];
                }
            }
        }
    }

    // Right
    for (int i = CELLS_X_COUNT; i <= CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_x[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i - 1, j, pt.x, pt.y, &par_m);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GAS_CONSTANT / m_mol;
                    double gam = cp / cv;

                    par_p.r = par_m.r;
                    par_p.u = -par_m.u;
                    par_p.v = par_m.v;
                    par_p.p = par_m.p;
                    par_p.e = par_p.p / (par_p.r * (gam - 1.0));
                    par_p.e_tot = par_p.e + (par_p.u * par_p.u + par_p.v * par_p.v) * 0.5;
                    par_p.cz = sqrt(gam * par_p.p / par_p.r);
                    par_p.t = par_p.e / cv;

                }
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    flx[i_comp] = 0.5 * ((par_p.r * par_p.c[i_comp] * par_p.u + par_m.r * par_m.c[i_comp] * par_m.u)
                                         - alpha * (par_p.r * par_p.c[i_comp] - par_m.r * par_m.c[i_comp]));
                }

                flx[COMPONENTS_COUNT] =
                        0.5 * ((par_p.r * par_p.u * par_p.u + par_p.p + par_m.r * par_m.u * par_m.u + par_m.p)
                               - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[COMPONENTS_COUNT + 1] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                                   - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[COMPONENTS_COUNT + 2] = 0.5 *
                                            (((par_p.r * par_p.e_tot + par_p.p) * par_p.u +
                                              (par_m.r * par_m.e_tot + par_m.p) * par_m.u)
                                             - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        int_m.fields[i_fld][k] += gw_edge_x[i][j][i_gp] * flx[i_fld] * bf(k, i - 1, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    int_m.fields[i_fld][k] *= gj_edge_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i - 1][j].fields[i_fld][k] -= int_m.fields[i_fld][k];
                }
            }
        }
    }

    // Bottom
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j <= 0; j++) {
            data_t int_m, int_p;

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_y[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GAS_CONSTANT / m_mol;
                    double gam = cp / cv;

                    par_m.r = 12.090;
                    par_m.p = 2.152e+5;
                    par_m.u = 0.0;
                    par_m.v = 97.76;
                    par_m.e = par_m.p / (par_m.r * (gam - 1.0));
                    par_m.e_tot = par_m.e + (par_m.u * par_m.u + par_m.v * par_m.v) * 0.5;
                    par_m.cz = sqrt(gam * par_m.p / par_m.r);
                    par_m.t = par_m.e / cv;

                }
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    flx[i_comp] = 0.5 * ((par_p.r * par_p.c[i_comp] * par_p.v + par_m.r * par_m.c[i_comp] * par_m.v)
                                         - alpha * (par_p.r * par_p.c[i_comp] - par_m.r * par_m.c[i_comp]));
                }

                flx[COMPONENTS_COUNT] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                               - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[COMPONENTS_COUNT + 1] =
                        0.5 * ((par_p.r * par_p.v * par_p.v + par_p.p + par_m.r * par_m.v * par_m.v + par_m.p)
                               - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[COMPONENTS_COUNT + 2] = 0.5 *
                                            (((par_p.r * par_p.e_tot + par_p.p) * par_p.v +
                                              (par_m.r * par_m.e_tot + par_m.p) * par_m.v)
                                             - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        int_p.fields[i_fld][k] += gw_edge_y[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    int_p.fields[i_fld][k] *= gj_edge_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i][j].fields[i_fld][k] += int_p.fields[i_fld][k];
                }
            }
        }
    }

    // Top
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = CELLS_Y_COUNT; j <= CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_y[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j - 1, pt.x, pt.y, &par_m);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GAS_CONSTANT / m_mol;
                    double gam = cp / cv;

                    par_p.r = par_m.r;
                    par_p.u = par_m.u;
                    par_p.v = -par_m.v;
                    par_p.p = par_m.p;
                    par_p.e = par_p.p / (par_p.r * (gam - 1.0));
                    par_p.e_tot = par_p.e + (par_p.u * par_p.u + par_p.v * par_p.v) * 0.5;
                    par_p.cz = sqrt(gam * par_p.p / par_p.r);
                    par_p.t = par_p.e / cv;

                }
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    flx[i_comp] = 0.5 * ((par_p.r * par_p.c[i_comp] * par_p.v + par_m.r * par_m.c[i_comp] * par_m.v)
                                         - alpha * (par_p.r * par_p.c[i_comp] - par_m.r * par_m.c[i_comp]));
                }

                flx[COMPONENTS_COUNT] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                               - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[COMPONENTS_COUNT + 1] =
                        0.5 * ((par_p.r * par_p.v * par_p.v + par_p.p + par_m.r * par_m.v * par_m.v + par_m.p)
                               - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[COMPONENTS_COUNT + 2] = 0.5 *
                                            (((par_p.r * par_p.e_tot + par_p.p) * par_p.v +
                                              (par_m.r * par_m.e_tot + par_m.p) * par_m.v)
                                             - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                    for (int k = 0; k < BASE_FN_COUNT; k++) {
                        int_m.fields[i_fld][k] += gw_edge_y[i][j][i_gp] * flx[i_fld] * bf(k, i, j - 1, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    int_m.fields[i_fld][k] *= gj_edge_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    r_int[i][j - 1].fields[i_fld][k] -= int_m.fields[i_fld][k];
                }
            }
        }
    }
}

void calc_new() {
    double tmp[BASE_FN_COUNT];
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                mult_matr_vec(matr_a[i][j], r_int[i][j].fields[i_fld], tmp);
                for (int k = 0; k < BASE_FN_COUNT; k++) {
                    data[i][j].fields[i_fld][k] += TAU * tmp[k];
                }
            }
        }
    }
}

double _minmod(double a, double b, double c) {
    if ((_SIGN_(a) == _SIGN_(b)) && (_SIGN_(b) == _SIGN_(c))) {
        return _SIGN_(a) * _MIN_(_MIN_(fabs(a), fabs(b)), fabs(c));
    } else {
        return 0.0;
    }
}

void calc_lim() {

    for (int i = 1; i < CELLS_X_COUNT - 1; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                double fm = data[i - 1][j].fields[i_fld][0];
                double fc = data[i][j].fields[i_fld][0];
                double fp = data[i + 1][j].fields[i_fld][0];
                data[i][j].fields[i_fld][1] = _minmod(data[i][j].fields[i_fld][1], LIMITER_ALPHA * (fc - fm),
                                                      LIMITER_ALPHA * (fp - fc));
            }
        }
    }
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 1; j < CELLS_Y_COUNT - 1; j++) {
            for (int i_fld = 0; i_fld < 3 + COMPONENTS_COUNT; i_fld++) {
                double fm = data[i][j - 1].fields[i_fld][0];
                double fc = data[i][j].fields[i_fld][0];
                double fp = data[i][j + 1].fields[i_fld][0];
                data[i][j].fields[i_fld][2] = _minmod(data[i][j].fields[i_fld][2], LIMITER_ALPHA * (fc - fm),
                                                      LIMITER_ALPHA * (fp - fc));
            }
        }
    }
}

double bf(int i_func, int i, int j, double x, double y) {
    switch (i_func) {
        case 0:
            return 1.0;
            break;
        case 1:
            return (x - centerCell[i][j].x) / HX;
            break;
        case 2:
            return (y - centerCell[i][j].y) / HY;
            break;
    }
}

double bf_dx(int i_func, int i, int j, double x, double y) {
    switch (i_func) {
        case 0:
            return 0.0;
            break;
        case 1:
            return 1.0 / HX;
            break;
        case 2:
            return 0.0;
            break;
    }
}

double bf_dy(int i_func, int i, int j, double x, double y) {
    switch (i_func) {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 1.0 / HY;
            break;
    }
}

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

    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        prim->c[i] = get_field_rc(i, j, x, y, i) / prim->r;
    }
}

double get_component_cp(int id) {
    double cps[COMPONENTS_COUNT] = {1014.16, 1014.16, 1014.16, 1014.16, 1014.16};
    return cps[id];
}

double get_component_M(int id) {
    double Ms[COMPONENTS_COUNT] = {0.02869409, 0.02869409, 0.02869409, 0.02869409, 0.02869409};
    return Ms[id];
}

double get_cell_cp(int i, int j) {
    double cp = 0;
    double gjw;

    for (int i_com = 0; i_com < COMPONENTS_COUNT; i_com++) {
        for (int igp = 0; igp < GP_CELL_COUNT; igp++) {
            gjw = gj_cell[i][j] * gw_cell[i][j][igp] / (HX * HY);
            cp += get_field_rc(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y, i_com) * get_component_cp(i_com) *
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
            M += get_field_rc(i, j, gp_cell[i][j][igp].x, gp_cell[i][j][igp].y, i_com) / get_component_M(i_com) *
                 gjw;
        }
    }

    return 1 / M;
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

//todo change
double reactionSpeed0(double concentration[]) {
    printf("Function 0\n");
    return 0;
};

double reactionSpeed1(double concentration[]) {
    printf("Function 1\n");
    return 0;
};

double reactionSpeed2(double concentration[]) {
    printf("Function 2\n");
    return 0;
};

double reactionSpeed3(double concentration[]) {
    printf("Function 3\n");
    return 0;
};

double reactionSpeed4(double concentration[]) {
    printf("Function 4\n");
    return 0;
};

void save_vtk(int num) {
    char fName[50];
    prim_t pr;
    sprintf(fName, "res_%010d.vtk", num);
    FILE *fp = fopen(fName, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "GASDIN data file\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", CELLS_X_COUNT + 1, CELLS_Y_COUNT + 1, 1);
    fprintf(fp, "POINTS %d float\n", (CELLS_X_COUNT + 1) * (CELLS_Y_COUNT + 1));
    for (int j = 0; j <= CELLS_Y_COUNT; j++) {
        for (int i = 0; i <= CELLS_X_COUNT; i++) {
            double x = X_MIN + i * HX;
            double y = Y_MIN + j * HY;
            fprintf(fp, "%f %f %f\n", x, y, 0.0);
        }
    }
    fprintf(fp, "CELL_DATA %d\n", CELLS_X_COUNT * CELLS_Y_COUNT);
    fprintf(fp, "SCALARS Density float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.r);
        }
    }

    fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.p);
        }
    }

    fprintf(fp, "SCALARS TotEnergy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e_tot);
        }
    }

    fprintf(fp, "SCALARS Energy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e);
        }
    }

    fprintf(fp, "VECTORS Velosity float\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f %f %f\n", pr.u, pr.v, 0.0);
        }
    }

    fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.t);
        }
    }

    fprintf(fp, "VECTORS Concentrations float\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            for (int k = 0; k < COMPONENTS_COUNT - 1; k++) {
                fprintf(fp, "%f ", concentrations[i][j][k]);
            }
            fprintf(fp, "%f\n", concentrations[i][j][COMPONENTS_COUNT - 1]);
        }
    }

    fclose(fp);
    printf("File '%s' saved...\n", fName);

}