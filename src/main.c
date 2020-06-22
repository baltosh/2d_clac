#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "include/global.h"
#include "include/basis.h"
#include "include/utils.h"
#include "include/chemistry.h"
#include "include/physics.h"
#include "include/parameters.h"

int step = 0;
int main(int argc, char **argv) {
    double t = 0.0;

    init();
    save_vtk(step);
    while (t < TIME_MAX) {
        t += TAU;
        step++;
        printf("%i\n", step);
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

            center_cell[i][j].x = 0.5 * (xmin + xmax);
            center_cell[i][j].y = 0.5 * (ymin + ymax);
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
            double rc[COMPONENTS_COUNT];

//            if (y < -0.005) {
//                r = 12.090;
//                p = 2.152e+5;
//                u = 0.0;
//                v = 97.76;
//            } else if (y > HY * sin(M_PI * x / 0.01)) {
//                r = 1.198;
//                p = 1.e+5;
//                u = 0.0;
//                v = 0.0;
//            } else {
//                r = 6.037;
//                p = 1.e+5;
//                u = 0.0;
//                v = 0.0;
//            }

            r = 12.293;
            p = 30.e+5;
            u = 0.0;
            v = 0.0;

            rc[0] = r;
            rc[1] = 0;
            rc[2] = 0;
            rc[3] = 0;

            r = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                r += rc[i_comp];
            }

            cp_mix = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                cp_mix += rc[i_comp] / r * get_component_cp(i_comp);
            }

            double tmp = 0.0;
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                tmp += rc[i_comp] / r / get_component_M(i_comp);
            }
            m_mol_mix = 1 / tmp;

            //p = r * GAS_CONSTANT * T / m_mol_mix;

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

            double quad_rc[4][COMPONENTS_COUNT];

            point_t center = center_cell[i][j];
            //первый квадрат x < xc, y > yc
            q_point[0][0].x = center.x - HX / 2;
            q_point[0][0].y = center.y;
            q_point[0][1].x = center.x;
            q_point[0][1].y = center.y + HY / 2;
            //второй квадрат x > xc, y > yc
            q_point[1][0].x = center.x;
            q_point[1][0].y = center.y;
            q_point[1][1].x = center.x + HX / 2;
            q_point[1][1].y = center.y + HY / 2;
            //третий квадрат x < xc, y < yc
            q_point[2][0].x = center.x - HX / 2;
            q_point[2][0].y = center.y - HY / 2;
            q_point[2][1].x = center.x;
            q_point[2][1].y = center.y;
            //четвертый квадрат x > xc, y > yc
            q_point[3][0].x = center.x;
            q_point[3][0].y = center.y - HY / 2;
            q_point[3][1].x = center.x + HX / 2;
            q_point[3][1].y = center.y;
            //определяем числа для квадратурных формул Гаусса
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
            //считаем интегральные средние для концентраций на начальный момент
            for (int i_quad = 0; i_quad < 4; i_quad++) {
                double old_rc[COMPONENTS_COUNT], new_rc[COMPONENTS_COUNT];
                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    old_rc[i_comp] = 0.0;

                    for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                        old_rc[i_comp] +=
                                get_field_rc(i, j, q_gp_cell[i_quad][i_gp].x, q_gp_cell[i_quad][i_gp].y, i_comp) *
                                q_gw_cell[i_quad][i_gp];
                    }
                    old_rc[i_comp] *= q_gj_cell[i_quad] / (HX * HY / 4);
                }
                double temperature = 0.0;
                for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                    temperature += get_field_T(i, j, q_gp_cell[i_quad][i_gp].x, q_gp_cell[i_quad][i_gp].y) * q_gw_cell[i_quad][i_gp];
                }
                temperature *= q_gj_cell[i_quad] / (HX * HY / 4);
                //решаем систему диффуров относительно интегральных средних для концентраций
                double t = 0.0;
                double dt = TAU;// / 10;

                while (t < TAU) {
                    t += dt;
                    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                        new_rc[i_comp] = old_rc[i_comp] + dt * reactionSpeeds[i_comp](old_rc, temperature);
                    }

                    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                        old_rc[i_comp] = new_rc[i_comp];
                    }
                }
                //интегральные средние для концентраций в новый момент времени в переменной old_rc
                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    quad_rc[i_quad][i_comp] = old_rc[i_comp];
                }
            }
            //находим координаты разложение функций концентраций на базисные функции
            double rc_bf[BASE_FN_COUNT], b[BASE_FN_COUNT];
            for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                for (int i_bf = 0; i_bf < BASE_FN_COUNT; i_bf++) {
                    b[i_bf] = 0.0;

                    for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                        point_t gauss_point = gp_cell[i][j][i_gp];
                        double rc_cap;
                        if (gauss_point.x < center.x) {
                            if (gauss_point.y < center.y) {
                                rc_cap = quad_rc[2][i_comp];
                            } else {
                                rc_cap = quad_rc[0][i_comp];
                            }
                        } else {
                            if (gauss_point.y < gauss_point.y) {
                                rc_cap = quad_rc[3][i_comp];
                            } else {
                                rc_cap = quad_rc[1][i_comp];
                            }
                        }

                        b[i_bf] += rc_cap * bf(i_bf, i, j, gauss_point.x, gauss_point.y) * gw_cell[i][j][i_gp];
                    }
                    b[i_bf] *= gj_cell[i][j];
                }

                mult_matr_vec(matr_a[i][j], b, rc_bf);
                for (int i_bf = 0; i_bf < BASE_FN_COUNT; i_bf++) {
                    data[i][j].rc[i_comp][i_bf] = rc_bf[i_bf];
                }
            }
        }
    }
}

void calc_int() {
    data_t fint;
    memset(&fint, 0, sizeof(data_t));
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            for (int i_gp = 0; i_gp < GP_CELL_COUNT; i_gp++) {
                point_t pt = gp_cell[i][j][i_gp];
                prim_t par;
                memset(&par, 0, sizeof(prim_t));
                cons_to_prim(i, j, pt.x, pt.y, &par);
                double f[3 + COMPONENTS_COUNT], g[3 + COMPONENTS_COUNT];

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    f[i_comp] = par.r * par.c[i_comp] * par.u;
                }

                f[COMPONENTS_COUNT] = par.r * par.u * par.u + par.p;
                f[COMPONENTS_COUNT + 1] = par.r * par.u * par.v;
                f[COMPONENTS_COUNT + 2] = par.u * (par.r * par.e_tot + par.p);

                for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
                    g[i_comp] = par.r * par.c[i_comp] * par.v;
                }
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
            memset(&int_m, 0, sizeof(data_t));
            memset(&int_p, 0, sizeof(data_t));

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_x[i][j][i_gp];
                prim_t par_m, par_p;
                memset(&par_m, 0, sizeof(prim_t));
                memset(&par_p, 0, sizeof(prim_t));
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i - 1, j, pt.x, pt.y, &par_m);
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                calc_x_flx(&par_m, &par_p, flx);

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
            memset(&int_m, 0, sizeof(data_t));
            memset(&int_p, 0, sizeof(data_t));

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_y[i][j][i_gp];
                prim_t par_m, par_p;
                memset(&par_m, 0, sizeof(prim_t));
                memset(&par_p, 0, sizeof(prim_t));
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j - 1, pt.x, pt.y, &par_m);
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                calc_y_flx(&par_m, &par_p, flx);

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

    // Left Чистый этан втекает в трубу слева
    for (int i = 0; i <= 0; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;
            memset(&int_m, 0, sizeof(data_t));
            memset(&int_p, 0, sizeof(data_t));

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_x[i][j][i_gp];
                prim_t par_m, par_p;
                memset(&par_m, 0, sizeof(prim_t));
                memset(&par_p, 0, sizeof(prim_t));
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                get_left_boundary(&par_m, &par_p);

                calc_x_flx(&par_m, &par_p, flx);

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

    // Right Правый конец трубы является стоком
    for (int i = CELLS_X_COUNT; i <= CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            data_t int_m, int_p;
            memset(&int_m, 0, sizeof(data_t));
            memset(&int_p, 0, sizeof(data_t));

            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_x[i][j][i_gp];
                prim_t par_m, par_p;
                memset(&par_m, 0, sizeof(prim_t));
                memset(&par_p, 0, sizeof(prim_t));
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i - 1, j, pt.x, pt.y, &par_m);
                get_right_boundary(&par_m, &par_p);

                calc_x_flx(&par_m, &par_p, flx);


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
//Поток через верхнюю и нижнюю границы не учитывается, так как на них опредлено условие непротекания
    // Bottom
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j <= 0; j++) {
            data_t int_m, int_p;
            memset(&int_m, 0, sizeof(data_t));
            memset(&int_p, 0, sizeof(data_t));
            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_y[i][j][i_gp];
                prim_t par_m, par_p;
                memset(&par_m, 0, sizeof(prim_t));
                memset(&par_p, 0, sizeof(prim_t));
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                get_bottom_boundary(&par_m, &par_p);

                calc_y_flx(&par_m, &par_p, flx);

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
            memset(&int_m, 0, sizeof(data_t));
            memset(&int_p, 0, sizeof(data_t));
            for (int i_gp = 0; i_gp < GP_EDGE_COUNT; i_gp++) {
                point_t pt = gp_edge_y[i][j][i_gp];
                prim_t par_m, par_p;
                memset(&par_m, 0, sizeof(prim_t));
                memset(&par_p, 0, sizeof(prim_t));
                double flx[3 + COMPONENTS_COUNT];
                cons_to_prim(i, j - 1, pt.x, pt.y, &par_m);
                get_top_boundary(&par_m, &par_p);

                calc_y_flx(&par_m, &par_p, flx);

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
//            prim_t test_prim;
//            test(i, j, &test_prim);
//            if (test_prim.r != test_prim.r) {
//                printf("Density in %i %i Is NAN!!!11\n", i, j);
//            }

        }
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

void calc_y_flx(prim_t *par_m, prim_t *par_p, double *flx) {
    double alpha = _MAX_(par_m->cz + sqrt(par_m->u * par_m->u + par_m->v * par_m->v),
                         par_p->cz + sqrt(par_p->u * par_p->u + par_p->v * par_p->v));

    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        flx[i_comp] = 0.5 * ((par_p->c[i_comp] * par_p->r * par_p->v + par_m->c[i_comp] * par_m->r * par_m->v)
                             - alpha * (par_p->c[i_comp] * par_p->r - par_m->c[i_comp] * par_m->r));
    }

    flx[COMPONENTS_COUNT] = 0.5 * ((par_p->r * par_p->u * par_p->v + par_m->r * par_m->u * par_m->v)
                                   - alpha * (par_p->r * par_p->u - par_m->r * par_m->u));
    flx[COMPONENTS_COUNT + 1] =
            0.5 * ((par_p->r * par_p->v * par_p->v + par_p->p + par_m->r * par_m->v * par_m->v + par_m->p)
                   - alpha * (par_p->r * par_p->v - par_m->r * par_m->v));
    flx[COMPONENTS_COUNT + 2] = 0.5 *
                                (((par_p->r * par_p->e_tot + par_p->p) * par_p->v +
                                  (par_m->r * par_m->e_tot + par_m->p) * par_m->v)
                                 - alpha * (par_p->r * par_p->e_tot - par_m->r * par_m->e_tot));
}

void calc_x_flx(prim_t *par_m, prim_t *par_p, double *flx) {
    double alpha = _MAX_(par_m->cz + sqrt(par_m->u * par_m->u + par_m->v * par_m->v),
                         par_p->cz + sqrt(par_p->u * par_p->u + par_p->v * par_p->v));

    for (int i_comp = 0; i_comp < COMPONENTS_COUNT; i_comp++) {
        flx[i_comp] = 0.5 * ((par_p->c[i_comp] * par_p->r * par_p->u + par_m->c[i_comp] * par_m->r * par_m->u)
                             - alpha * (par_p->c[i_comp] * par_p->r - par_m->c[i_comp] * par_m->r));
    }

    flx[COMPONENTS_COUNT] =
            0.5 * ((par_p->r * par_p->u * par_p->u + par_p->p + par_m->r * par_m->u * par_m->u + par_m->p)
                   - alpha * (par_p->r * par_p->u - par_m->r * par_m->u));
    flx[COMPONENTS_COUNT + 1] = 0.5 * ((par_p->r * par_p->u * par_p->v + par_m->r * par_m->u * par_m->v)
                                       - alpha * (par_p->r * par_p->v - par_m->r * par_m->v));
    flx[COMPONENTS_COUNT + 2] = 0.5 *
                                (((par_p->r * par_p->e_tot + par_p->p) * par_p->u +
                                  (par_m->r * par_m->e_tot + par_m->p) * par_m->u)
                                 - alpha * (par_p->r * par_p->e_tot - par_m->r * par_m->e_tot));
}

void test(int i, int j, prim_t *pr) {
    memset(pr, 0, sizeof(prim_t));
    cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, pr);
}