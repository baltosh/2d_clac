#include <stdio.h>
#include <math.h>
#include <string.h>

#define M_PI        3.14159265358979323846
#define GasConstant 8.314472

#define NumberBasisFunctions    3
#define NumberCellsX            80
#define NumberCellsY            800
#define NumberComponents        5

#define NumberGaussPointsEdge 2
#define NumberGaussPointsCell 4

#define XMin  0.0
#define XMax  0.04
#define YMin -0.1
#define YMax  0.3

#define Tau     5.0e-8
#define TimeMax 2.5e-3

#define SaveStep 100

#define LimiterAlpha 2.0

#define _SIGN_(X) (fabs(X)/(X))
#define _MIN_(X, Y) ((X)<(Y) ? (X) : (Y))
#define _MAX_(X, Y) ((X)>(Y) ? (X) : (Y))

typedef struct {
    union {
        double fields[4][NumberBasisFunctions];
        struct {
            double ro[NumberBasisFunctions];
            double ru[NumberBasisFunctions];
            double rv[NumberBasisFunctions];
            double re[NumberBasisFunctions];
        };
    };
} data_t;

typedef struct {
    double r, p, u, v, e, e_tot, t, cz;
} prim_t;

typedef struct {
    double x, y;
} point_t;

typedef struct {
    int id[NumberComponents];
    double concentration[NumberComponents];
} components_t;

double HX, HY;

data_t data[NumberCellsX][NumberCellsY], r_int[NumberCellsX][NumberCellsY];
double concentrations[NumberCellsX][NumberCellsY][NumberComponents];

point_t gaussPointEdgeX[NumberCellsX + 1][NumberCellsY][NumberGaussPointsEdge];
point_t gaussPointEdgeY[NumberCellsX][NumberCellsY + 1][NumberGaussPointsEdge];
double gaussWeightEdgeX[NumberCellsX + 1][NumberCellsY][NumberGaussPointsEdge];
double gaussWeightEdgeY[NumberCellsX][NumberCellsY + 1][NumberGaussPointsEdge];
double gj_edg_x[NumberCellsX + 1][NumberCellsY];
double gj_edg_y[NumberCellsX][NumberCellsY + 1];

point_t gaussPointCell[NumberCellsX][NumberCellsY][NumberGaussPointsCell];
double gaussWeightCell[NumberCellsX][NumberCellsY][NumberGaussPointsCell];
double gj_cell[NumberCellsX][NumberCellsY];
point_t centerCell[NumberCellsX][NumberCellsY];


double matr_a[NumberCellsX][NumberCellsY][NumberBasisFunctions][NumberBasisFunctions];

void init();

void calc_cnc();

void calc_flx();

void calc_int();

void calc_new();

void calc_lim();

void zero_r();

void save_vtk(int step);

double bf(int i_func, int i, int j, double x, double y);

double bf_dx(int i_func, int i, int j, double x, double y);

double bf_dy(int i_func, int i, int j, double x, double y);

double get_fld(int i_fld, int i, int j, double x, double y);

void cons_to_prim(int i, int j, double x, double y, prim_t *prim);

double get_component_cp(int id);

double get_component_M(int id);

double get_cell_cp(double concentration[NumberComponents]);

double get_cell_M(double concentration[NumberComponents]);

double reactionSpeed0(double[]);

double reactionSpeed1(double[]);

double reactionSpeed2(double[]);

double reactionSpeed3(double[]);

double reactionSpeed4(double[]);

double (*reactionSpeeds[NumberComponents])(double[]) = {reactionSpeed0,
                                                        reactionSpeed1,
                                                        reactionSpeed2,
                                                        reactionSpeed3,
                                                        reactionSpeed4};


int main(int argc, char **argv) {
    double t = 0.0;
    int step = 0;

    init();
    calc_cnc();
    calc_cnc();
    save_vtk(step);
    while (t < TimeMax) {
        t += Tau;
        step++;
        zero_r();
        calc_cnc();
        calc_int();
        calc_flx();
        calc_new();
        calc_lim();
        if (step % SaveStep == 0) save_vtk(step);
    }
    return 0;
}

void inverse_matr(double a_src[NumberBasisFunctions][NumberBasisFunctions],
                  double am[NumberBasisFunctions][NumberBasisFunctions], int N) {
    double a[NumberBasisFunctions][NumberBasisFunctions];
    for (int i = 0; i < NumberBasisFunctions; i++) {
        for (int j = 0; j < NumberBasisFunctions; j++) {
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

void mult_matr_vec(double matr[NumberBasisFunctions][NumberBasisFunctions], double vec[NumberBasisFunctions],
                   double res[NumberBasisFunctions]) {
    for (int i = 0; i < NumberBasisFunctions; i++) {
        res[i] = 0.0;
        for (int j = 0; j < NumberBasisFunctions; j++) {
            res[i] += matr[i][j] * vec[j];
        }
    }
}

void zero_r() {
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            memset(&r_int[i][j], 0, sizeof(data_t));
        }
    }
}

void init() {
    HX = (XMax - XMin) / NumberCellsX;
    HY = (YMax - YMin) / NumberCellsY;


    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            point_t p[NumberGaussPointsCell];
            p[0].x = -1.0 / sqrt(3.0);
            p[0].y = -1.0 / sqrt(3.0);
            p[1].x = 1.0 / sqrt(3.0);
            p[1].y = -1.0 / sqrt(3.0);
            p[2].x = 1.0 / sqrt(3.0);
            p[2].y = 1.0 / sqrt(3.0);
            p[3].x = -1.0 / sqrt(3.0);
            p[3].y = 1.0 / sqrt(3.0);
            double xmin = XMin + i * HX;
            double ymin = YMin + j * HY;
            double xmax = xmin + HX;
            double ymax = ymin + HY;
            for (int i_gp = 0; i_gp < NumberGaussPointsCell; i_gp++) {
                gaussPointCell[i][j][i_gp].x = 0.5 * (xmin + xmax) + p[i_gp].x * (xmax - xmin) * 0.5;
                gaussPointCell[i][j][i_gp].y = 0.5 * (ymin + ymax) + p[i_gp].y * (ymax - ymin) * 0.5;
                gaussWeightCell[i][j][i_gp] = 1.0;
                gj_cell[i][j] = 0.25 * (xmax - xmin) * (ymax - ymin);
            }

            centerCell[i][j].x = 0.5 * (xmin + xmax);
            centerCell[i][j].y = 0.5 * (ymin + ymax);
        }
    }

    for (int i = 0; i <= NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            double xmin = XMin + i * HX;
            double ymin = YMin + j * HY;
            double ymax = ymin + HY;
            double p[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                gaussPointEdgeX[i][j][i_gp].x = xmin;
                gaussPointEdgeX[i][j][i_gp].y = 0.5 * (ymin + ymax) + p[i_gp] * (ymax - ymin) * 0.5;
                gaussWeightEdgeX[i][j][i_gp] = 1.0;
                gj_edg_x[i][j] = 0.5 * (ymax - ymin);
            }
        }
    }

    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j <= NumberCellsY; j++) {
            double xmin = XMin + i * HX;
            double ymin = YMin + j * HY;
            double xmax = xmin + HX;
            double p[] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                gaussPointEdgeY[i][j][i_gp].x = 0.5 * (xmin + xmax) + p[i_gp] * (xmax - xmin) * 0.5;;
                gaussPointEdgeY[i][j][i_gp].y = ymin;
                gaussWeightEdgeY[i][j][i_gp] = 1.0;
                gj_edg_y[i][j] = 0.5 * (xmax - xmin);
            }
        }
    }

    double matr[NumberBasisFunctions][NumberBasisFunctions];
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            for (int m = 0; m < NumberBasisFunctions; m++) {
                for (int l = 0; l < NumberBasisFunctions; l++) {
                    matr[m][l] = 0.0;
                    for (int i_gp = 0; i_gp < NumberGaussPointsCell; i_gp++) {
                        matr[m][l] += gaussWeightCell[i][j][i_gp] *
                                      bf(m, i, j, gaussPointCell[i][j][i_gp].x, gaussPointCell[i][j][i_gp].y)
                                      * bf(l, i, j, gaussPointCell[i][j][i_gp].x, gaussPointCell[i][j][i_gp].y);
                    }
                    matr[m][l] *= gj_cell[i][j];
                }
            }
            inverse_matr(matr, matr_a[i][j], NumberBasisFunctions);
        }
    }
    //todo change
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            double x = XMin + (i + 0.5) * HX;
            double y = YMin + (j + 0.5) * HY;
            memset(&data[i][j], 0, sizeof(data_t));
            double r, p, u, v, cp, m_mol;
            if (y < -0.005) {
                r = 12.090;
                p = 2.152e+5;
                u = 0.0;
                v = 97.76;
                cp = 1014.16;
                m_mol = 0.02869409;
            } else if (y > HY * sin(M_PI * x / 0.01)) {
                r = 1.198;
                p = 1.e+5;
                u = 0.0;
                v = 0.0;
                cp = 1014.16;
                m_mol = 0.02869409;
            } else {
                r = 6.037;
                p = 1.e+5;
                u = 0.0;
                v = 0.0;
                cp = 1014.16;
                m_mol = 0.02869409;
            }
            double cv = cp - GasConstant / m_mol;
            double gam = cp / cv;
            data[i][j].ro[0] = r;
            data[i][j].ru[0] = r * u;
            data[i][j].rv[0] = r * v;
            data[i][j].re[0] = p / (gam - 1.0) + r * (u * u + v * v) * 0.5;
        }
    }

    //concentrations initialization
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            memset(&concentrations[i][j], 0, sizeof(double) * NumberComponents);
            concentrations[i][j][0] = 0.15;
            concentrations[i][j][1] = 0.2;
            concentrations[i][j][2] = 0.3;
            concentrations[i][j][3] = 0.2;
            concentrations[i][j][4] = 0.15;
        }
    }
}

void calc_cnc() {

    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            point_t center = centerCell[i][j];

            point_t centerQuad[4];
            memset(centerQuad, 0, sizeof(point_t) * 4);

            centerQuad[0].x = center.x + HX / 4;
            centerQuad[0].y = center.y + HY / 4;
            centerQuad[1].x = center.x + HX / 4;
            centerQuad[1].y = center.y - HY / 4;
            centerQuad[2].x = center.x - HX / 4;
            centerQuad[2].y = center.y - HY / 4;
            centerQuad[3].x = center.x - HX / 4;
            centerQuad[3].y = center.y + HY / 4;

            double newConcentration[NumberComponents] = {0.0, 0.0, 0.0, 0.0, 0.0};

            for (int q = 0; q < 4; q++) {

                double oldCnc[NumberComponents], newCnc[NumberComponents];

                memset(oldCnc, 0, sizeof(double) * NumberComponents);
                memset(newCnc, 0, sizeof(double) * NumberComponents);

                prim_t par;
                cons_to_prim(i, j, centerQuad[q].x, centerQuad[q].y, &par);
                double *concentration = concentrations[i][j];

                double t = 0.0;
                double dt = Tau / 10;

                for (int i_com = 0; i_com < NumberComponents; i_com++) {
                    oldCnc[i_com] = par.r * concentration[i_com];
                }

                while (t < Tau) {
                    t += dt;
                    for (int i_com = 0; i_com < NumberComponents; i_com++) {
                        newCnc[i_com] = oldCnc[i_com] + dt * reactionSpeeds[i_com](oldCnc);
                    }

                    for (int i_com = 0; i_com < NumberComponents; i_com++) {
                        oldCnc[i_com] = newCnc[i_com];
                    }
                }

                for (int i_com = 0; i_com < NumberComponents; i_com++) {
                    newConcentration[i_com] += oldCnc[i_com] - Tau * (
                            (newCnc[i_com] - concentration[i_com] * par.r) * par.u / HX +
                            (newCnc[i_com] - concentration[i_com] * par.r) * par.v / HY);
                }
            }
            for (int i_com = 0; i_com < NumberComponents; i_com++) {
                newConcentration[i_com] /= 4;
            }
        }
    }
};

void calc_int() {
    double fint[4][NumberBasisFunctions];
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(fint[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsCell; i_gp++) {
                point_t pt = gaussPointCell[i][j][i_gp];
                prim_t par;
                cons_to_prim(i, j, pt.x, pt.y, &par);
                double f[4], g[4];
                f[0] = par.r * par.u;
                f[1] = f[0] * par.u + par.p;
                f[2] = f[0] * par.v;
                f[3] = par.u * (par.r * par.e_tot + par.p);

                g[0] = par.r * par.v;
                g[1] = g[0] * par.u;
                g[2] = g[0] * par.v + par.p;
                g[3] = par.v * (par.r * par.e_tot + par.p);

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        fint[i_fld][k] += gaussWeightCell[i][j][i_gp] * f[i_fld] * bf_dx(k, i, j, pt.x, pt.y);
                        fint[i_fld][k] += gaussWeightCell[i][j][i_gp] * g[i_fld] * bf_dy(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    fint[i_fld][k] *= gj_cell[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i][j].fields[i_fld][k] += fint[i_fld][k];
                }
            }


        }
    }
}

void calc_flx() {
    // X direction
    for (int i = 1; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            double int_m[4][NumberBasisFunctions], int_p[4][NumberBasisFunctions];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double) * NumberBasisFunctions);
                memset(int_p[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                point_t pt = gaussPointEdgeX[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i - 1, j, pt.x, pt.y, &par_m);
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));
                flx[0] = 0.5 * ((par_p.r * par_p.u + par_m.r * par_m.u) - alpha * (par_p.r - par_m.r));
                flx[1] = 0.5 * ((par_p.r * par_p.u * par_p.u + par_p.p + par_m.r * par_m.u * par_m.u + par_m.p)
                                - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[2] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[3] = 0.5 *
                         (((par_p.r * par_p.e_tot + par_p.p) * par_p.u + (par_m.r * par_m.e_tot + par_m.p) * par_m.u)
                          - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        int_m[i_fld][k] += gaussWeightEdgeX[i][j][i_gp] * flx[i_fld] * bf(k, i - 1, j, pt.x, pt.y);
                        int_p[i_fld][k] += gaussWeightEdgeX[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    int_m[i_fld][k] *= gj_edg_x[i][j];
                    int_p[i_fld][k] *= gj_edg_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i - 1][j].fields[i_fld][k] -= int_m[i_fld][k];
                    r_int[i][j].fields[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }

    // Y direction
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 1; j < NumberCellsY; j++) {
            double int_m[4][NumberBasisFunctions], int_p[4][NumberBasisFunctions];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double) * NumberBasisFunctions);
                memset(int_p[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                point_t pt = gaussPointEdgeY[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j - 1, pt.x, pt.y, &par_m);
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                double alpha = _MAX_(par_m.cz + sqrt(par_m.u * par_m.u + par_m.v * par_m.v),
                                     par_p.cz + sqrt(par_p.u * par_p.u + par_p.v * par_p.v));
                flx[0] = 0.5 * ((par_p.r * par_p.v + par_m.r * par_m.v) - alpha * (par_p.r - par_m.r));
                flx[1] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[2] = 0.5 * ((par_p.r * par_p.v * par_p.v + par_p.p + par_m.r * par_m.v * par_m.v + par_m.p)
                                - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[3] = 0.5 *
                         (((par_p.r * par_p.e_tot + par_p.p) * par_p.v + (par_m.r * par_m.e_tot + par_m.p) * par_m.v)
                          - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        int_m[i_fld][k] += gaussWeightEdgeY[i][j][i_gp] * flx[i_fld] * bf(k, i, j - 1, pt.x, pt.y);
                        int_p[i_fld][k] += gaussWeightEdgeY[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    int_m[i_fld][k] *= gj_edg_y[i][j];
                    int_p[i_fld][k] *= gj_edg_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i][j - 1].fields[i_fld][k] -= int_m[i_fld][k];
                    r_int[i][j].fields[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }

    // Left
    for (int i = 0; i <= 0; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            double int_m[4][NumberBasisFunctions], int_p[4][NumberBasisFunctions];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_p[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                point_t pt = gaussPointEdgeX[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GasConstant / m_mol;
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
                flx[0] = 0.5 * ((par_p.r * par_p.u + par_m.r * par_m.u) - alpha * (par_p.r - par_m.r));
                flx[1] = 0.5 * ((par_p.r * par_p.u * par_p.u + par_p.p + par_m.r * par_m.u * par_m.u + par_m.p)
                                - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[2] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[3] = 0.5 *
                         (((par_p.r * par_p.e_tot + par_p.p) * par_p.u + (par_m.r * par_m.e_tot + par_m.p) * par_m.u)
                          - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        int_p[i_fld][k] += gaussWeightEdgeX[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    int_p[i_fld][k] *= gj_edg_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i][j].fields[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }

    // Right
    for (int i = NumberCellsX; i <= NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            double int_m[4][NumberBasisFunctions], int_p[4][NumberBasisFunctions];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                point_t pt = gaussPointEdgeX[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i - 1, j, pt.x, pt.y, &par_m);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GasConstant / m_mol;
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
                flx[0] = 0.5 * ((par_p.r * par_p.u + par_m.r * par_m.u) - alpha * (par_p.r - par_m.r));
                flx[1] = 0.5 * ((par_p.r * par_p.u * par_p.u + par_p.p + par_m.r * par_m.u * par_m.u + par_m.p)
                                - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[2] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[3] = 0.5 *
                         (((par_p.r * par_p.e_tot + par_p.p) * par_p.u + (par_m.r * par_m.e_tot + par_m.p) * par_m.u)
                          - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        int_m[i_fld][k] += gaussWeightEdgeX[i][j][i_gp] * flx[i_fld] * bf(k, i - 1, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    int_m[i_fld][k] *= gj_edg_x[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i - 1][j].fields[i_fld][k] -= int_m[i_fld][k];
                }
            }
        }
    }

    // Bottom
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j <= 0; j++) {
            double int_m[4][NumberBasisFunctions], int_p[4][NumberBasisFunctions];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_p[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                point_t pt = gaussPointEdgeY[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j, pt.x, pt.y, &par_p);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GasConstant / m_mol;
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

                flx[0] = 0.5 * ((par_p.r * par_p.v + par_m.r * par_m.v) - alpha * (par_p.r - par_m.r));
                flx[1] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[2] = 0.5 * ((par_p.r * par_p.v * par_p.v + par_p.p + par_m.r * par_m.v * par_m.v + par_m.p)
                                - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[3] = 0.5 *
                         (((par_p.r * par_p.e_tot + par_p.p) * par_p.v + (par_m.r * par_m.e_tot + par_m.p) * par_m.v)
                          - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        int_p[i_fld][k] += gaussWeightEdgeY[i][j][i_gp] * flx[i_fld] * bf(k, i, j, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    int_p[i_fld][k] *= gj_edg_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i][j].fields[i_fld][k] += int_p[i_fld][k];
                }
            }
        }
    }


    // Top
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = NumberCellsY; j <= NumberCellsY; j++) {
            double int_m[4][NumberBasisFunctions], int_p[4][NumberBasisFunctions];
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                memset(int_m[i_fld], 0, sizeof(double) * NumberBasisFunctions);
            }
            for (int i_gp = 0; i_gp < NumberGaussPointsEdge; i_gp++) {
                point_t pt = gaussPointEdgeY[i][j][i_gp];
                prim_t par_m, par_p;
                double flx[4];
                cons_to_prim(i, j - 1, pt.x, pt.y, &par_m);
                {
                    double cp = 1014.16;
                    double m_mol = 0.02869409;
                    double cv = cp - GasConstant / m_mol;
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

                flx[0] = 0.5 * ((par_p.r * par_p.v + par_m.r * par_m.v) - alpha * (par_p.r - par_m.r));
                flx[1] = 0.5 * ((par_p.r * par_p.u * par_p.v + par_m.r * par_m.u * par_m.v)
                                - alpha * (par_p.r * par_p.u - par_m.r * par_m.u));
                flx[2] = 0.5 * ((par_p.r * par_p.v * par_p.v + par_p.p + par_m.r * par_m.v * par_m.v + par_m.p)
                                - alpha * (par_p.r * par_p.v - par_m.r * par_m.v));
                flx[3] = 0.5 *
                         (((par_p.r * par_p.e_tot + par_p.p) * par_p.v + (par_m.r * par_m.e_tot + par_m.p) * par_m.v)
                          - alpha * (par_p.r * par_p.e_tot - par_m.r * par_m.e_tot));

                for (int i_fld = 0; i_fld < 4; i_fld++) {
                    for (int k = 0; k < NumberBasisFunctions; k++) {
                        int_m[i_fld][k] += gaussWeightEdgeY[i][j][i_gp] * flx[i_fld] * bf(k, i, j - 1, pt.x, pt.y);
                    }
                }
            }
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    int_m[i_fld][k] *= gj_edg_y[i][j];
                }
            }

            for (int i_fld = 0; i_fld < 4; i_fld++) {
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    r_int[i][j - 1].fields[i_fld][k] -= int_m[i_fld][k];
                }
            }
        }
    }
}

void calc_new() {
    double tmp[NumberBasisFunctions];
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                mult_matr_vec(matr_a[i][j], r_int[i][j].fields[i_fld], tmp);
                for (int k = 0; k < NumberBasisFunctions; k++) {
                    data[i][j].fields[i_fld][k] += Tau * tmp[k];
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

    for (int i = 1; i < NumberCellsX - 1; i++) {
        for (int j = 0; j < NumberCellsY; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                double fm = data[i - 1][j].fields[i_fld][0];
                double fc = data[i][j].fields[i_fld][0];
                double fp = data[i + 1][j].fields[i_fld][0];
                data[i][j].fields[i_fld][1] = _minmod(data[i][j].fields[i_fld][1], LimiterAlpha * (fc - fm),
                                                      LimiterAlpha * (fp - fc));
            }
        }
    }
    for (int i = 0; i < NumberCellsX; i++) {
        for (int j = 1; j < NumberCellsY - 1; j++) {
            for (int i_fld = 0; i_fld < 4; i_fld++) {
                double fm = data[i][j - 1].fields[i_fld][0];
                double fc = data[i][j].fields[i_fld][0];
                double fp = data[i][j + 1].fields[i_fld][0];
                data[i][j].fields[i_fld][2] = _minmod(data[i][j].fields[i_fld][2], LimiterAlpha * (fc - fm),
                                                      LimiterAlpha * (fp - fc));
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

double get_fld(int i_fld, int i, int j, double x, double y) {
    double res = 0.0;
    for (int i_bf = 0; i_bf < NumberBasisFunctions; i_bf++) {
        res += data[i][j].fields[i_fld][i_bf] * bf(i_bf, i, j, x, y);
    }
    return res;
}

void cons_to_prim(int i, int j, double x, double y, prim_t *prim) {
    double *concentration = concentrations[i][j];
    double cp = get_cell_cp(concentration);
    double m_mol = get_cell_M(concentration);

    double cv = cp - GasConstant / m_mol;
    double gam = cp / cv;

    double ro = get_fld(0, i, j, x, y);
    double ru = get_fld(1, i, j, x, y);
    double rv = get_fld(2, i, j, x, y);
    double re = get_fld(3, i, j, x, y);

    prim->r = ro;
    prim->u = ru / ro;
    prim->v = rv / ro;
    prim->e_tot = re / ro;
    prim->e = prim->e_tot - (prim->u * prim->u + prim->v * prim->v) * 0.5;
    prim->p = prim->r * prim->e * (gam - 1.0);
    prim->cz = sqrt(gam * prim->p / prim->r);
    prim->t = prim->e / cv;
}

double get_component_cp(int id) {
    double cps[NumberComponents] = {1014.16, 1014.16, 1014.16, 1014.16, 1014.16};
    return cps[id];
}

double get_component_M(int id) {
    double Ms[NumberComponents] = {0.02869409, 0.02869409, 0.02869409, 0.02869409, 0.02869409};
    return Ms[id];
}

double get_cell_cp(double concentration[NumberComponents]) {
    double cp = 0;
    for (int i_com = 0; i_com < NumberComponents; i_com++) {
        cp += concentration[i_com] * get_component_cp(i_com);
    }
    return cp;
}

double get_cell_M(double concentration[NumberComponents]) {
    double M = 0;
    for (int i_com = 0; i_com < NumberComponents; i_com++) {
        M += concentration[i_com] / get_component_M(i_com);
    }
    return 1 / M;
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
    fprintf(fp, "DIMENSIONS %d %d %d\n", NumberCellsX + 1, NumberCellsY + 1, 1);
    fprintf(fp, "POINTS %d float\n", (NumberCellsX + 1) * (NumberCellsY + 1));
    for (int j = 0; j <= NumberCellsY; j++) {
        for (int i = 0; i <= NumberCellsX; i++) {
            double x = XMin + i * HX;
            double y = YMin + j * HY;
            fprintf(fp, "%f %f %f\n", x, y, 0.0);
        }
    }
    fprintf(fp, "CELL_DATA %d\n", NumberCellsX * NumberCellsY);
    fprintf(fp, "SCALARS Density float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.r);
        }
    }

    fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.p);
        }
    }

    fprintf(fp, "SCALARS TotEnergy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e_tot);
        }
    }

    fprintf(fp, "SCALARS Energy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e);
        }
    }

    fprintf(fp, "VECTORS Velosity float\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f %f %f\n", pr.u, pr.v, 0.0);
        }
    }

    fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            cons_to_prim(i, j, centerCell[i][j].x, centerCell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.t);
        }
    }

    fprintf(fp, "VECTORS Concentrations float\n");
    for (int j = 0; j < NumberCellsY; j++) {
        for (int i = 0; i < NumberCellsX; i++) {
            for (int k = 0; k < NumberComponents - 1; k++) {
                fprintf(fp, "%f ", concentrations[i][j][k]);
            }
            fprintf(fp, "%f\n", concentrations[i][j][NumberComponents - 1]);
        }
    }

    fclose(fp);
    printf("File '%s' saved...\n", fName);

}