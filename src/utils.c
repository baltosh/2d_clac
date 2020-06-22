//
// Created by nzvi on 19.05.2020.
//

#include "include/utils.h"

void zero_r() {
    for (int i = 0; i < CELLS_X_COUNT; i++) {
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            memset(&r_int[i][j], 0, sizeof(data_t));
        }
    }
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

double _minmod(double a, double b, double c) {
    if ((_SIGN_(a) == _SIGN_(b)) && (_SIGN_(b) == _SIGN_(c))) {
        return _SIGN_(a) * _MIN_(_MIN_(fabs(a), fabs(b)), fabs(c));
    } else {
        return 0.0;
    }
}

void save_vtk(int num) {
    char fName[50];
    prim_t pr;
    memset(&pr, 0, sizeof(prim_t));
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
            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.r);
        }
    }

    fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.p);
        }
    }

    fprintf(fp, "SCALARS TotEnergy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e_tot);
        }
    }

    fprintf(fp, "SCALARS Energy float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.e);
        }
    }

    fprintf(fp, "VECTORS Velosity float\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
            fprintf(fp, "%f %f %f\n", pr.u, pr.v, 0.0);
        }
    }

    fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n");
    for (int j = 0; j < CELLS_Y_COUNT; j++) {
        for (int i = 0; i < CELLS_X_COUNT; i++) {
            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
            fprintf(fp, "%f\n", pr.t);
        }
    }

//    fprintf(fp, "VECTORS Concentrations float\n");
//    for (int j = 0; j < CELLS_Y_COUNT; j++) {
//        for (int i = 0; i < CELLS_X_COUNT; i++) {
//            cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
//            for (int k = 0; k < COMPONENTS_COUNT - 1; k++) {
//                fprintf(fp, "%f ", pr.c[k]);
//            }
//            fprintf(fp, "%f\n", pr.c[COMPONENTS_COUNT - 1]);
//        }
//    }
    for (int k = 0; k < COMPONENTS_COUNT; k++) {
        fprintf(fp, "SCALARS Concentration%i float 1\nLOOKUP_TABLE default\n", k);
        for (int j = 0; j < CELLS_Y_COUNT; j++) {
            for (int i = 0; i < CELLS_X_COUNT; i++) {
                cons_to_prim(i, j, center_cell[i][j].x, center_cell[i][j].y, &pr);
                fprintf(fp, "%f\n", pr.c[k]);
            }
        }
    }

    fclose(fp);
    printf("File '%s' saved...\n", fName);
}
