//
// Created by nzvi on 19.05.2020.
//
#include "include/basis.h"

double bf(int i_func, int i, int j, double x, double y) {
    switch (i_func) {
        case 0:
            return 1.0;
            break;
        case 1:
            return (x - center_cell[i][j].x) / HX;
            break;
        case 2:
            return (y - center_cell[i][j].y) / HY;
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
