//
// Created by nzvi on 16.05.2020.
//

#ifndef INC_2D_CALC_GLOBAL_H
#define INC_2D_CALC_GLOBAL_H


#define M_PI        3.14159265358979323846
#define GAS_CONSTANT 8.31446261815324

#define BASE_FN_COUNT    3
#define CELLS_X_COUNT            80
#define CELLS_Y_COUNT            800
#define COMPONENTS_COUNT        4

#define GP_EDGE_COUNT 2
#define GP_CELL_COUNT 4

#define X_MIN  0.0
#define X_MAX  0.04
#define Y_MIN -0.1
#define Y_MAX  0.3

#define TAU     5.0e-8
#define TIME_MAX 2.5e-3

#define SAVE_STEP 100

#define LIMITER_ALPHA 2.0

#define _SIGN_(X) (fabs(X)/(X))
#define _MIN_(X, Y) ((X)<(Y) ? (X) : (Y))
#define _MAX_(X, Y) ((X)>(Y) ? (X) : (Y))


typedef struct prim {
    double r;             /**< density        */
    double u;             /**< velocity       */
    double v;             /**< velocity       */
    double e;             /**< energy         */
    double e_tot;         /**< total energy   */
    double p;             /**< pressure       */
    double t;             /**< temperature    */
    double cz;            /**< sound velocity */
    double gam;
    double cp;
    double cv;
    double m;
    double c[COMPONENTS_COUNT]; /**< concentrations */
} prim_t;


typedef struct cons {
    union {
        double fields[3 + COMPONENTS_COUNT][BASE_FN_COUNT];
        struct {
            double rc[COMPONENTS_COUNT][BASE_FN_COUNT];
            double ru[BASE_FN_COUNT];
            double rv[BASE_FN_COUNT];
            double re[BASE_FN_COUNT];
        };
    };

} data_t;

typedef struct {
    double x, y;
} point_t;

double HX, HY;

data_t data[CELLS_X_COUNT][CELLS_Y_COUNT], r_int[CELLS_X_COUNT][CELLS_Y_COUNT];

point_t gp_edge_x[CELLS_X_COUNT + 1][CELLS_Y_COUNT][GP_EDGE_COUNT];
point_t gp_edge_y[CELLS_X_COUNT][CELLS_Y_COUNT + 1][GP_EDGE_COUNT];
double gw_edge_x[CELLS_X_COUNT + 1][CELLS_Y_COUNT][GP_EDGE_COUNT];
double gw_edge_y[CELLS_X_COUNT][CELLS_Y_COUNT + 1][GP_EDGE_COUNT];
double gj_edge_x[CELLS_X_COUNT + 1][CELLS_Y_COUNT];
double gj_edge_y[CELLS_X_COUNT][CELLS_Y_COUNT + 1];

point_t gp_cell[CELLS_X_COUNT][CELLS_Y_COUNT][GP_CELL_COUNT];
double gw_cell[CELLS_X_COUNT][CELLS_Y_COUNT][GP_CELL_COUNT];
double gj_cell[CELLS_X_COUNT][CELLS_Y_COUNT];
point_t centerCell[CELLS_X_COUNT][CELLS_Y_COUNT];

point_t q_point[4][2];
point_t q_gp_cell[4][GP_CELL_COUNT];
double q_gw_cell[4][GP_CELL_COUNT];
double q_gj_cell[4];



double matr_a[CELLS_X_COUNT][CELLS_Y_COUNT][BASE_FN_COUNT][BASE_FN_COUNT];

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

double get_field_ru(int i, int j, double x, double y);

double get_field_rv(int i, int j, double x, double y);

double get_field_re(int i, int j, double x, double y);

double get_field_rc(int i, int j, double x, double y, int e);

void cons_to_prim(int i, int j, double x, double y, prim_t *prim);

double get_component_cp(int id);

double get_component_M(int id);

double get_cell_cp(int i, int j);

double get_cell_M(int i, int j);

double reactionSpeed0(double[]);

double reactionSpeed1(double[]);

double reactionSpeed2(double[]);

double reactionSpeed3(double[]);

double reactionSpeed4(double[]);

double (*reactionSpeeds[COMPONENTS_COUNT])(double[]) = {reactionSpeed0,
                                                        reactionSpeed1,
                                                        reactionSpeed2,
                                                        reactionSpeed3,
                                                        reactionSpeed4};

#endif //INC_2D_CALC_GLOBAL_H
