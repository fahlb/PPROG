#include <math.h>
#include <stdio.h>

// global so it only needs to be read once in main.c and can be used in ODEs()
extern double m_input[];
extern double t_f_input;

double r_ij_d(int i, int j, int d, int D, const double var[]);

double distance_ij(int i, int j, int D, const double var[]);

double T_i(int i, int D, double m[], const double var[]);

double V_ij(int i, int j, int D, double m[], const double var[]);

double E_tot(int N, int D, double m[], const double var[]);

double distance2cm(int i, int N, int D, double m[], const double var[]);

double r_polar_i(int i, int N, int D, const double var[]);

double ODEs(int k, double t, const double var[]);