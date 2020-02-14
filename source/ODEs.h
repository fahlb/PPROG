#include <math.h>
#include <stdio.h>

double r_ij_d(int i, int j, int d, int D, const double var[]);

double distance_ij(int i, int j, int D, const double var[]);

double T_i(int i, int D, double m[], const double var[]);

double V_ij(int i, int j, int D, double m[], const double var[]);

double E_tot(int N, int D, double m[], const double var[]);

double distance2cm(int i, int N, int D, double m[], const double var[]);

double r_polar_i(int i, int N, int D, const double var[]);

double ODEs(int k, double t, const double var[]);