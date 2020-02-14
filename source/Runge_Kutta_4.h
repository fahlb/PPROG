#include <math.h>
#include <stdio.h>
#include <stdbool.h>

void NextEU(double (*f)(int, double, const double*),  // function pointer
            double h, double *t, int nVar, double x[]);

void EU(double (*f)(int, double, const double[]),  // function pointer
        double h, double *t, int nVar, double x[nVar], 
        FILE** output, double h_min);

void NextRK4(double (*f)(int, double, const double*),  // function pointer
             double h, double *t, int nVar, double x[]);

void RK4(double (*f)(int, double, const double[]),  // function pointer
        double h, double *t, int nVar, double x[nVar], 
        FILE** output, double h_min);


void AdaptiveRK4(double (*f)(int, double, const double[]),  // function pointer
                double h, double *t, int nVar, double x[nVar], 
                FILE** output, double h_min, double precission);

// 4-th order RK with adaptive step size based on conservatino of total energy:
void energyAdaptiveRK4(double (*f)(int, double, const double[]),
                double h, double *t, int nVar, double x[nVar], 
                FILE** output, double h_min, double precision);