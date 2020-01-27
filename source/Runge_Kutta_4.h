#include <math.h>
#include <stdio.h>

void EU(int nSteps, double h, double *t, int nVar, double x[][nVar]);
void RK4(int nSteps, double h, double *t, int nVar, double x[][nVar]);