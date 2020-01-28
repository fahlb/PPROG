// testing environment for Runge Kutta

#include <math.h>           // ??
#include <stdio.h>          // IO: fprintf(), etc.
#include "Runge_Kutta_4.c"  // RK4()
#include "ODEs.c"

int main(){
    // set stepsize h, start and finish time
    const double h = 0.1;         // [s]
    const double h_min = 0.05; // minimum step size to be written

    const double t_0 = 0.0;       // [s]
    const double t_f = 20.0;      // [s]

    // calculate number of steps N
    const int nVar = 12;

    // for results
    double x[nVar];
    double t;

    // set initial conditions
    x[0] = 0;   // [m]
    x[1] = 1;   // [m/s]
    t = t_0;    // [s]
  
    // write head and initial conditions
    FILE* output;
    char path[50]; 
    sprintf(path, "../output/data.dat",h);
    output = fopen(path, "w+");
    fprintf(output, "# t[s]         ");
    fprintf(output, "r_00[m]      ");
    fprintf(output, "r_01[m]      ");
    fprintf(output, "v_00[m/s]    ");
    fprintf(output, "v_01[m/s]    ");
    fprintf(output, "r_10[m]      ");
    fprintf(output, "r_11[m]      ");
    fprintf(output, "v_10[m/s]    ");
    fprintf(output, "v_11[m/s]    ");
    fprintf(output, "r_20[m]      ");
    fprintf(output, "r_21[m]      ");
    fprintf(output, "v_20[m/s]    ");
    fprintf(output, "v_21[m/s]    ");
    fprintf(output, "\n");
    // 
    fprintf(output, "  %.3e", t); 
    for(int i=0; i<nVar;i++){
        fprintf(output, "    %.3e", x[i]);
    } 
    fprintf(output, "\n");


    // calculate
    RK4(&ODEs, h, &t,nVar, x, &output, h_min); 

    // don't forget to close file!
    fclose(output);

    return 0;
}
