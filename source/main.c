// testing environment for Runge Kutta

#include <math.h>           // compile with -lm
#include <stdio.h>          // IO: fprintf(), etc.
#include "Runge_Kutta_4.h"  // RK4()
#include "ODEs.h"
#include "read_parameters.h"

int main(){
    const int N = 3;  // number of objects
    const int D = 2;  // number of dimensions

    double parameters[2*N*D+N+5];
    read_parameters("input.dat", N, D, parameters);
    
    double m[N];// masses of objects
    for(int i=0;i<N;i++){
        m[i] = parameters[2*N*D +i];
    }

    // set stepsize h, start and finish time
    const double h = parameters[2*N*D+N +0];         // [s]
    const double h_min = parameters[2*N*D+N +1]; // minimum step size written
    const double precission = parameters[2*N*D+N +2];  
    // ^ needs to be sufficiently small; wont work otherwise (no idea why..)
    const double t_0 = parameters[2*N*D+N +3];      // [s]
    const double t_f = parameters[2*N*D+N +4];      // [s]

    int nVar = 2*N*D; 
    // for results
    double var[nVar];
    double t;

    // set initial conditions
    for(int i=0;i<2*N*D;i++){
        var[i] = parameters[i];
    }

    t = t_0;    // [s]

    // check input
    for(int i=0;i<2*N*D+N+5;i++){
        printf("%f \n", parameters[i]);
    }
  
    // write head and initial conditions
    FILE* output;
    char path[50]; 
    sprintf(path, "../output/data.dat",h);
    output = fopen(path, "w+");
    fprintf(output, "# t[s]            ");
    fprintf(output, "x_00[m]         ");
    fprintf(output, "x_01[m]         ");
    fprintf(output, "v_00[m/s]       ");
    fprintf(output, "v_01[m/s]       ");
    fprintf(output, "x_10[m]         ");
    fprintf(output, "x_11[m]         ");
    fprintf(output, "v_10[m/s]       ");
    fprintf(output, "v_11[m/s]       ");
    fprintf(output, "x_20[m]         ");
    fprintf(output, "x_21[m]         ");
    fprintf(output, "v_20[m/s]       ");
    fprintf(output, "v_21[m/s]       ");
    fprintf(output, "\n");
    // 
    fprintf(output, " %+.6e", t); 
    for(int i=0; i<nVar;i++){
        fprintf(output, "   %+.6e", var[i]);
    } 
    fprintf(output, "\n");

    // calculate
    RK4(&ODEs, h, &t,nVar, var, &output, h_min); 

    // don't forget to close file!
    fclose(output);

    // now lets find the time 

    return 0;
}
