// testing environment for Runge Kutta

#include <math.h>           // compile with -lm
#include <stdio.h>          // IO: fprintf(), etc.
#include "Runge_Kutta_4.h"  // RK4()
#include "ODEs.h"

int main(){
    int N = 3;  // number of objects
    int D = 2;  // number of dimensions
    
    double m[N];// masses of objects
    m[0] = 5.688e26;    // M_Saturn     [kg]
    m[1] = 1.98e18;     // M_Janus      [kg]
    m[2] = 5.5e17;      // M_Epimetheus [kg]

    // set stepsize h, start and finish time
    const double h = 60;         // [s]
    const double h_min = 1*(24*60*60); // minimum step size to be written
    const double precission = 1e-6;  // precission when using adaptive step size
    // ^ needs to be sufficiently small; wont work otherwise (no idea why..)
    const double t_0 = 0.0;       // [s]
    const double t_f = 10*365.25*60*60;      // [s]

    const int nVar = 12; 
    // for results
    double var[nVar];
    double t;

    // set initial conditions
    var[0] = 0;         // x_01(t) Saturn
    var[1] = 0;         // x_02(t)
    var[2] = 0;         // v_01(t)
    var[3] = 0;         // v_02(t)
    var[4] = 151472e3;  // x_11(t) Janus
    var[5] = 0;         // x_12(t)
    var[6] = 0;         // v_11(t)
    var[7] = 15830.8;   // v_12(t)
    var[8] = -151422e3; // x_21(t) Epimetheus
    var[9] = 0;         // x_22(t)
    var[10]= 0;         // v_21(t)
    var[11]= -15833.4;  // v_22(t)


    t = t_0;    // [s]
  
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
