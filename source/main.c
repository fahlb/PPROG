// testing environment for Runge Kutta

#include <math.h>           // ??
#include <stdio.h>          // IO: fprintf(), etc.
#include "Runge_Kutta_4.c"  // RK4()
#include "ODEs.c"

int main(){
    // set stepsize h, start and finish time
    const double h = 60;         // [s]
    const double h_min = 24*3600;      // minimum step size to be written

    const double t_0 = 0.0;       // [s]
    const double t_f = 31536000;      // [s]

    // calculate number of steps N
//    const int nVar = 12;
    const int nVar = 12; // testing with just one moon

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
    var[7] = 15831.3;   // v_12(t)
    var[8] = -151422e3; // x_21(t) Epimetheus
    var[9] = 0;         // x_22(t)
    var[10]= 0;         // v_21(t)
    var[11]= -15833.9;  // v_22(t)


    t = t_0;    // [s]
  
    // write head and initial conditions
    FILE* output;
    char path[50]; 
    sprintf(path, "../output/data.dat",h);
    output = fopen(path, "w+");
    fprintf(output, "# t[s]         ");
    fprintf(output, "x_00[m]      ");
    fprintf(output, "x_01[m]      ");
    fprintf(output, "v_00[m/s]    ");
    fprintf(output, "v_01[m/s]    ");
    fprintf(output, "x_10[m]      ");
    fprintf(output, "x_11[m]      ");
    fprintf(output, "v_10[m/s]    ");
    fprintf(output, "v_11[m/s]    ");
    fprintf(output, "x_20[m]      ");
    fprintf(output, "x_21[m]      ");
    fprintf(output, "v_20[m/s]    ");
    fprintf(output, "v_21[m/s]    ");
    fprintf(output, "\n");
    // 
    fprintf(output, "  %.3e", t); 
    for(int i=0; i<nVar;i++){
        fprintf(output, "    %.3e", var[i]);
    } 
    fprintf(output, "\n");

    // calculate
    RK4(&ODEs, h, &t,nVar, var, &output, h_min); 

    // don't forget to close file!
    fclose(output);

    return 0;
}
