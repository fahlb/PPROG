// testing environment for Runge Kutta

#include <math.h>           // ??
#include <stdio.h>          // IO: fprintf(), etc.
#include "../source/Runge_Kutta_4.h"  // RK4()

/* TODO:
    - write in rk functions, pass FILE* if necessary and just 
      wirte head of file in main()
      don't forget to close in main()

*/

double ODEs(int i, double t, const double var[]){
    // var[0]: x(t)
    // var[1]: v(t)

    double x = var[0];  // [m]
    double v = var[1];  // [m/s]
    double omega = 1;   // [1/s]

    switch(i){
        case 0:
            // dx/dt = v(t)
            return v;
        case 1: 
            // dv/dt = -omega^2 * x
            return - omega*omega*x;
        case -1:
            // use as doIf() ?
            if(t<20.0){return 1;}
            if(t>=20.0){return 0;}
        default: 
            printf("invalid argument i=%d for ODEs()\n", i);
            return 0;
    }
}

int main(){
    // set stepsize h, start and finish time
    const double h = 0.1;         // [s]
    const double h_min = 0.05; // minimum step size to be written
 
    const double t_0 = 0.0;       // [s]
    const double t_f = 20.0;      // [s]
    
    // calculate number of steps N
    const int nVar   = 2;

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
//    sprintf(path, "output/IE2_RK_h=%.3f.dat",h);
    sprintf(path, "output/IE2_RK_adaptive.dat",h);
  
    output = fopen(path, "w+");
//    output = fopen("../output/test.dat", "w+"); 
    fprintf(output, "# t[s]         ");
    fprintf(output, "x[m]         "  );
    fprintf(output, "v[m/s]       "  );
    fprintf(output, "\n");
    // 
    fprintf(output, "  %.3e", t); 
    for(int i=0; i<nVar;i++){
        fprintf(output, "    %.3e", x[i]);
    } 
    fprintf(output, "\n");


    // calculate
    AdaptiveRK4(&ODEs, h, &t,nVar, x, &output, h_min,0.001); 

    // don't forget to close file! - Maybe not necessary; closes automatically
    fclose(output);

    return 0;
}
