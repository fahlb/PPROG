#include <math.h>           // ??
#include <stdio.h>          // IO: fprintf(), etc.

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
        default: 
            printf("invalid argument i=%d for ODEs()\n", i);
            return 0;
    }
}

//bool stopIf();  // function for breaking condition of integration?

void EU(int nSteps, int nVar, double h, double *t, double x[][nVar]){

    // calculate steps
    int n = 1;
    while(n<nSteps){
        for(int i=0;i<nVar;i++){
            x[n][i] = x[n-1][i] + h*ODEs(i,*t,x[n-1]);
        }
        *t += h;

        n += 1;   // increase step counter
    }
}

int main(){
    // set stepsize h, start and finish time
    const double h = 0.001;         // [s]
    const double t_0 = 0.0;       // [s]
    const double t_f = 20.0;      // [s]
    
    // calculate number of steps N
    const int nSteps = (t_f-t_0)/h;   
    const int nVar   = 2;

    // for results
    double x[nSteps][nVar];
    double t;

    // set initial conditions
    x[0][0] = 0; // [m]
    x[0][1] = 1; // [m/s]
    t = t_0;     // [s]

    // calculate
    EU(nSteps, nVar, h, &t, x);

    // write results into file

    double h_min = 0.05; // minimum step size writen
    //double nSteps_write = nSteps*(h/h_min);
    int n_skip = h_min/h;    // n skipped to stick to h_min

    FILE* output;
    char path[50];
    sprintf(path, "output/IE2_h=%.3f.dat",h);
    output = fopen(path, "w+");
    fprintf(output, "# t[s]         ");
    fprintf(output, "x[m]         "  );
    fprintf(output, "v[m/s]       "  );
    fprintf(output, "\n");

    for(int n=0;n<nSteps;n++){
        if(h<h_min && (n%n_skip)){continue;} //skip if not multiple of n_skip
        fprintf(output, "  %.3e", n*h);
        for(int i=0; i<nVar;i++){
            fprintf(output, "    %.3e", x[n][i]);
        }
        fprintf(output, "\n");
    }
    fclose(output);

    return 0;
}