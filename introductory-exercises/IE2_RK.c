// testing environment for Runge Kutta

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

void EU(int nSteps, double h, double *t, int nVar, double x[][nVar]){

    // calculate steps
    int n = 0;
    while(n<nSteps-1){
        for(int i=0;i<nVar;i++){
            x[n+1][i] = x[n][i] + h*ODEs(i,*t,x[n]);
        }
        *t += h;

        n += 1;   // increase step counter
    }
}


void RK4(int nSteps, double h, double *t, int nVar, double x[][nVar]){

    double k1[nVar], k2[nVar], k3[nVar], k4[nVar];  // k_j for each variable

    int n = 0;  // start with first step (n=0 : initial conditions)
    while(n<nSteps-1){
        double temp_x [nVar]; // temporar copies of x[n][i] to calculate k_j[i]
        
        int i;
        //k1:
        for(i=0;i<nVar;i++){
            temp_x[i] = x[n][i];
            k1[i] = ODEs(i,*t,temp_x);
        }
        //k2:
        for(i=0;i<nVar;i++){
            temp_x[i] = x[n][i] + h/2*k1[i];
            k2[i] = ODEs(i,*t+h/2, temp_x);
        }
        //k3:
        for(i=0;i<nVar;i++){
            temp_x[i] = x[n][i] + h/2*k2[i];
            k3[i] = ODEs(i,*t+h/2, temp_x);
        }
        //k4:
        for(i=0;i<nVar;i++){
            temp_x[i] = x[n][i] + h*k3[i];
            k4[i] = ODEs(i,*t+h, temp_x);
        }

        // make final step
        for(i=0;i<nVar;i++){
            x[n+1][i] = x[n][i] + h/6*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
        }
        *t += h;

        n += 1;   // increase step counter
    }
}


int main(){
    // set stepsize h, start and finish time
    const double h = 0.1;         // [s]
    const double t_0 = 0.0;       // [s]
    const double t_f = 20.0;      // [s]
    
    // calculate number of steps N
    const int nSteps = (t_f-t_0)/h;   
    const int nVar   = 2;

    // for results
    double x[nSteps+1][nVar];
    double t;

    // set initial conditions
    x[0][0] = 0; // [m]
    x[0][1] = 1; // [m/s]
    t = t_0;     // [s]

    // calculate
    RK4(nSteps, h, &t, nVar, x);

    // write results into file

    double h_min = 0.05; // minimum step size writen
    //double nSteps_write = nSteps*(h/h_min);
    int n_skip = h_min/h;    // n skipped to stick to h_min

    FILE* output;
    char path[50];
    sprintf(path, "output/IE2_RK_h=%.3f.dat",h);
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