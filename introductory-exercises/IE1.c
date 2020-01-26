/*  TODO:
    - don't overcomplicate IE1...
    - function pointer to ODE function?
    - write results within EU() ?
    - write no smaller stepsize than h=0.05s
*/

#include <math.h>           // ??
#include <stdio.h>          // IO: fprintf(), etc.
#include <stdbool.h>        // bool

double ODEs(int i, double t, const double var[1]){
    // var[0]: x(t)
    double x = var[0];  // [m]
    double alpha = 1;   // [s^-1]

    switch(i){
        case 0:
            // xdot = -alpha*x
//            printf("%f\n",-alpha*x);
            return -alpha*x;
        default: 
            printf("invalid argument i=%d for ODEs()\n", i);
            return 0;
    }
}

//bool stopIf();  // function for breaking condition of integration?

void EU(int nSteps, int nVar, double t, double h, double var[nVar], 
        double tArr[nVar], double results[nVar][nSteps]){

    // save initial conditions
    tArr[0] = t;
    for(int i=0;i<nVar;i++){
        results[i][0] = var[i];
    }


    // calculate steps
    for(int n=1;n<=nSteps;n++){
        for(int i=0;i<nVar;i++){
            var[i] += h*ODEs(i,t,var);
        }
        t += h;

        // save results
        tArr[n] = t;
        for(int i=0;i<nVar;i++){
            results[i][n] = var[i];
        }
    }
}

int main(){
    // set stepsize h, start and finish time
    double h = 0.01;         // [s]
    double t_0 = 0.0;       // [s]
    double t_f = 20.0;      // [s]
    
    // calculate number of steps N
    int nSteps = (t_f-t_0)/h + 1;   
    int nVar   = 1;

    // for results
    double results[nVar][nSteps];
    double tArr[nSteps];

    // initial conditions
    double var0[nVar]; 
    var0[0] = 1.0;        // [m]

    // calculate
    EU(nSteps, nVar, t_0, h, var0, tArr, results);

    // write results into file
    FILE * output;
    output = fopen("output/IE1.dat", "w+");
    fprintf(output, "# t[s]         x[m]         \n");
    for(int n=0;n<=nSteps;n++){
        double h_min = 0.5;  // write no smaller step size than h_min
        if(h<h_min && n%((int)(h/h_min)+1)!=1){
            n += (int)(h/h_min)+1;
            printf("%i\n",n);
        }
        fprintf(output, "  %.3e", tArr[n]);
        for(int i=0; i<nVar;i++){
            fprintf(output, "    %.3e\n", results[i][n]);
        }
    }
    fclose(output);

    return 0;
}