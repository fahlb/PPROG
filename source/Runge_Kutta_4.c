#include <math.h>
#include <stdio.h>
#include <stdbool.h>
//#include "Runge_Kutta_4.h"

//double ODEs(int i, double t, const double var[]);

/* TODO:
    - time is broken, see below:
    - In order to make the AdaptiveRK() work I need to change t* to tArr
      and save each time step. Furthermore I need to keep track of
      the number of steps taken or, alternatively, check the size of
      tArr for example to get the steps it took in the end.
      How do I ensure the calculation runs for the time specified?
      * use static nSteps to get time to stop?
      * use completely different method, maybe like in bachelor code?
      * hardcode ...
    - Just fuck adaptive step size for now? (but would be useful for BA, too..)

*/

void EU(double (*f)(int, double, const double[]),  // function pointer
        int nSteps, double h, double *t, int nVar, double x[][nVar]){

    // calculate steps
    int n = 0;
    while(n<nSteps-1){
        for(int i=0;i<nVar;i++){
            x[n+1][i] = x[n][i] + h*(*f)(i,*t,x[n]);
        }
        *t += h;

        n += 1;   // increase step counter
    }
}
void NextRK4(double (*f)(int, double, const double*),  // function pointer
             int n, double h, double *t, int nVar, double x[]){

    double k1[nVar], k2[nVar], k3[nVar], k4[nVar];  // k_j for each variable
    double temp_x [nVar]; // temporar copies of x[n][i] to calculate k_j[i]
    int i;

    //k1:
    for(i=0;i<nVar;i++){
        temp_x[i] = x[i];
        k1[i] = (*f)(i,*t,temp_x);
    }
    //k2:
    for(i=0;i<nVar;i++){
        temp_x[i] = x[i] + h/2*k1[i];
        k2[i] = (*f)(i,*t+h/2, temp_x);
    }
    //k3:
    for(i=0;i<nVar;i++){
        temp_x[i] = x[i] + h/2*k2[i];
        k3[i] = (*f)(i,*t+h/2, temp_x);
    }
    //k4:
    for(i=0;i<nVar;i++){
        temp_x[i] = x[i] + h*k3[i];
        k4[i] = (*f)(i,*t+h, temp_x);
    }

    // make final step
    for(i=0;i<nVar;i++){
        x[i] = x[i] + h/6*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
    }
//    *t += h;


}


void RK4(double (*f)(int, double, const double*),  // function pointer
         int nSteps, double h, double *t, int nVar, double x[][nVar]){

    int n = 0;  // start with first step (n=0 : initial conditions)
//    while(n<nSteps-1){
    while((*f)(-1,*t,x[n])){    // breaks code :(   


        double temp_copy_x[nVar];
        for(int i=0;i<nVar;i++){
            temp_copy_x[i] = x[n][i];
        }
        
        NextRK4(f,n,h,t,nVar,temp_copy_x);
        *t += h;
        for(int i=0;i<nVar;i++){
            x[n+1][i] = temp_copy_x[i];
        }
        printf("last n: %d \n", n);
        n += 1;   // increase step counter
    }
}

void AdaptiveRK4(double (*f)(int, double, const double*),  // function pointer
         int nSteps, double h, double *t, int nVar, double x[][nVar]){

    double max_error[nVar]; // make argument; maximum error for each variable
    for(int i=0;i<nVar;i++){
        max_error[i] = 1e-2;
    }   // temporary!!!

    int n = 0;  // start with first step (n=0 : initial conditions)
    double x_1full[nVar];   // temporary x[i] for 1 full step
    double x_2half[nVar];   // temporary x[i] for 2 half steps
    
//    while(n<nSteps-1){
    while((*f)(-1,*t,x[n])){    // breaks code :(   
        // save copy of x[i] for both advances
        for(int i=0;i<nVar;i++){
            x_1full[i] = x[n][i];
            x_2half[i] = x[n][i];
        }
        // make one full step:
        NextRK4(f,n,h,t,nVar,x_1full);
        // make two half steps:
        NextRK4(f,n,h/2,t,nVar,x_2half);
        double newt = *t+h/2;
        NextRK4(f,n+1,h/2,&(newt),nVar,x_2half);
        
        double difference[nVar];    // difference of 1 full and 2 half steps
        bool take_step = true;
        for(int i=0;i<nVar;i++){
            difference[i] = fabs(x_1full[i]-x_2half[i]);
            if(difference[i] > max_error[i]){
                take_step = false;
                printf("max_error exceeded! difference[%d] = %f\n",i,difference[i]);
                break;
            }
        }
        if(take_step){
            for(int i=0;i<nVar;i++){
                x[n+1][i] = x_1full[i];
            }
            h=h*2;  // double step size
        }
        else{
            h=h/2;  // halve step size
            continue;   // skip n++ -> redo with new h
        }
        printf("last n: %d \n", n);
        n += 1;   // increase step counter
    }
}


