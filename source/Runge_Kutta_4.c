#include <math.h>
#include <stdio.h>
#include <stdbool.h>
//#include "Runge_Kutta_4.h"

//double ODEs(int i, double t, const double var[]);

/* TODO:
    - EU is broken; why the fuck is it more precise than RK?!
    - ???
*/

void NextEU(double (*f)(int, double, const double*),  // function pointer
            double h, double *t, int nVar, double x[]){
    
    for(int i=0;i<nVar;i++){
        x[i] += h*(*f)(i,*t,x);
    }
}

void EU(double (*f)(int, double, const double[]),  // function pointer
        double h, double *t, int nVar, double x[nVar], 
        FILE** output, double h_min){
    
    // calculate steps
    int n = 0;
    while((*f)(-1,*t,x)){
        NextEU(f,h,t,nVar,x);
        *t += h;
        
        printf("last n: %d \n", n);
        n += 1;   // increase step counter

        // write step
        if(h<h_min && (n%(int)(h_min/h))){continue;} //skip if not multiple of n_skip
        fprintf(*output, "  %.3e", *t);
        for(int i=0; i<nVar;i++){
            fprintf(*output, "    %.3e", x[i]);
        }
        fprintf(*output, "\n");
    }
}

void NextRK4(double (*f)(int, double, const double*),  // function pointer
             double h, double *t, int nVar, double x[]){

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
        x[i] += h/6*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
    }
//    *t += h;


}

void RK4(double (*f)(int, double, const double[]),  // function pointer
        double h, double *t, int nVar, double x[nVar], 
        FILE** output, double h_min){

    int n = 0;  // start with first step (n=0 : initial conditions)
//    while(n<nSteps-1){
    while((*f)(-1,*t,x)){    // breaks code :(   
        
        NextRK4(f,h,t,nVar,x);
        *t += h;

        // increase step counter
        printf("last n: %d \n", n);
        n += 1;
        // write step
        if(h<h_min && (n%(int)(h_min/h))){continue;} //skip if not multiple of n_skip
        fprintf(*output, "  %.3e", *t);
        for(int i=0; i<nVar;i++){
            fprintf(*output, "    %.3e", x[i]);
        }
        fprintf(*output, "\n");
    }
}


void AdaptiveRK4(double (*f)(int, double, const double[]),  // function pointer
                double h, double *t, int nVar, double x[nVar], 
                FILE** output, double h_min){

    double max_error[nVar]; // make argument; maximum error for each variable
    for(int i=0;i<nVar;i++){
        max_error[i] = 1e-5;
    }   // temporary!!!

    int n = 0;  // start with first step (n=0 : initial conditions)
    double x_1full[nVar];   // temporary x[i] for 1 full step
    double x_2half[nVar];   // temporary x[i] for 2 half steps
    
//    while(n<nSteps-1){
    while((*f)(-1,*t,x)){    // breaks code?! :(   
        // save copy of x[i] for both advances
        for(int i=0;i<nVar;i++){
            x_1full[i] = x[i];
            x_2half[i] = x[i];
        }
        // make one full step:
        NextRK4(f,h,t,nVar,x_1full);
        // make two half steps:
        NextRK4(f,h/2,t,nVar,x_2half);
        double newt = *t+h/2;
        NextRK4(f,h/2,&(newt),nVar,x_2half);
        
        double difference[nVar];    // difference of 1 full and 2 half steps
        bool take_step = true;
        for(int i=0;i<nVar;i++){
            difference[i] = fabs(x_1full[i]-x_2half[i]);
            if(difference[i] > max_error[i]){
                take_step = false;
            //  printf("max_error exceeded! difference[%d] = %f\n",i,difference[i]);
                break;
            }
        }
        if(take_step){
            for(int i=0;i<nVar;i++){
                x[i] = x_1full[i];
            }
            *t+=h;
            h=h*2;  // double step size
        }
        else{
            h=h/2;  // halve step size
            continue;   // skip n++ -> redo with new h
        }

        printf("last n: %d \n", n);
        n += 1;   // increase step counter

        // write step
        if(h<h_min && (n%(int)(h_min/h))){continue;} //skip if not multiple of n_skip
        fprintf(*output, "  %.3e", *t);
        for(int i=0; i<nVar;i++){
            fprintf(*output, "    %.3e", x[i]);
        }
        fprintf(*output, "\n");
    }
}


