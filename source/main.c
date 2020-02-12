// testing environment for Runge Kutta

#include <math.h>           // compile with -lm
#include <stdio.h>          // IO: fprintf(), etc.
#include "Runge_Kutta_4.c"  // RK4()
#include "ODEs.c"

/*
int line2array(FILE* myFile,int N, double myArr[]){
    int rt; // return value of fscanf()
    for(int i=0;i<N;i++){
        rt=fscanf(myFile, "%lf", &myArr[i]);
        if(rt==EOF){break;}
        if(rt!=1  ){break;}
    }
    return rt;
}


void readfile(FILE* file, int N){
    printf("readfile: \n");
    
    while(1){
        int rt; // retrun value for fscanf
        double dArray[N];
        fscanf(file, "%*[^\n]\n");  // skip lines without values (?) 
                                    // (https://stackoverflow.com/a/16108311)
        rt=line2array(file, N, dArray);
        // if rt is not 1 as expected:
        if(rt==EOF){break;}
        if(rt!=1  ){break;}

        for(int i=0;i<N;i++){
            printf("%+.3e", dArray[i]);
            printf("\n");
        }
    }
}

void write_dist2cm(char* inPath, char* outPath, int N, int D, double m[]){
    // loads saved date at given path and calculates distance 
    // of moons to center of mass for every (saved) time step
    // N is the number of floats in one line

    FILE* input;
    input = fopen(inPath, "r");

    FILE* output;
    output = fopen(outPath, "w+");

    // write head for output file
    fprintf(output, "# t[s]         ");
    fprintf(output, "r_cm_0[m]    ");
    fprintf(output, "r_cm_1[m]    ");
    fprintf(output, "r_cm_2[m]    ");
    fprintf(output, "\n");

    // check if file open
    if(input == NULL){
        printf("unable to open %s \n", inPath);
        return;
    }
    int nVar = 2*N*D;
    double dArr[nVar+1];    // array for nVar variables + time

    while(1){
        int rt; // retrun value for fscanf
        fscanf(input, "%*[^\n]\n");  // skip lines without values (?) 
                                    // (https://stackoverflow.com/a/16108311)
        rt=line2array(input, nVar+1, dArr);
        // if rt is not 1 as expected:
        if(rt==EOF){break;}
        if(rt!=1  ){break;}

        // separate time from variables (familar form)
        double t = dArr[0];
        double var[nVar];

        for(int i=1;i<nVar+1;i++){
            var[i] = dArr[i];
        }
        
        // find distance of each object to center of mass:

        // find coordinatess of center of mass (cm):
        double cm_coord[D];
        double x_cm;
        double y_cm;
        // get total mass of system
        double m_tot = 0;   // total mass of system
        for(int i=0;i<N;i++){
            m_tot += m[i];
        }
        for(int i=0;i<N;i++){
            for(int d=0;d<D;d++){
                cm_coord[d] = m[i] / m_tot * var[2*i*D +d];
            }
        }

        // find distances of objects to cm:
        double distance2cm[N];  // distance of objects i to center of mass (cm)

        double dist_square[N];  // temp square of distance to cm
        // get components of distance:
        for(int i=0;i<N;i++){
            for(int d=0;d<D;d++){
                dist_square[i] += var[2*i*D +d]*cm_coord[d];
            }
            distance2cm[i] = sqrt(dist_square[i]);
        }

        // write calculated distances:

        fprintf(output, " %+.3e", t);
        for(int i=0;i<N;i++){
            fprintf(output, "   %+.3e", distance2cm[i]);
        }
        fprintf(output,"\n");
    }
    fclose(input);
    fclose(output);
}
*/

int main(){
    int N = 3;  // number of objects
    int D = 2;  // number of dimensions
    
    double m[N];// masses of objects
    m[0] = 5.688e26;    // M_Saturn     [kg]
    m[1] = 1.98e18;     // M_Janus      [kg]
    m[2] = 5.5e17;      // M_Epimetheus [kg]

    // set stepsize h, start and finish time
    const double h = 1e2;         // [s]
    const double h_min = 7*(24*60*60); // minimum step size to be written
    const double precission = 1e-6;  // precission when using adaptive step size
    // ^ needs to be sufficiently small; wont work otherwise (no idea why..)
    const double t_0 = 0.0;       // [s]
    const double t_f = 10*365.25*60*60;      // [s]

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
    AdaptiveRK4(&ODEs, h, &t,nVar, var, &output, h_min,precission); 

    // don't forget to close file!
    fclose(output);

    // now lets find the time 

    return 0;
}
