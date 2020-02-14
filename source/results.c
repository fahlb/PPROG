#include <math.h>           // compile with -lm
#include <stdio.h>          // IO: fprintf(), etc.
#include "ODEs.h"


int line2array(FILE* myFile,int N, double myArr[]){
    int rt; // return value of fscanf()
    for(int i=0;i<N;i++){
        rt=fscanf(myFile, "%lf", &myArr[i]);
        if(rt==EOF){printf("");break;}
        if(rt!=1  ){printf("");break;}
    }
    return rt;
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
    fprintf(output, "#t[s]            ");
    fprintf(output, "r_cm_0[m]       ");
    fprintf(output, "r_cm_1[m]       ");
    fprintf(output, "r_cm_2[m]       ");
    fprintf(output, "\n");

    // check if file open
    if(input == NULL){
        printf("unable to open %s \n", inPath);
        return;
    }
    int nVar = 2*N*D;
    double dArr[nVar+1];    // array for nVar variables + time

    int first_loop = 1;
    int initial_order;

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
            var[i-1] = dArr[i];
        }

        // write calculated distances:

        fprintf(output, " %+.6e", t);
        for(int i=0;i<N;i++){
            fprintf(output, "   %+.6e", distance2cm(i,N,D,m,var));
        }
        // do some stupid workaround because shit is fucked somehow..
        // for some odd reason distance2cm(i,..)-distance2cm(j,..) returns nonsense
//        double d_cm_1 = distance2cm(1,N,D,m,var);
//        double d_cm_2 = distance2cm(2,N,D,m,var);
        double d_cm_1 = r_polar_i(1,N,D,var);
        double d_cm_2 = r_polar_i(2,N,D,var);
        
        double difference = d_cm_1 - d_cm_2;
        fprintf(output, "   %+.6e", difference);
    
        fprintf(output,"\n");

    }
    fclose(input);
    fclose(output);
}

int main(){
    int N = 3;
    int D = 2;

    double m[N];// masses of objects
    m[0] = 5.688e26;    // M_Saturn     [kg]
    m[1] = 1.98e18;     // M_Janus      [kg]
    m[2] = 5.5e17;      // M_Epimetheus [kg]

    write_dist2cm("../output/data.dat", "../output/cm_data.dat",3,2,m);

    return 0;
}