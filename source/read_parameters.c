#include "read_parameters.h"

void read_parameters(char* myPath, int N, int D, double parameters[2*N*D+N+5]){
    
    FILE* input;
    //sprintf(myPath, "input_alt.txt");
    input = fopen(myPath, "r");

    if(input == NULL){
        printf("unable to open %s \n", myPath);
        return;
    }

    fscanf(input, "%*[^\n]\n");  // skip lines without values (?) 
    // read initial values
    for(int i=0;i<2*N*D;i++){
        fscanf(input, "%lf" , &parameters[i]);
    }
    fscanf(input, "%*[^\n]\n");  // skip lines without values (?) 
    fscanf(input, "%*[^\n]\n");  // skip lines without values (?) 
    // read masses
    for(int i=2*D*N;i<2*D*N+N;i++){
        fscanf(input, "%lf", &parameters[i]);
    }

    fscanf(input, "%*[^\n]\n");  // skip lines without values (?) 
    fscanf(input, "%*[^\n]\n");  // skip lines without values (?) 
    // read integration parameters
    fscanf(input, "%lf", &parameters[2*D*N+N +0]);
    fscanf(input, "%lf", &parameters[2*D*N+N +1]);
    fscanf(input, "%lf", &parameters[2*D*N+N +2]);
    fscanf(input, "%lf", &parameters[2*D*N+N +3]);
    fscanf(input, "%lf", &parameters[2*D*N+N +4]);

/*
    printf("\n\n");
    // test reading  
    for(int i=0;i<2*N*D+N+5;i++){
        printf("%f \n", parameters[i]);
    }
*/
    //read_parameters("input_alt.txt", N,D,parameters);
    fclose(input);

}
