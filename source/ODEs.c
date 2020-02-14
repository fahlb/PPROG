#include "ODEs.h"
#include "read_parameters.h"

const double G = 6.67384e-11;   // [N m^2 2^-2]

double r_ij_d(int i, int j, int d, int D, const double var[]){
    // returns d component of of vector from j to i
    // d=0:x-component, d=1:y-component
    return var[2*i*D +d] - var[2*j*D +d];   // 
}

double distance_ij(int i, int j, int D, const double var[]){
    // returns absolute value of distance between objects i and j

    // get vector components of r_ij:
    double r_ij[D];
    for(int d=0;d<D;d++){
        r_ij[d] = r_ij_d(i,j,d,D,var);
    }
    // calculate square of distance vector r_ij:
    double r_squared = 0;
    for(int d=0;d<D;d++){
        r_squared += r_ij[d]*r_ij[d];
    }
    return sqrt(r_squared);
}

double F_ij(int i, int j, int D, double m[], const double var[]){
    // returns directionsless value of force between bodies
    //const double G = 6.67384e-11;    // [N m^2 2^-2]
    return -G*m[i]*m[j]/pow(distance_ij(i,j,D,var),2);
}

double T_i(int i, int D, double m[], const double var[]){
    // returns kinetic energy of body i

    // calculate square of abolute value of velocity:
    double v_squared = 0;
    for(int d=0;d<D;d++){
        v_squared += pow(var[(2*i+1)*D + d],2);
    }
    return 0.5 * m[i] * v_squared;
}

double V_ij(int i, int j, int D, double m[], const double var[]){
    // return potential between bodies i and j
    //const double G = 6.67384e-11;    // [N m^2 2^-2]
    return -G*m[i]*m[j]/distance_ij(i,j,D,var);
}

double E_tot(int N, int D, double m[], const double var[]){
    // returns total energy; could be put in ODEs() for k=-2 !
    double V = 0;
    double T = 0;

    for(int i=0;i<N;i++){
        T += T_i(i,D,m,var);

        // calculate V_i:
        for(int j=0;j<N;j++){
            if(i==j || i<j){continue;}  // skip i=j and duplicates
            V += V_ij(i,j,D,m,var);
        }
    }
    //printf("E_tot: V = %e \n",V);
    //printf("E_tot: T = %e \n",T);
    return (V + T);
}

double distance2cm(int i, int N, int D, double m[], const double var[]){
    double distance2cm_i;

    // find coordinatess of center of mass (cm):
    double cm_coord[D];

    // get total mass of system
    double m_tot = 0;   // total mass of system
    for(int k=0;k<N;k++){
        m_tot += m[k];
    }
    for(int d=0;d<D;d++){
        cm_coord[d] = 0;
        for(int k=N;k<N;k++){
            cm_coord[d] += m[k] / m_tot * var[2*k*D +d];
        }
    }
    //printf("distance2cm: i=%d\n", i);
    //printf("distance2cm: x_cm = %+.3e \n", cm_coord[0]);
    //printf("distance2cm: y_cm = %+.3e \n", cm_coord[1]);
    //printf("\n");


    // find distances of objects to cm:
    double dist_square;  // temp square of distance to cm
    // get components of distance:
    for(int d=0;d<D;d++){
        dist_square += pow(var[2*i*D +d] - cm_coord[d], 2);
    }
    distance2cm_i = sqrt(dist_square);

    return distance2cm_i;
}

double r_polar_i(int i, int N, int D, const double var[]){
    double temp = 0;
    for(int d=0;d<D;d++){
        temp += pow(var[2*i*D +d],2);
    }
    return sqrt(temp);
}


double ODEs(int k, double t, const double var[]){
    /* x_nd, v_nd
     *
     * Indix nd:
     * n = 0,..,N-1 with N=3 number of bodies (index for body)
     *  n=0: Saturn
     *  n=1: Janus
     *  n=2: Epimetheus
     * 
     * d = 0,..,D-1 with D=2 number of dimensions (index of dimension)
     *  d=0: x coordinate
     *  d=0: y coordinate
     * 
     * Table of variables coresponding to index i
     *  var[0] : x_01(t)
     *  var[1] : x_02(t)
     *  var[2] : v_01(t)
     *  var[3] : v_02(t)
     *  var[4] : x_11(t)
     *  var[5] : x_12(t)
     *  var[6] : v_11(t)
     *  var[7] : v_12(t)
     *  var[8] : x_21(t)
     *  var[9] : x_22(t)
     *  var[10]: v_21(t)
     *  var[11]: v_22(t)
     */

//    const int N = 3; // do in main() ?
    const int N = 3; // testing with only one moon
    const int D = 2; // ..?

    // read from input file:
    double parameters[2*N*D+N+5];
    read_parameters("input.dat", N, D, parameters);
    // masses
    double m[N];    // do also in main() ? + read from file?
    m[0] = parameters[2*N*D +0];    // M_Saturn     [kg]
    m[1] = parameters[2*N*D +1];    // M_Janus      [kg]
    m[2] = parameters[2*N*D +2];    // M_Epimetheus [kg]
    double t_f = parameters[2*N*D +N+4];

    // calculate every F_i (dimensions seperately) (basically F_nd !!!)
    double F_i[N][D];  // d komponent of force acting on body i (=n, F_nd)
    for(int d=0;d<D;d++){
        for(int i=0;i<N;i++){   // n=i
            F_i[i][d] = 0;
            for(int j=0;j<N;j++){   // n=j
                if(i==j){continue;}
                // sum up forces
                F_i[i][d] += F_ij(i,j,D,m,var)
                             *(r_ij_d(i,j,d,D,var)/distance_ij(i,j,D,var));
            }             // ^ multiply with unit vector in dimension d !
        }
    }

    //printf("i = %d", k);

    // dx_nd/dt = v_nd
    // dv_nd/dt = F_nd/m_n

    switch(k){
        case 0: // dx_00/dt = v_01(t)
            return var[2];

        case 1: // dx_01/dt = v_02(t)
            return var[3];

        case 2: // dv_00/dt = F_00/m_0
            return F_i[0][0]/m[0];

        case 3: // dv_01/dt = F_01/m_0
            return F_i[0][1]/m[0];

        case 4: // dx_10/dt = v_11(t)
            return var[6];

        case 5: // dx_11/dt = v_12(t)
            return var[7];

        case 6: // dv_10/dt = F_10/m_1
            return F_i[1][0]/m[1];

        case 7: // dv_11/dt = F_11/m_1
            return F_i[1][1]/m[1];

        case 8: // dx_20/dt = v_21(t)
            return var[10];

        case 9: // dx_21/dt = v_22(t)
            return var[11];

        case 10:// dv_20/dt = F_20/m_2
            return F_i[2][0]/m[2];

        case 11:// dv_21/dt = F_21/m_2
            return F_i[2][1]/m[2];

        case -1:    // return 0, if conditio to stop integration is reached
                    // return 1, otherwise (continue integration)
            {
            // time to stop at [s] <- should not be hardcoded!
            if(t< t_f){return 1;}
            if(t>=t_f){return 0;}
            }

        case -2:    // returns total energy of system
            return E_tot(N,D,m,var);

        default: 
            printf("invalid argument i=%d for ODEs()\n", k);
            return 0;
    }
}
