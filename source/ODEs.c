#include <math.h>
#include <stdio.h>

/* TODO:
    - fix F_ij()
    - do just Saturn + Janus first? 
*/

double r_ij_d(int i, int j, int d, int D, const double var[]){
    // returns d component of of vector from j to i
    // d=0:x-component, d=1:y-component
    return var[2*i*D +d] - var[2*j*D +d];   // 
}

double distance_ij(int i, int j, int D, const double var[]){
    // returns absolute value of distance between objects i and j
    double r_ij[D];
    for(int d=0;d<D;d++){
        r_ij[d] = r_ij_d(i,j,d,D,var);
    }

    double result_squared = 0;
    for(int d=0;d<D;d++){
        result_squared += r_ij[d]*r_ij[d];
    }
    return sqrt(result_squared);
}

double F_ij(int i, int j, int D, double m[], const double var[]){
    // returns absolute value of force between objects with correct sign!!
    const double G = 6.67384e-11;    // [N m^2 2^-2]
    return -G*m[i]*m[j]/pow(distance_ij(i,j,D,var),2);
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

    // masses
    double m[N];    // do also in main() ? + read from file?
    m[0] = 5.688e26;    // M_Saturn     [kg]
    m[1] = 1.98e18;     // M_Janus      [kg]
    m[2] = 5.5e17;      // M_Epimetheus [kg]

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
            double t_f = 315360000; // 10yrs
            if(t< t_f){return 1;}
            if(t>=t_f){return 0;}
            }
    
        default: 
            printf("invalid argument i=%d for ODEs()\n", k);
            return 0;
    }
}
