/*

INPUT:
   x,y,z  positions given in [kpc]
   n_sub  number f particles in array
   mass   given in solar masses
   softening [kpc]

OUTPUT:
   U_P    potential energy array

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void potential_nbody_ccc(const float *x, const float *y, const float *z, const int n_sub, const float softening, float *U_p){

  int i,j;
  float U_temp=0;

  float soft2 = softening*softening;

  //--- Loop over all particles
  for (i=0; i<n_sub; i++) {

    U_temp = 0.0; 
    //--- Loop again over all particles
    for (j=0;j<n_sub;j++){
      
      //--- Avoid itself
      if (i == j) continue;

      U_temp += 1.0 / sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]) + (z[j] - z[i])*(z[j] - z[i]) + soft2);

    } // end for j

    //--- Store particle's potential energy. 
    //    Multiply constants outside summatory
    U_p[i] = U_temp;

  } // end for i


}

