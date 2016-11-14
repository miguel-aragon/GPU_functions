/*

Get closest distances between two sets of particles. Use mark array to
flag both sets.

NOTES:
  - Pass a mask array with one number per species of particles. Only two species
  - Particle positions MUST be different between sets because of a dirty fix to
    avoid self-distances.

   Miguel Aragon Calvo  Oct/2010
   

"This software contains source code provided by NVIDIA Corporation."

"Glue c code based on galaxy collision demo"

History:

   - 21/10/2010 Code based on the GPU potential code

*/

/*
 * Copyright 1993-2006 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.  This source code is a "commercial item" as
 * that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer software" and "commercial computer software
 * documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 */


#include <math.h>
#include <stdio.h>

#define BLOCKDIM 256

__constant__ float softeningSquared;

// Macros to simplify shared memory addressing
#define SX(i) sharedPos[i+blockDim.x*threadIdx.y]

//==================================
//  Find the distance between two particles in different sets
//==================================
__device__ float bodyBodyInteraction(float minDis, float4 pos_j, float4 pos_i) {


  //--- Only do pairs between both sets
  if (pos_i.w == pos_j.w) return minDis;

  //--- Distance vector
  float3 r;
  r.x = pos_i.x - pos_j.x;
  r.y = pos_i.y - pos_j.y;
  r.z = pos_i.z - pos_j.z;
  
  //--- Squared distance
  float distSqr = r.x*r.x + r.y*r.y + r.z*r.z;
  
  //--- Avoid itself
  if (distSqr != 0) {
    
    float invDist = distSqr;
    minDis += invDist;
  }
  
  return minDis;
}


//==================================
// This is the "tile_calculation" function from the GPUG3 article.
//==================================
__device__ float tile_potential(float4 myPos, float pot) {

    extern __shared__ float4 sharedPos[];    

    unsigned long i = 0;

    //--- Here we unroll the loop: LOOP_UNROLL = 4
    for (unsigned int counter = 0; counter < blockDim.x; ) {
        pot = bodyBodyInteraction(pot, SX(i++), myPos); 
        pot = bodyBodyInteraction(pot, SX(i++), myPos); 
        pot = bodyBodyInteraction(pot, SX(i++), myPos); 
        pot = bodyBodyInteraction(pot, SX(i++), myPos); 
	counter += 4;
    }

    return pot;
}

//==================================
// WRAP is used to force each block to start working on a different 
// chunk (and wrap around back to the beginning of the array) so that
// not all multiprocessors try to read the same memory locations at once.
//==================================
#define WRAP(x,m) (((x)<m)?(x):(x-m))  // Mod without divide, works on values from 0 up to 2m

__device__ float computePotential(float4 bodyPos, float4* positions, int numBodies){

    extern __shared__ float4 sharedPos[];

    float pot = 0.0f;
    
    int p = blockDim.x;
    int q = blockDim.y;
    int n = numBodies;
    int numTiles = n / (p * q);

    for (int tile = blockIdx.y; tile < numTiles + blockIdx.y; tile++) {
		
      sharedPos[threadIdx.x+blockDim.x*threadIdx.y] = positions[WRAP(blockIdx.x + tile, gridDim.x) * p + threadIdx.x];
       
      __syncthreads();

      //--- This is the "tile_calculation" function from the GPUG3 article.
      pot = tile_potential(bodyPos, pot);
        
      __syncthreads();
    }

    return pot;
}

__global__ void integrateBodies(float4* Pos, float* poten, int numBodies){

  //--- Get the index of this thread ?
  int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
  
  float4 pos_i = Pos[index];
  
  //--- Return potential
  float pot = computePotential(pos_i, Pos, numBodies);

  //--- Put potential in fourth field (mass)
  //Pos[index].w = pot;

  poten[index] = pot;
	
}

//============================================================
//  
//============================================================



//---  Memory buffers on the GPU
float4 *pos1;
float  *pote;
int old_buf;
int np;
int np_nopad;


//==============================================
//---  Interface routines...
//==============================================

#include <unistd.h>

extern "C"
{

  #define BLOCKSIZE 256

  //==================================
  /* Allocate GPU memory and set initial positions */
  //==================================
  void init_nbody(int _n, float *_pos, float *_mass){
    int i;

    /* Pad with zero mass particles if not a multiple of blocksize */
    np= (_n/BLOCKSIZE)*BLOCKSIZE;
    if(np<_n)
      np = np + BLOCKSIZE;

    /* Allocate GPU arrays */
    cudaMalloc((void **) &pos1, sizeof(float4)*np);
    cudaMalloc((void **) &pote, sizeof(float) *np);

    /* Prepare initial conditions */
    float *posbuf = (float *) malloc(4*sizeof(float)*np);
    float *potbuf = (float *) malloc(  sizeof(float)*np);
    for(i=0; i<_n; i++){
		posbuf[4*i+0] = _pos[3*i+0];
		posbuf[4*i+1] = _pos[3*i+1];
		posbuf[4*i+2] = _pos[3*i+2];
		posbuf[4*i+3] = _mass[i];
		potbuf[i]     = 0.0f;
      }

    //--- Pad particles
    for(i=_n; i<np; i++){
		posbuf[4*i+0] = 0.0;
		posbuf[4*i+1] = 0.0;
		posbuf[4*i+2] = 0.0;
		posbuf[4*i+3] = 0.0;
		potbuf[i]     = 0.0;
		//--- Initialize to a very large value for the bubble minimum finder
		pote[i]       = 1000000000.0;
      }

    /* Copy to GPU */
    old_buf = 1;
    cudaMemcpy(pos1, posbuf, sizeof(float4)*np, cudaMemcpyHostToDevice);  
    cudaMemcpy(pote, potbuf, sizeof(float )*np, cudaMemcpyHostToDevice);  
    free(posbuf);
    free(potbuf);

    np_nopad = _n;

  }

  cudaEvent_t evt;
  int underway = 0;

  //==================================
  /* Do the actual potential  */
  //==================================
  void compute_potential(void) {

    /* Execute the kernel */
    dim3 Dg(np/BLOCKSIZE);
    dim3 Db(BLOCKSIZE);
    size_t Ns = 4 * sizeof(float) * BLOCKSIZE;

    cudaEventCreate(&evt);

    integrateBodies <<< Dg, Db, Ns >>> (pos1, pote, np);

    cudaEventRecord(evt, 0);

    underway = 1;
  }

  //==================================
  /* Check whether the calculation is done */
  //==================================
  int nbody_finished() {

    if(cudaEventQuery(evt) == cudaErrorNotReady){
	return 0;
      } else {
	cudaEventDestroy(evt);
	underway = 0;
	return 1;
      }
  }

  //==================================
  /* Shut down and deallocate */
  //==================================
  void dealloc_nbody(int dump){
    if(underway==1)
      while(nbody_finished()==0);
    cudaFree(pos1);
    cudaFree(pote);
  }

  //==================================
  /*  Set softening */
  //==================================
  void set_softening(float eps)
  {
    float eps2 = eps*eps;
    cudaMemcpyToSymbol("softeningSquared", &eps2, sizeof(float), 0, cudaMemcpyHostToDevice);
  }

  //==================================
  //--- Retrieve positions
  //==================================
  void get_positions(float *buf){

    if(underway==1)
      while(nbody_finished()==0);

    float *pos = (float *) malloc(4*sizeof(float)*np);
    cudaMemcpy(pos, pos1, sizeof(float)*4*np, cudaMemcpyDeviceToHost);
   
    int i;
    for(i=0;i<np_nopad;i++){
		buf[4*i+0] = pos[4*i+0];
		buf[4*i+1] = pos[4*i+1];
		buf[4*i+2] = pos[4*i+2];
		buf[4*i+3] = pos[4*i+3];
      }
    free(pos);
  }

  //==================================
  //--- Retrieve potential
  //==================================
   void get_potential(float *buf){

    //--- Wait until computation is finish
    if(underway==1)
      while(nbody_finished()==0);
    
    float *pot = (float *) malloc(sizeof(float)*np);
    cudaMemcpy(pot, pote, sizeof(float)*np, cudaMemcpyDeviceToHost);
    
    int i;
    for(i=0;i<np_nopad;i++)
      {
	buf[i] = pot[i];
      }
    free(pot);
  }
  
} //--- end extern "C"
