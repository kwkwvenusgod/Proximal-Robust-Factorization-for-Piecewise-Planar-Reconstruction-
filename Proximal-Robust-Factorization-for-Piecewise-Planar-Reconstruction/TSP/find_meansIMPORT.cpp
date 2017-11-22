// =============================================================================
// == calc_sdfIMPORT.cpp
// == --------------------------------------------------------------------------
// == The MEX interface file to calculate a signed distance function
// == --------------------------------------------------------------------------
// == Copyright 2011. MIT. All Rights Reserved.
// == Written by Jason Chang 06-13-2011
// =============================================================================

#include "mex.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "debugMEX.h"
#include "helperMEX.h"
#include "matrix.h"


#define NUMARGS 3
#define NUMOUT 1



// find_valid_noiseIMPORT(phi, mask, cc, logpim_p, logpim_m, Eright_same, Eright_diff, Edown_same, Edown_diff, FG_topologies, BG_topologies, allowInsert, allowRemove)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs[0], MATRIX, "int32"); // z
   checkInput(prhs[1], SCALAR, "double"); // Nz
   checkInput(prhs[2], MATRIX, "double"); // im

   arr(int) z = getArrayInput<int>(prhs[0]);
   int Nz = getInput<double>(prhs[1]);
   arr(double) im = getArrayInput<double>(prhs[2]);

   const int *dims = mxGetDimensions(prhs[2]);
   int xdim = dims[0];
   int ydim = dims[1];
   int N = xdim*ydim;

   plhs[0] = mxCreateNumericMatrix(Nz, 1, mxDOUBLE_CLASS, mxREAL);
   arr(double) output = getArrayInput<double>(plhs[0]);

   arr(double) count = allocate_memory<double>(Nz);
   memset(count, 0, sizeof(double)*Nz);
   memset(output, 0, sizeof(double)*Nz);

   for (int i=0; i<N; i++)
   {
      if (z[i]>=0)
      {
         output[z[i]] += im[i];
         count[z[i]]++;
      }
   }

   for (int i=0; i<Nz; i++)
      if (count[i]>0)
         output[i] /= count[i];

   deallocate_memory(count);
}
