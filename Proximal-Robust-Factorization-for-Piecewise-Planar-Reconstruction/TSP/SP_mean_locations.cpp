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
#include "linkedList.cpp"
#include "utils.h"

#define NUMARGS 2
#define NUMOUT 1



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, SCALAR ,  "double"); //K
   checkInput(prhs, 1, MATRIX ,  "uint32"); //label

   int K = getInput<double>(prhs[0]);
   arr(unsigned int) label = getArrayInput<unsigned int>(prhs[1]);
   int X = mxGetM(prhs[1]);
   int Y = mxGetN(prhs[1]);

   plhs[0] = mxCreateNumericMatrix(K,2,mxDOUBLE_CLASS,mxREAL);
   arr(double) total_X = getArrayInput<double>(plhs[0]);
   arr(double) total_Y = total_X + K;
   arr(double) count = allocate_memory<double>(K,0);

   for (int x=0; x<X; x++) for (int y=0; y<Y; y++)
   {
      int index = x+y*X;
      unsigned int k = label[index];
      total_X[k] += x;
      total_Y[k] += y;
      count[k]++;
   }

   for (int k=0; k<K; k++)
   {
      if (count[k]>0)
      {
         total_X[k] /= count[k];
         total_Y[k] /= count[k];
      }
      else
      {
         total_X[k] = -1;
         total_Y[k] = -1;
      }
   }

   deallocate_memory(count);
}
