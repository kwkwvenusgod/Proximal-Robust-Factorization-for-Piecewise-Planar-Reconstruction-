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

#define NUMARGS 5
#define NUMOUT 4

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, MATRIX3 ,  "double"); // data
   checkInput(prhs, 1, CELL); // sp_indices
   checkInput(prhs, 2, VECTOR  ,  "double"); // var
   checkInput(prhs, 3, CELL); // neighbors
   checkInput(prhs, 4, SCALAR  ,  "logical"); // useFlow

   const int *dims = mxGetDimensions(prhs[0]);
   int N = dims[0]; // number of data points
   int T = dims[1]; // number of dimensions
   int D = dims[2]; // number of frames
   int NT = N*T;
   arr(double) data = getArrayInput<double>(prhs[0]);
   arr(double) var = getArrayInput<double>(prhs[2]);
   arr(arr(int)) neighbors = allocate_memory< arr(int) >(N*T);
   arr(int) Nneighbors = allocate_memory<int>(N*T);
   for (int i=0; i<NT; i++)
   {
      if (mxGetCell(prhs[3],i))
      {
         neighbors[i] = getArrayInput<int>(mxGetCell(prhs[3],i));
         Nneighbors[i] = mxGetNumberOfElements(mxGetCell(prhs[3],i));
      }
      else
         Nneighbors[i] = 0;
   }

   bool useFlow = getInput<bool>(prhs[4]);

   plhs[0] = mxCreateNumericMatrix(N, N, mxDOUBLE_CLASS, mxREAL);
   plhs[1] = mxCreateNumericMatrix(N, N, mxDOUBLE_CLASS, mxREAL);
   plhs[2] = mxCreateNumericMatrix(N, 1, mxINT32_CLASS, mxREAL);
   plhs[3] = mxCreateNumericMatrix(N, 1, mxINT32_CLASS, mxREAL);
   arr(double) dist = getArrayInput<double>(plhs[0]);
   arr(double) ncount = getArrayInput<double>(plhs[1]);
   arr(int) alive_born = getArrayInput<int>(plhs[2]);
   arr(int) alive_dead = getArrayInput<int>(plhs[3]);

   for (int i=0; i<N; i++)
   {
      alive_born[i] = T;
      alive_dead[i] = T-1;
      bool found = false;
      for (int t=0; t<T; t++)
      {
         bool alive = !mxIsInf(data[i+t*N]) && !mxIsNaN(data[i+t*N]) && mxGetCell(prhs[1],i+t*N) && mxGetNumberOfElements(mxGetCell(prhs[1],i+t*N))>0;
         if (!found && alive)
         {
            found = true;
            alive_born[i] = t;
         }
         else if (found && !alive)
         {
            alive_dead[i] = t-1;
            break;
         }
      }
   }

   double inf = mxGetInf();
   for (int i=0; i<N*N; i++)
   {
      dist[i] = inf;
      ncount[i] = 0;
   }

   for (int i=0; i<N; i++)
   {
      for (int j=i+1; j<N; j++)
      {
         double total = 0;
         int tstart = max(alive_born[i], alive_born[j]);
         int tstop = min(alive_dead[i], alive_dead[j]);

         if ((tstop >= tstart && !useFlow) || (tstop>tstart && useFlow))
         {
            double neighborCount = 0;
            for (int t=tstart; t<=tstop; t++)
            {
               bool tempFound = false;
               for (int n=0; n<Nneighbors[i+t*N]; n++)
               {
                  //if (i==0 && j==261)
                     //mexPrintf("%d / %d neighbor=%d\n", n, Nneighbors[i+t*N], neighbors[i+t*N][n]);
                  if (neighbors[i+t*N][n]==j)
                  {
                     neighborCount++;
                     tempFound = true;
                     break;
                  }
               }
            }
            ncount[i+j*N] = neighborCount / (tstop-tstart+1);
            ncount[j+i*N] = neighborCount / (tstop-tstart+1);
            ncount[i+j*N] = neighborCount;// / (tstop-tstart+1);
            ncount[j+i*N] = neighborCount;// / (tstop-tstart+1);

            if (true)
            {
               double time;
               double maxDist = 0;
               double total2 = 0;

               if (useFlow) // flow only defined for first 1:time-1 frames
               {
                  time = tstop-tstart;
                  for (int t=tstart; t<tstop; t++)
                  {
                     double temp = 0;
                     for (int d=0; d<D-2; d++)
                        temp += pow(data[i+t*N+d*NT]-data[j+t*N+d*NT],2)/(var[d]);
                     total += temp;
                     maxDist = max(maxDist, temp);
                  }
                  for (int d=D-2; d<D; d++)
                  {
                     double tempi = 0;
                     double tempj = 0;
                     int tflowstart = max(tstop-5,tstart);
                     for (int t=tflowstart; t<tstop; t++)
                     {
                        tempi += data[i+t*N+d*NT];
                        tempj += data[j+t*N+d*NT];
                     }
                     if (tflowstart!=tstop)
                     {
                        tempi = tempi / (tstop-tflowstart);
                        tempj = tempj / (tstop-tflowstart);
                     }
                     total += pow(tempi-tempj,2)/var[d];
                  }
               }
               else
               {
                  time = tstop-tstart+1;
                  for (int t=tstart; t<=tstop; t++)
                  {
                     double temp = 0;
                     for (int d=0; d<D-2; d++)
                        temp += pow(data[i+t*N+d*NT]-data[j+t*N+d*NT],2)/(var[d]);
                     total += temp;
                     maxDist = max(maxDist, temp);
                  }
               }

   //            if (exp(-total / time)==1)
   //               mexPrintf("%d - %d\n", tstart, tstop);

               dist[i+j*N] = total / time;
               //dist[i+j*N] = maxDist;
               dist[j+i*N] = dist[i+j*N];
            }
         }
      }
   }

   deallocate_memory(neighbors);
}
