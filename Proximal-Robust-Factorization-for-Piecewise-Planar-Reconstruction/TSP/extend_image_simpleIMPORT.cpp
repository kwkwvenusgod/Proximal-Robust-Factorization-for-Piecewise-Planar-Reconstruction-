#include "mex.h"
#include "helperMEX.h"
#include "extend_image2.cpp"

#define NUMARGS 3
#define NUMOUT 1
// [extended, sdf] = extend_image2IMPORT(image, mask)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
      mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
      mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, MATRIX,  "double"); //image
   checkInput(prhs, 1, MATRIX,  "logical"); //mask
   checkInput(prhs, 2, SCALAR,  "double"); //band_size

   //Create matrix for the return argument
   const int *dims = mxGetDimensions(prhs[0]);
   const int numDims = mxGetNumberOfDimensions(prhs[0]);
   int xdim = dims[0];
   int ydim = dims[1];
   int D = 1;
   if (numDims>2)
      D = dims[2];
   int N = xdim*ydim;
   
   plhs[0] = mxCreateNumericArray(numDims, dims, mxDOUBLE_CLASS, mxREAL);

   //Assign pointers to input and output
   arr(double) im = getArrayInput<double>(prhs[0]);
   arr(bool) mask = getArrayInput<bool>(prhs[1]);
   double band_size = getInput<double>(prhs[2]);

   arr(double) extended = getArrayInput<double>(plhs[0]);

   double farthest_size;
   arr(double) sdf = allocate_memory<double>(N);
   arr(double) closest_curve_x = allocate_memory<double>(N);
   arr(double) closest_curve_y = allocate_memory<double>(N);
   arr(double) farthest_x = allocate_memory<double>(N);
   arr(double) farthest_y = allocate_memory<double>(N);
   calc_closest_curve(xdim, ydim, mask, band_size, sdf, closest_curve_x, closest_curve_y, farthest_size, farthest_x, farthest_y);
   
   for (int i=0; i<N; i++)
   {
      int index = closest_curve_x[i] + closest_curve_y[i]*xdim;
      for (int d=0; d<D; d++)
         extended[i + d*N] = im[index + d*N];
   }

   deallocate_memory(closest_curve_x);
   deallocate_memory(closest_curve_y);
   deallocate_memory(farthest_x);
   deallocate_memory(farthest_y);
   deallocate_memory(sdf);
}
