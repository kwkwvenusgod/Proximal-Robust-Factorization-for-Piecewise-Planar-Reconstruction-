
#include "mex.h"
#include "IMG.h"
void mexFunction( int nlhs, mxArray *plhs[],
                                  int nrhs, const mxArray*prhs[] )
{
             srand (0);
             srand48(7);
             
             IMG sp_img;
             sp_img.mxReadIMG(prhs[0]);
             int its = getInput<double>(prhs[1]);
             arr(double) thres = getArrayInput<double>(prhs[2]);
             //sp_img.Test_gating(its,thres); 
             sp_img.Match_BIP(thres); 
             sp_img.mxWriteIMG(plhs,prhs[0]);
            }
