
#include "mex.h"
#include "IMG.h"
void mexFunction( int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray*prhs[] )
{
   if (nrhs!=2) mexErrMsgTxt("Only the IMG structure and # of its expected");
   //if (nlhs!=5) mexErrMsgTxt("A returned IMG structure isn't present");
   if (mxGetNumberOfElements(prhs[0])!=1) mexErrMsgTxt("One structure expected.");

//srand (0);
//srand48(7);
   IMG sp_img;
   sp_img.mxReadIMG(prhs[0]);

   int its = getInput<double>(prhs[1]);
   for (int i=0; i<its; i++)
      if (sp_img.move_switch_IMG())
         break;

   sp_img.mxWriteIMG(plhs,prhs[0]);
}
