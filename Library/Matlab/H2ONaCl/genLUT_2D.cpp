#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include "H2ONaCl.H"

double add(double b, double c)
{
  std::cout<<"test mex: "<<b+c<<std::endl;
  return b+c;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check arguments
    if(nrhs<7 || nrhs >9){
    mexErrMsgTxt("Usage: createLUT_2D(xmin, xmax, ymin, ymax, const_which_var, TorH, min_level=4, max_level=6)");
    }
    double xmin = *(mxGetPr(prhs[0]));
    double xmax = *(mxGetPr(prhs[1]));
    double ymin = *(mxGetPr(prhs[2]));
    double ymax = *(mxGetPr(prhs[3]));
    double constZ = *(mxGetPr(prhs[4]));
    LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var = (LOOKUPTABLE_FOREST::CONST_WHICH_VAR)(mxGetScalar(prhs[5]));
    LOOKUPTABLE_FOREST::EOS_ENERGY TorH = (LOOKUPTABLE_FOREST::EOS_ENERGY)(mxGetScalar(prhs[6]));

    H2ONaCl::cH2ONaCl sw;
  #if USE_OMP == 1
    sw.set_num_threads(8);
  #endif
    sw.createLUT_2D_TPX(xmin, xmax, ymin, ymax, constZ, const_which_var, TorH);
    sw.save_lut_to_binary("lut.bin");
    sw.save_lut_to_vtk("lut.vtu");
    // double *a;

    // double b, c;

    // plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    // plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    // plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    // a = mxGetPr(plhs[0]);

    // b = *(mxGetPr(prhs[0]));

    // c = *(mxGetPr(prhs[1]));

    // *mxGetPr(plhs[0])=add(b,c);
    // *mxGetPr(plhs[1]) = b;
    // *mxGetPr(plhs[2]) = c;
}
