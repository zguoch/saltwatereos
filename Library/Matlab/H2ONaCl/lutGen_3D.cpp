#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include "H2ONaCl.H"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check arguments
    if(nrhs<8 || nrhs >11){
    mexErrMsgTxt("Usage: createLUT_2D(xmin, xmax, ymin, ymax, zmin, zmax, TorH, outputFile_without_ext, min_level=4, max_level=6, num_threads=1)");
    }

    int min_level = 4, max_level = 6, num_threads = 1;
    double xmin = *(mxGetPr(prhs[0]));
    double xmax = *(mxGetPr(prhs[1]));
    double ymin = *(mxGetPr(prhs[2]));
    double ymax = *(mxGetPr(prhs[3]));
    double zmin = *(mxGetPr(prhs[4]));
    double zmax = *(mxGetPr(prhs[5]));
    LOOKUPTABLE_FOREST::EOS_ENERGY TorH = (LOOKUPTABLE_FOREST::EOS_ENERGY)(mxGetScalar(prhs[6]));
    string filename = mxArrayToString(prhs[7]);
    if(nrhs>=9)min_level = (int)(*(mxGetPr(prhs[8])));
    if(nrhs>=10)max_level = (int)(*(mxGetPr(prhs[9])));
    if(nrhs==11)num_threads = (int)(*(mxGetPr(prhs[10])));

    H2ONaCl::cH2ONaCl sw;
  #if USE_OMP == 1
    sw.set_num_threads(num_threads);
  #endif
    sw.createLUT_3D(xmin, xmax, ymin, ymax, zmin, zmax, TorH, min_level, max_level, Update_prop_rho | Update_prop_h | Update_prop_T);
    sw.save_lut_to_binary(filename+".bin");
    sw.save_lut_to_vtk(filename+".vtu");
}
