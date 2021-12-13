#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include "H2ONaCl.H"

// Print information of a binary file of LUT.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check arguments
    if(nrhs!=1){
    mexErrMsgTxt("Usage: infoLUT(filename)");
    }
    string filename = mxArrayToString(prhs[0]);
    H2ONaCl::cH2ONaCl sw;
    sw.loadLUT(filename);
}
