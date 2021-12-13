#include "mex.h"
#include <matrix.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;
#include "H2ONaCl.H"

#define MAX(a, b) ((a)>(b) ? (a):(b))
#define MIN(a, b) ((a)<(b) ? (a):(b))

void mGetVector(const mxArray *prhs, double **out, const char *varname, mwSize *out_m)
{
  mwSize m, n;
  double *temp;

  m = mxGetM(prhs);
  n = mxGetN(prhs);

  if(!(m==1 || n==1)){
    char buff[256];
    sprintf(buff, "'%s' must be a vector.\n", varname);
    mexErrMsgTxt(buff);
  }

  if(!mxIsDouble(prhs)){    char buff[256];
    sprintf(buff, "'%s' must be of type 'double'.\n", varname);
    mexErrMsgTxt(buff);
  }

  temp = (double*)mxGetData(prhs);
  *out = temp;
  *out_m = MAX(m,n);
}

void mGetMatrix(const mxArray *prhs, double **out, const char *varname, mwSize *out_m, mwSize *out_n)
{
  mwSize m, n;
  double *temp;

  m = mxGetM(prhs);
  n = mxGetN(prhs);

  if(!mxIsDouble(prhs)){
    char buff[256];
    sprintf(buff, "'%s' must be of type 'double'.\n", varname);
    mexErrMsgTxt(buff);
  }

  temp   = (double*)mxGetData(prhs);
  *out   = temp;
  *out_m = m;
  *out_n = n;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check arguments
    if(nrhs!=3){
    mexErrMsgTxt("Usage: lookupLUT_2D(filename, x, y)\nNote that the x, y order MUST BE consistent with the lookup table file.");
    }
    /* allocate input variables */
    double *xMat, *yMat;
    mwSize  m, n;
    /* allocate output variables */
    double *Rho_out, *T_out;
    LOOKUPTABLE_FOREST::NeedRefine *needRefine_out;

    string filename = mxArrayToString(prhs[0]);
    mGetMatrix(prhs[1], &xMat, "px", &m, &n);
    mGetMatrix(prhs[2], &yMat, "py", &m, &n);

    plhs[0] = mxCreateNumericMatrix(m, n, mxINT32_CLASS, mxREAL);
    needRefine_out = (LOOKUPTABLE_FOREST::NeedRefine*)mxGetData(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    Rho_out = (double*)mxGetData(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);
    T_out = (double*)mxGetData(plhs[2]);
    // maybe make a safety check, input matrix dimension consistency
    clock_t start = clock();
    STATUS("Passing data end, start loading data ...");
    H2ONaCl::cH2ONaCl sw;
    H2ONaCl::PROP_H2ONaCl prop_lookup;
    sw.loadLUT(filename); //\Todo pass this sw pointer as input arg from outside, don't always load the LUT file.
    STATUS_time("Load LUT end. ", clock()-start);
    start = clock();
    STATUS("Start looking up ...");
    const int dim = 2;
    LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf = NULL;
    int index = 0;
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        index = i*n +j;
        targetLeaf = sw.lookup_only(prop_lookup, xMat[index], yMat[index]); 
        needRefine_out[index] = targetLeaf->user_data->need_refine;
        T_out[index] = prop_lookup.T;
        Rho_out[index] = prop_lookup.Rho;
      }
    }
    STATUS_time("Lookup end. ", clock()-start);
}
